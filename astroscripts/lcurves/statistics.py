

from abc import ABCMeta, abstractmethod
import numpy as np
from math import factorial as fact
from scipy.stats import poisson, gamma, rv_continuous
from scipy.optimize import fsolve, fmin
from .lcurves import _bkgration_msg

def rate_uplimit(rawcnt: int, bkgcnt_r: float, exp: float, prob: float =0.9) -> float:
    """Calculate upper limit to count rate.

    Calculate quantile of the count rate corresponding to the given probability.
    Poisson's distribution is assumed.

    :param rawcnt: Number of raw counts accumulated in the source aperture.
    :param bkgcnt_r: Number of background counts rescaled(!) to the source aperture.
    :param exp: Exposure time.
    :param prob: Probability. The default is 0.9 (90% upper limit).
    :returns: Rate upper limit, float.
    """
    return (fsolve(lambda x: poisson.cdf(rawcnt, bkgcnt_r+x)
           - ( 1 -prob), rawcnt - bkgcnt_r)[0] ) /exp


class _Loredo(rv_continuous, metaclass=ABCMeta):
    """Abstract base class for count rate statistical distributions."""

    @abstractmethod
    def __init__(self, **kwargs):
        self._rawcnt = self._exp = None
        parargs = {}   # Arguments to pass to the parent constructor
        for key in ('b', 'xtol'):
            if key in kwargs:
                parargs[key] = kwargs.pop(key)

        for argname, argval in kwargs.items():
            setattr(self, '_' + str(argname), argval)

        super().__init__(a=0, **parargs)

    def _pdf(self, srate: float, *args):
        """Probability density function of the source count rate.

        :param srate: source count rate
        """
        vectorized_pdf = np.vectorize(self._pdf_scalar)
        return vectorized_pdf(srate)

    @abstractmethod
    def _pdf_scalar(self, srate: float):
        """The function to be vectorized in _pdf"""

    def interval_mode(self, alpha=0.683) -> (float, float):
        """Return credible interval symmetric probability around the mode.

        :param alpha: Confidence level
        :returns: a pair of values, tuple
        """
        xmod = self.mode()
        Mex = self.expect()
        cdf_at_xmod = self.cdf(xmod)
        # How close is the mode to the median?
        if cdf_at_xmod < alpha / 2:  # There is no left tail
            xmin = 0.0
            xmax = fsolve(lambda x: self.cdf(x) - alpha, 1.5 * Mex)[0]
        else:  # Get tails of equal probabilities
            xmin = fsolve(lambda x: cdf_at_xmod - self.cdf(x) - alpha/2, Mex/2)[0]
            xmax = fsolve(lambda x: self.cdf(x) - cdf_at_xmod - alpha/2, 1.5*Mex)[0]
        assert xmin <= xmod <= xmax
        return xmin, xmax

    def mode(self):
        """Return the distribution mode."""
        return fmin(lambda x: -1 * self.pdf(x), self.expect(), disp=False)[0]

    def _interval_shortest_get_B(self, A, alpha, Me) -> float:
        """Return B for a particular A.

        Int(pdf(x), x=[A, B]) = alpha
        or
        cdf(B) - cdf(A) = alpha
        """
        cdfA=self.cdf(A)
        B = fsolve(lambda x: self.cdf(x) - cdfA - alpha, Me)[0]
        return B

    def interval_shortest(self, alpha: float) -> (float, float):
        """Return the shortest credible interval."""
        Me = self.expect()
        A = fmin(lambda x: self._interval_shortest_get_B(x, alpha, Me) - x, 0.0, disp=False,
                 xtol=Me / 10000)[0]
        B = self._interval_shortest_get_B(A, alpha, Me)
        assert A <= self.mode() <= B
        return A, B


class Loredo_rate_zero_background(_Loredo):
    """
    Posterior probability distribution for the count rate in the case
    then the background is known to be zero

    :param rawcnt: Number of accumulated counts
    :param exp: Exposure time
    """

    def __init__(self, rawcnt: int, exp: float):
        super().__init__(rawcnt=rawcnt, exp=exp)

    def _pdf_scalar(self, srate):
        """Probability density function of the source count rate.

        :param srate: source count rate
        """
        T = self._exp
        n = self._rawcnt
        # Gamma distribution multiplied by exposure time
        return T * gamma.pdf(T * srate, n + 1, scale=1)  # Loredo (5.5)


class Loredo_rate_known_background(_Loredo):
    """
    Posterior probability distribution for count rate in the case of
    the background is known a priori (not measured empirically)

    :param rawcnt: Number of total counts accumulated in the source aperture
    :param exp: Exposure time
    :param brate: Known background count rate
    """

    def __init__(self, rawcnt: int, exp: float, *, brate: float):
        self._brate = None
        super().__init__(rawcnt=rawcnt, exp=exp, brate=brate)

    def _pdf_scalar(self, srate: float, *args):
        """Probability density function of the source count rate.

        :param srate: source count rate
        """
        T = self._exp
        n = self._rawcnt
        brate = self._brate
        C = 1 / poisson.cdf(n, brate * T)  # Normalization, Loredo (5.7)
        # Gamma distribution normalized and multiplied by exposure time
        return C * T * gamma.pdf(
            T * (srate + brate), n + 1, scale=1)  # Loredo (5.6)


class Loredo_rate_measured_background(_Loredo):
    """
    Posterior probability distribution for count rate in the case
    when the background

    :param rawcnt: Number of total counts accumulated in the source aperture
    :param exp: Exposure time
    :param bkgcnt: Number of counts accumulated in the background aperture
    :param bkgratio: Ratio of apertures  S(bkg)/S(src)
    :param b: rv_continues's b parameter (uplimit to replace infinity)
    """

    def __init__(self, rawcnt: int, exp: float , *, bkgcnt: int, bkgratio: float,
                 b: float = None):
        self._bkgcnt = self._bkgratio = None
        _bkgration_msg(bkgratio)
        b = b or 100*(rawcnt+bkgcnt+1)/exp  # Upper limit for the probability density function argument
        super().__init__(rawcnt=rawcnt, exp=exp, bkgcnt=bkgcnt,
                         bkgratio=bkgratio, b=b, xtol=1e-14)

        Non, Noff, S, T = rawcnt, bkgcnt, bkgratio, exp
        term = lambda i: (1 + S) ** i * fact(Non + Noff - i) / fact(Non - i)
        denominator = np.sum([term(j) for j in range(Non + 1)])
        self._Ci = np.zeros(shape=(Non+1, 1))
        for i in range(Non + 1):
            self._Ci[i] = term(i) / denominator  # Normalization, Loredo (5.14)

    def _pdf_scalar(self, srate: float):
        T = self._exp
        Ci = self._Ci

        return sum([Ci[i ] * T *(srate *T )**i * np.exp(- 1 *srate *T ) /fact(i)
                    for i in range(len(Ci))])  # Loredo (5.13)


