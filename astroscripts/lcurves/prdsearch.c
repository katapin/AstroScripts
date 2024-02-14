#include <stdio.h>
#include <math.h>

/*Calculate raw (non-normalized) Lomb-Scargle power at
arbitrary frequencies. See Horne & Baliunas, 1986, 
%ApJ, 302, 757
http://adsabs.harvard.edu/abs/1986ApJ...302..757H
The arguments are:
  time  - time stamps of the light curve;
  data  - count rate or flux;
  freq  - desired frequencies;
  power - an array to write the result, must have
          the same length as freq;
  inlen - length of the time and data arrays 
  oulen - length of the freq and power arrays
  prog  - it true, prints progress

*/
int lspower(float *time, float *data, float *freq, 
            float *power, int inlen, int oulen)
{
    int i, j;
    float sins, coss;  //sum of sin and cos in eq. (2) 
    float omtau;       // omega*tau in eq(2)
    float arg;
    float ysin, sinsq, ycos, cossq;

    for (i=0; i<oulen; i++)   // Frequencies
    {
        sins = coss = omtau = 0.0f;
        for (j=0; j<inlen; j++)   // Calculate phase shift
        {
            sins += sinf(4.0f*M_PI*freq[i]*time[j]);
            coss += cosf(4.0f*M_PI*freq[i]*time[j]);
        }
        omtau = 0.5f*atan2f(sins, coss);
        arg = 2.0*M_PI*freq[i]*time[j] -omtau;

        ycos = cossq = ysin = sinsq = 0.0f;
        for (j=0; j<inlen; j++)
        {
            ycos += data[j]*cosf(arg);
            ysin += data[j]*sinf(arg);
            cossq = powf(cosf(arg), 2);
            sinsq = powf(sinf(arg), 2);
        }
        power[i] = powf(ycos, 2)/cossq +
            powf(ysin, 2)/sinsq;
        printf("%f %f\n", arg, power[i]);
    }
}