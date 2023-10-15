#!/usr/bin/python
import os, sys
from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn
import myfits
from myfits import FilePathAbs
from myfits.plot import PlotVector
import numpy as np
from enum import IntEnum


class DataVector(MaskedColumn):  # PLot vector

	_vectortype = IntEnum('VectorType', ['WithoutErrors', 'SymmetricErrors'])

	TypeWithoutErrors = _vectortype.WithoutErrors
	TypeSymmetricErrors = _vectortype.SymmetricErrors

	def __new__(cls, values: Column, *args, **kwargs):
		return super().__new__(cls, values)

	def __init__(self, values: Column, *, symmetric_errors=None):
		self._type = self.TypeWithoutErrors
		if symmetric_errors is not None:
			if type(symmetric_errors) != MaskedColumn:
				self._errors = symmetric_errors
				self._type = self.TypeSymmetricErrors

	def isgood(self):
		return not np.all(self.data.mask)  # All values masked, it's an empty array

	@property
	def errors(self):
		if self._type == self.TypeSymmetricErrors:
			return self._errors
		else:
			np.zeros_like(self)

	@property
	def type(self):
		return self._type


# def __getitem__(self,):
# pass


class TableWithVectors(Table):  # QDP table
	def __init__(self, table: Table):
		super().__init__(table)
		vectors = []
		colnames = iter(self.colnames)
		for col in colnames:
			if col + '_err' in self.colnames:
				curvector = DataVector(self[col], symmetric_errors=self[next(colnames)])
			elif col + '_perr' in self.colnames:
				raise NotImplementedError()
			else:
				curvector = DataVector(self[col])
			if curvector.isgood() == False:
				continue
			vectors.append(curvector)
		self.vectors = tuple(vectors)

	@property
	def width(self):
		return len(self.colnames)

	@property
	def real_width(self):
		"""Rutern number of columns with taking into accont error columns and empty vectors."""
		vec_width = {DataVector.TypeWithoutErrors: 1, DataVector.TypeSymmetricErrors: 2}
		return sum(vec_width[vec.type] for vec in self.vectors)


class QDPfile():
	"""Class to store the content of a QPD file."""

	def __init__(self, filepath: FilePathAbs):
		self.tables = tuple(self._read(filepath))
		self.vectors = tuple(vec for tbl in self.tables for vec in tbl.vectors)

	@property
	def width(self):
		return self.tables[0].width

	@staticmethod
	def _read(filepath: FilePathAbs):
		"""Return list of TableWithVectors object.

		It's the internal read functions needed because astropy.io.ascii.read()
		doesn't understand string starting from '@'. So we intend to 
		remove such a string before feeding ascii.read()"""
		with open(filepath) as f:
			lines = f.readlines()
			if len(lines) < 2:
				print(f'Wrong file {path}')
				sys.exit(1)

		# Remove string not supported by astropy's ascii.read()
		if lines[1][0] == '@':
			del lines[1]
		tables = []
		i = 0
		try:
			while True:  # Read all available tables
				tables.append(TableWithVectors(ascii.read(lines, table_id=i, guess=False, format='qdp')))
				i += 1
		except IndexError:
			pass
		return tables

	def __getitem__(self, index):
		if type(index) == int:
			return self.tables[index]
		elif type(index) == tuple:
			return self.tables[index[0]].vectors[index[1]]
		else:
			raise TypeError('Unsupported type of index')


class SingleObsData():
	"""Class to store the (datapoints + smooth model) pair for one observation."""

	def __init__(self, data_file: FilePathAbs, model_file: FilePathAbs):
		self.dataqdp = self._load_and_check(data_file, 'data')
		self.modelqdp = self._load_and_check(model_file, 'model')

		self.datX = self.dataqdp[(0, 0)]  # Data, X values
		self.datY = self.dataqdp[(0, 1)]  # Data, Y values
		self.datResid = self.dataqdp[(1, 1)]  # D ata, Residulas
		self.datStM = self.dataqdp[(0, 2)]
		modX = self.modelqdp[(0, 0)]

		mask = np.logical_and(modX > self.datX[0] / 2, modX < 2 * self.datX[-1])
		# self.modX = modX[mask]
		# self.modY = self.modelqdp[(0,1)][mask]
		self.modX = self.modelqdp[(0, 0)]
		self.modY = self.modelqdp[(0, 1)]
		self.modComps = []
		for vec in self.modelqdp[0].vectors[2:]:
			# self.modComps.append(vec[mask])
			self.modComps.append(vec)

	def _load_and_check(self, filepath: FilePathAbs, filetype):
		qdp = QDPfile(filepath)
		if len(qdp.tables) > (2 if filetype == 'data' else 1):
			print(f"Error: file '{filepath}' contains too many tables for type {filetype}")
			raise ValueError(f"Wrong content of '{filepath}'")
		if filetype == 'data':
			if len(qdp.tables) > 1:
				if len(qdp.tables[0]) != len(qdp.tables[0]):
					print(f"Error: Spectrum and uncertainties in '{filepath}' have different length")
					raise ValueError(f"Wrong content of '{filepath}'")
		return qdp

	@property
	def nmodelcomp(self):
		"""Return number of independent model components."""
		return len(self.modelqdp.vectors) - 2  # Exclude X-axis and the model itself

	@property
	def residuals_present(self):
		return True if len(self.dataqdp.tables) > 1 else False

	def determine_required_ncols(self, with_residuals):
		ndatcols = 7 if with_residuals else 5  # X Xerr Y Yerr [Uncert Uncert_err ] StepModel
		nmodcols = self.nmodelcomp + 4 + (1 if with_residuals else 0)  # X Xerr Mod NO Comp1 Comp2  ...
		return max(ndatcols, nmodcols)  # New number of columns

	# @staticmethod
	# def _generate_table(ncols, *args):
	# lines=[]
	# curcols=len(args)
	# for vals in zip(*args):
	# lines.append(' '.join(str(x) for x in vals) + (' {}'.format(' '.join(['NO']*(ncols-curcols))) if ncols != curcols else '') )
	# return lines

	@staticmethod
	def _generate_table(*args):
		lines = []
		for vals in zip(*args):
			lines.append(' '.join(str(x) for x in vals))
		return lines

	# def reformat(self, ncols, with_residuals, with_row_delimiter=False):
	# lines=[]
	# if with_residuals:
	# resid_vectors = (self.datResid, self.datResid.errors)

	# lines.extend(self._generate_table(ncols, self.datX, self.datX.errors, self.datY, self.datY.errors, *resid_vectors))
	# lines.append(' '.join(['NO']*ncols))
	# lines.extend(self._generate_table(ncols, self.modX, self.modX.errors, self.modY, *self.modComps))
	# if with_row_delimiter:
	# lines.append(' '.join(['NO']*ncols))
	# return lines

	def reformat(self, ncols, with_residuals, with_row_delimiter=False):
		lines = []
		esid_vectors = []
		if with_residuals:
			resid_vectors = [self.datResid, self.datResid.errors]

		dummy_dat = ['No'] * len(self.datX)  # Dummy column for data table
		args = [self.datX, self.datX.errors, self.datY, self.datY.errors, *resid_vectors, self.datStM]
		args.extend([dummy_dat] * (ncols - len(args)))  # Append dummy columns
		lines.extend(self._generate_table(*args))
		lines.append(' '.join(['NO'] * ncols))

		dummy_mod = ['No'] * len(self.modX)  # Dummy column for model table
		args = [self.modX, self.modX.errors, self.modY, dummy_mod]
		modComps = self.modComps.copy()  # Copy to do pop()
		if len(modComps) > 0:
			args.append(modComps.pop(0))
			if with_residuals:
				args.append(dummy_mod)
			args.extend(modComps)
		assert ncols >= len(args)
		args.extend([dummy_mod] * (ncols - len(args)))
		lines.extend(self._generate_table(*args))
		if with_row_delimiter:
			lines.append(' '.join(['NO'] * ncols))
		return lines

	def generate_header(self, pcofile, with_residuals):
		first_line = 'READ SERR 1 2'
		if with_residuals: first_line += ' 3'
		return [first_line, '@' + pcofile, '!']


def _main():
	import argparse
	parser = argparse.ArgumentParser(
		description="Combine QPD files from XSPEC "
					"to product more pretty spectra")
	parser.add_argument('-o', '--output', nargs=1, required=True,
						help="Name of the QDP-file to save the result")
	parser.add_argument('-p', '--pcofile', nargs=1,
						help="Name of the pco-file", default='spectra.pco')
	parser.add_argument('--without_residuals', help="Don't include residuals",
						action='store_true')
	parser.add_argument('files', nargs='+',
						help="Pairs of QDP-files with the observed data "
							 "and the smooth model")

	argnspace = parser.parse_args(sys.argv[1:])

	if len(argnspace.files) % 2 != 0:
		myfits.die('Number of the input qdp-files must be even.')

	infiles = argnspace.files
	items = []
	for dat_file, mod_file in [(infiles[i], infiles[i + 1]) for i in range(0, len(infiles), 2)]:
		items.append(SingleObsData(dat_file, mod_file))

	if argnspace.without_residuals:
		with_residuals = False
	else:
		resid_check_list = [item.residuals_present for item in items]
		if any(resid_check_list) != all(resid_check_list):
			print("Error: Some files has residuals, but some doesn't. To disable residuals for" \
				  "all the files use option '--without_residuals'")
			sys.exit(1)
		with_residuals = all(resid_check_list)  # Use automatic value
	ncols = max([item.determine_required_ncols(with_residuals) for item in items])
	lines = items[0].generate_header(argnspace.pcofile, with_residuals)
	for item in items[:-1]:
		lines.extend(item.reformat(ncols, with_residuals, with_row_delimiter=True))
	lines.extend(items[-1].reformat(ncols, with_residuals, with_row_delimiter=False))
	with open(argnspace.output[0], 'w') as f:
		f.write('\n'.join(lines))


if __name__ == '__main__':
	_main()
