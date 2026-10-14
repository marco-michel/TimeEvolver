/**
* Round trip tests for the HDF5 output. saveHDF5 and loadHDF5 are each other's
* inverse, so they can be checked against one another without reference files.
*/

#include <complex>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "matrixDataTypes.h"
#include "krylovExceptions.h"

static int failures = 0;

static void check(const std::string& label, bool condition)
{
	if (condition)
		std::cout << "ok: " << label << std::endl;
	else
	{
		std::cerr << "FAIL: " << label << std::endl;
		failures++;
	}
}

#ifdef USE_HDF

/**
* A non square matrix, so that a confusion of the two dimensions cannot pass
* unnoticed, with entries that are neither purely real nor symmetric.
*/
static TE::smatrix buildMatrix()
{
	std::vector<std::complex<double>> values;
	std::vector<size_t> columns, rows;
	const size_t numColumns = 7, numRows = 5;

	for (size_t row = 0; row != numRows; row++)
	{
		values.push_back(std::complex<double>(1.0 + row, -0.5 * row));
		columns.push_back((2 * row) % numColumns);
		rows.push_back(row);
	}

	TE::smatrix mat(values.data(), columns.data(), rows.data(), values.size(),
		(unsigned int)numColumns, (unsigned int)numRows);
	mat.hermitian = true;
	mat.upperTri = true;
	return mat;
}

static void roundTrip()
{
	const std::string file = "testHDF5RoundTrip.h5";
	TE::smatrix original = buildMatrix();
	original.saveHDF5(file);

	TE::smatrix restored;
	restored.loadHDF5(file);

	check("column dimension survives", restored.n == original.n);
	check("row dimension survives", restored.m == original.m);
	check("number of entries survives", restored.numValues == original.numValues);
	check("hermitian flag survives", restored.hermitian == original.hermitian);
	check("upperTri flag survives", restored.upperTri == original.upperTri);
	check("sym flag survives", restored.sym == original.sym);

	bool entriesMatch = true;
	for (size_t i = 0; i != original.numValues; i++)
	{
		if (std::abs(restored.values[i] - original.values[i]) > 1e-15
			|| restored.columns[i] != original.columns[i]
			|| restored.rowIndex[i] != original.rowIndex[i])
			entriesMatch = false;
	}
	check("all entries survive unchanged", entriesMatch);

	std::remove(file.c_str());
}

/**
* A matrix large enough to cross the compression threshold, which exercises the
* chunked path rather than the contiguous one.
*/
static void roundTripCompressed()
{
	const std::string file = "testHDF5RoundTripLarge.h5";
	const size_t dim = 400000;

	std::vector<std::complex<double>> values(dim);
	std::vector<size_t> columns(dim), rows(dim);
	for (size_t i = 0; i != dim; i++)
	{
		values[i] = std::complex<double>(std::sin(0.001 * i), 0.0);
		columns[i] = i;
		rows[i] = i;
	}

	TE::smatrix original(values.data(), columns.data(), rows.data(), dim,
		(unsigned int)dim, (unsigned int)dim);
	original.saveHDF5(file);

	TE::smatrix restored;
	restored.loadHDF5(file);

	bool entriesMatch = (restored.numValues == dim);
	for (size_t i = 0; entriesMatch && i != dim; i++)
	{
		if (std::abs(restored.values[i] - original.values[i]) > 1e-15
			|| restored.columns[i] != original.columns[i])
			entriesMatch = false;
	}
	check("compressed matrix survives unchanged", entriesMatch);

	std::remove(file.c_str());
}

static void reportsMissingFile()
{
	TE::smatrix mat;
	try
	{
		mat.loadHDF5("this-file-does-not-exist.h5");
		check("missing file throws", false);
	}
	catch (const TE::krylovIOError&)
	{
		check("missing file throws krylovIOError", true);
	}
	catch (const std::exception&)
	{
		//An HDF5 error that was not translated would land here, or not be caught at all.
		check("missing file throws krylovIOError rather than another type", false);
	}
}

#endif

int main()
{
#ifdef USE_HDF
	roundTrip();
	roundTripCompressed();
	reportsMissingFile();
#else
	std::cout << "built without HDF5, nothing to check" << std::endl;
#endif

	if (failures != 0)
	{
		std::cerr << failures << " HDF5 check(s) failed" << std::endl;
		return 1;
	}
	std::cout << "All HDF5 round trip checks passed." << std::endl;
	return 0;
}
