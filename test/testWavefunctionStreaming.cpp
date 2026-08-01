/**
* Checks the streaming of the sampled wavefunction into the result file.
*
* The history of states is never held in memory, so the file is the only place
* it exists and these checks read it back.
*/

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "krylovTimeEvolver.h"
#include "krylovObservables.h"
#include "krylovExceptions.h"
#include "matrixDataTypes.h"
#include "krylovHelper.h"

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

#include "hdf5Support.h"

/**
* A hopping chain, which spreads the initial state over the whole basis and so
* gives every component of every sampled state something to be wrong about.
*/
static std::unique_ptr<TE::smatrix> hoppingChain(size_t dim)
{
	std::vector<std::complex<double>> values;
	std::vector<size_t> columns, rows;
	for (size_t i = 0; i + 1 < dim; i++)
	{
		values.push_back(1.0); columns.push_back(i + 1); rows.push_back(i);
		values.push_back(1.0); columns.push_back(i);     rows.push_back(i + 1);
	}
	return std::unique_ptr<TE::smatrix>(new TE::smatrix(values.data(), columns.data(), rows.data(),
		values.size(), (unsigned int)dim, (unsigned int)dim));
}

/**
* An observable that stops the time evolution once it has been sampled a given
* number of times, so that the early exit path can be exercised.
*/
class stoppingObservable : public krylovBasicObservable
{
public:
	stoppingObservable(size_t stopAfter) : krylovBasicObservable("stopper"), stopAfter(stopAfter), seen(0)
	{
		type = VOID_TYPE_OBS;
	}

	std::complex<double> expectation(std::complex<double>* vec, int len) override
	{
		(void)vec; (void)len;
		if (++seen > stopAfter)
			throw requestStopException();
		if (expectationValues != nullptr && sampleIndex < numSamples)
			expectationValues[sampleIndex++] = 0.0;
		return 0.0;
	}

private:
	size_t stopAfter;
	size_t seen;
};

/**
* Read the wavefunction dataset back into memory. Only ever called on the small
* systems used here, where doing so is affordable.
*/
static std::vector<std::complex<double>> readStates(const std::string& file, hsize_t& numRows, hsize_t& numColumns)
{
	H5::H5File handle(file, H5F_ACC_RDONLY);
	H5::DataSet dataset = handle.openDataSet("wavefunction");

	hsize_t dimensions[2];
	dataset.getSpace().getSimpleExtentDims(dimensions);
	numRows = dimensions[0];
	numColumns = dimensions[1];

	std::vector<TE::hdf5::complexType> buffer(numRows * numColumns);
	if (!buffer.empty())
		dataset.read(buffer.data(), TE::hdf5::complexDataType());

	std::vector<std::complex<double>> states(buffer.size());
	for (size_t i = 0; i != buffer.size(); i++)
		states[i] = std::complex<double>(buffer[i].r, buffer[i].i);
	return states;
}

/**
* Every sampling point ends up in the file, in order, and the last row is the
* state the time evolution reports as its result.
*/
static void streamsEverySample()
{
	const size_t dim = 32;
	const double maxT = 1.0, samplingStep = 0.1;
	const size_t expectedSamples = (size_t)std::floor(maxT / samplingStep) + 1;

	std::vector<std::complex<double>> initial(dim, 0.0);
	initial[dim / 2] = 1.0;

	std::vector<std::unique_ptr<krylovBasicObservable>> observables;
	krylovTimeEvolver evolver(maxT, initial.data(), samplingStep, std::move(observables), hoppingChain(dim),
		1.0, 1e-8, 20, false, false);

	parameter_list parameters;
	parameters.push_back(paraPush("chainLength", false, (int)dim));
	//A caller passing the sampling step in as a model parameter as well, as the
	//example does
	parameters.push_back(paraPush("samplingStep", false, samplingStep));

	hdf5ResultWriter writer(parameters, "testWavefunctionStreaming");
	const std::string file = writer.fileName();
	evolver.setSampleWriter(&writer);

	std::unique_ptr<krylovReturn> result(evolver.timeEvolve());
	writer.writeMetadata(*result);

	hsize_t numRows = 0, numColumns = 0;
	std::vector<std::complex<double>> states = readStates(file, numRows, numColumns);

	check("one row per sampling point", numRows == (hsize_t)expectedSamples);
	check("row length is the Hilbert space dimension", numColumns == (hsize_t)dim);
	check("the number of rows agrees with the reported sample count", numRows == (hsize_t)result->numSamples);

	bool initialMatches = (numRows != 0);
	for (hsize_t i = 0; initialMatches && i != numColumns; i++)
		initialMatches = std::abs(states[i] - initial[i]) < 1e-14;
	check("the first row is the initial state", initialMatches);

	//The last sampling point is the state the evolution ends on, so the two have
	//to agree; this is what ties the streamed rows to the computed result.
	bool finalMatches = (numRows != 0);
	for (hsize_t i = 0; finalMatches && i != numColumns; i++)
		finalMatches = std::abs(states[(numRows - 1) * numColumns + i] - result->evolvedState[i]) < 1e-12;
	check("the last row is the evolved state", finalMatches);

	//The same state repeated in row after row would pass every check above, so
	//make sure the rows really are different points in time.
	double largestChange = 0.0;
	for (hsize_t i = 0; numRows > 1 && i != numColumns; i++)
		largestChange = std::max(largestChange, std::abs(states[numColumns + i] - states[i]));
	check("consecutive rows differ", largestChange > 1e-6);

	//Rows are useless without a time axis
	H5::H5File handle(file, H5F_ACC_RDONLY);
	double storedStep = 0.0, storedTotal = 0.0, storedEvolved = 0.0;
	TE::hdf5::readAttribute(handle, "samplingStep", H5::PredType::NATIVE_DOUBLE, storedStep);
	TE::hdf5::readAttribute(handle, "totalTime", H5::PredType::NATIVE_DOUBLE, storedTotal);
	TE::hdf5::readAttribute(handle, "evolvedTime", H5::PredType::NATIVE_DOUBLE, storedEvolved);
	handle.close();

	check("the sampling step is recorded", std::abs(storedStep - samplingStep) < 1e-15);
	check("the requested evolution time is recorded", std::abs(storedTotal - maxT) < 1e-15);
	check("the time reached is recorded", std::abs(storedEvolved - maxT) < 1e-9);

	check("the last sample lands on the end of the interval",
		std::abs((double)(numRows - 1) * storedStep - maxT) < 1e-9);

	std::remove(file.c_str());
}

/**
* A time evolution stopped by an observable delivers fewer samples than were
* planned, and the dataset has to end where the sampling did rather than keep
* rows that were never written.
*/
static void truncatesOnEarlyStop()
{
	const size_t dim = 16;
	const double maxT = 1.0, samplingStep = 0.1;
	const size_t stopAfter = 4;

	std::vector<std::complex<double>> initial(dim, 0.0);
	initial[0] = 1.0;

	std::vector<std::unique_ptr<krylovBasicObservable>> observables;
	observables.push_back(std::make_unique<stoppingObservable>(stopAfter));

	krylovTimeEvolver evolver(maxT, initial.data(), samplingStep, std::move(observables), hoppingChain(dim),
		1.0, 1e-8, 10, false, false);

	parameter_list parameters;
	hdf5ResultWriter writer(parameters, "testWavefunctionStreamingStop");
	const std::string file = writer.fileName();
	evolver.setSampleWriter(&writer);

	std::unique_ptr<krylovReturn> result(evolver.timeEvolve());
	writer.writeMetadata(*result);

	hsize_t numRows = 0, numColumns = 0;
	readStates(file, numRows, numColumns);

	check("an interrupted evolution reports the stop", result->statusCode == 3 || result->numSamples == stopAfter);
	check("the dataset holds only the samples that were taken", numRows == (hsize_t)stopAfter);
	check("the dataset agrees with the reported sample count", numRows == (hsize_t)result->numSamples);

	std::remove(file.c_str());
}

/**
* Without a writer nothing is written, and in particular the time evolution does
* not start collecting states on its own.
*/
static void writesNothingWithoutAWriter()
{
	const size_t dim = 8;
	std::vector<std::complex<double>> initial(dim, 0.0);
	initial[0] = 1.0;

	std::vector<std::unique_ptr<krylovBasicObservable>> observables;
	observables.push_back(std::make_unique<krylovSpMatrixObservable>("hopping", hoppingChain(dim)));

	krylovTimeEvolver evolver(0.5, initial.data(), 0.1, std::move(observables), hoppingChain(dim),
		1.0, 1e-8, 8, false, false);

	std::unique_ptr<krylovReturn> result(evolver.timeEvolve());

	parameter_list parameters;
	saveResult(*result, parameters, "testWavefunctionStreamingNone");

	const std::string file = "testWavefunctionStreamingNone.h5";
	H5::H5File handle(file, H5F_ACC_RDONLY);

	//The observable is there, which distinguishes a file without a wavefunction
	//from a file that was never written
	check("the observables are written as before", handle.nameExists("hopping"));
	check("no wavefunction dataset is created unless one was asked for", !handle.nameExists("wavefunction"));

	handle.close();
	std::remove(file.c_str());
}

/**
* A model parameter whose name is one the writer uses for its own bookkeeping
* has to be refused with an explanation rather than with a bare HDF5 failure.
*/
static void rejectsReservedParameterNames()
{
	parameter_list parameters;
	parameters.push_back(paraPush("dim", false, 1));

	try
	{
		hdf5ResultWriter writer(parameters, "testWavefunctionStreamingReserved");
		check("a parameter colliding with a metadata attribute is refused", false);
	}
	catch (const TE::krylovIOError&)
	{
		check("a parameter colliding with a metadata attribute is refused", true);
	}

	std::remove("testWavefunctionStreamingReserved.h5");
}

#endif

int main()
{
#ifdef USE_HDF
	H5::Exception::dontPrint();
	streamsEverySample();
	truncatesOnEarlyStop();
	writesNothingWithoutAWriter();
	rejectsReservedParameterNames();
#else
	std::cout << "built without HDF5, nothing to check" << std::endl;
#endif

	if (failures != 0)
	{
		std::cerr << failures << " wavefunction streaming check(s) failed" << std::endl;
		return 1;
	}
	std::cout << "All wavefunction streaming checks passed." << std::endl;
	return 0;
}
