/**
* Checks that storing only the upper triangle of the Hermitian Hamiltonian
* describes the same physics as storing all of it: same norm, same matrix
* vector product, and the same time evolution to the last digit that the
* tolerance can promise.
*/

#include <cmath>
#include <complex>
#include <iostream>
#include <memory>
#include <vector>

#include "Basis.h"
#include "exampleHamiltonian.h"
#include "hamiltonian.h"
#include "krylovExceptions.h"
#include "krylovObservables.h"
#include "krylovTimeEvolver.h"
#include "matrixDataTypes.h"

using namespace TE;

static int failures = 0;

static void check(const std::string& label, bool condition, const std::string& detail = "")
{
	if (condition)
		std::cout << "ok: " << label << (detail.empty() ? "" : " (" + detail + ")") << std::endl;
	else
	{
		std::cerr << "FAIL: " << label << (detail.empty() ? "" : " (" + detail + ")") << std::endl;
		failures++;
	}
}

/**
* Builds the black hole example small enough to run quickly but large enough to
* have a genuinely non trivial off-diagonal structure.
*/
static std::unique_ptr<smatrix> buildHamiltonian(bool hermitianStorage, tensorBasis& space)
{
	exampleHamiltonian ham(8, 4, 12.0, 4, 4, 1, 0.0003, 1.0, 1.0, 0.065, 0.065);
	ham.createSimplifiedHamiltonian();
	return ham.createHamiltonMatrix(&space, hermitianStorage);
}

int main()
{
	tensorBasis space(8, 2, 4, 8, 1);

	std::unique_ptr<smatrix> full = buildHamiltonian(false, space);
	std::unique_ptr<smatrix> half = buildHamiltonian(true, space);

	check("triangular storage sets the flags", half->hermitian && half->upperTri);
	check("full storage leaves the flags alone", !full->hermitian && !full->upperTri);
	check("the triangle is roughly half of the matrix",
		half->numValues < full->numValues && 2 * half->numValues > full->numValues,
		std::to_string(half->numValues) + " of " + std::to_string(full->numValues) + " entries");

	//The norms drive the step size and the round-off estimate, so they have to
	//come out the same even though half of the entries are not there
	check("norm1 agrees", std::abs(full->norm1() - half->norm1()) < 1e-12,
		std::to_string(full->norm1()) + " vs " + std::to_string(half->norm1()));
	check("normInf agrees", std::abs(full->normInf() - half->normInf()) < 1e-12,
		std::to_string(full->normInf()) + " vs " + std::to_string(half->normInf()));

	//A matrix vector product with an arbitrary vector
	const size_t dim = full->m;
	std::vector<std::complex<double>> in(dim), outFull(dim), outHalf(dim);
	for (size_t i = 0; i != dim; i++)
		in[i] = std::complex<double>(std::sin(0.7 * i), std::cos(0.3 * i));

	full->initialize();
	half->initialize();
	full->spMV(std::complex<double>(1.0, 0.0), in.data(), outFull.data());
	half->spMV(std::complex<double>(1.0, 0.0), in.data(), outHalf.data());

	double maxProductDifference = 0;
	for (size_t i = 0; i != dim; i++)
		maxProductDifference = std::max(maxProductDifference, std::abs(outFull[i] - outHalf[i]));
	check("the matrix vector product agrees", maxProductDifference < 1e-12,
		"max difference " + std::to_string(maxProductDifference));

	//And finally the whole time evolution
	auto evolve = [&](bool hermitianStorage) {
		tensorBasis localSpace(8, 2, 4, 8, 1);
		std::unique_ptr<smatrix> ham = buildHamiltonian(hermitianStorage, localSpace);
		exampleHamiltonian description(8, 4, 12.0, 4, 4, 1, 0.0003, 1.0, 1.0, 0.065, 0.065);
		description.createSimplifiedHamiltonian();
		basisVector init = description.createInitState();

		std::vector<std::complex<double>> state(localSpace.numberElements, std::complex<double>(0.0, 0.0));
		state[localSpace.hashTable.find(init)->second] = std::complex<double>(1.0, 0.0);

		krylovTimeEvolver evolver(10.0, state.data(), 1.0, {}, std::move(ham),
			1.0, 1e-8, 40, false, false);
		evolver.changeLogLevel(krylovLogger::FATAL);
		return std::unique_ptr<krylovReturn>(evolver.timeEvolve());
	};

	std::unique_ptr<krylovReturn> resultFull = evolve(false);
	std::unique_ptr<krylovReturn> resultHalf = evolve(true);

	double maxStateDifference = 0;
	for (size_t i = 0; i != resultFull->dim; i++)
		maxStateDifference = std::max(maxStateDifference,
			std::abs(resultFull->evolvedState[i] - resultHalf->evolvedState[i]));

	check("the evolved state agrees", maxStateDifference < 1e-10,
		"max difference " + std::to_string(maxStateDifference));
	check("the same number of steps was taken", resultFull->n_steps == resultHalf->n_steps,
		std::to_string(resultFull->n_steps) + " vs " + std::to_string(resultHalf->n_steps));
	check("the error bound agrees", std::abs(resultFull->err - resultHalf->err) < 1e-14);

	if (failures == 0)
		std::cout << "All Hermitian storage checks passed." << std::endl;
	return failures == 0 ? 0 : 1;
}
