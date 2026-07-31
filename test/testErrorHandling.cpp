/**
* Checks that the library reports failures by throwing rather than by
* terminating the process. Every case below used to call exit(1), which made
* these errors unrecoverable for any embedding application.
*/

#include <complex>
#include <iostream>
#include <memory>
#include <vector>

#include "krylovTimeEvolver.h"
#include "krylovObservables.h"
#include "krylovExceptions.h"
#include "matrixDataTypes.h"

static int failures = 0;

template <typename Fn>
static void expectThrow(const std::string& label, Fn fn)
{
	try
	{
		fn();
		std::cerr << "FAIL: " << label << " did not throw" << std::endl;
		failures++;
	}
	catch (const TE::krylovError& e)
	{
		std::cout << "ok: " << label << " -> " << e.what() << std::endl;
	}
	catch (const std::exception& e)
	{
		std::cerr << "FAIL: " << label << " threw an unexpected type: " << e.what() << std::endl;
		failures++;
	}
}

/**
* A hopping chain. It does not produce a lucky breakdown, so the error bound is
* actually exercised.
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

int main()
{
	expectThrow("empty matrix", [] {
		std::complex<double> value = 1.0; size_t column = 0, row = 0;
		TE::smatrix empty(&value, &column, &row, 1, 0, 0);
	});

	expectThrow("empty vector", [] { TE::vector empty(0); });

	expectThrow("unnormalized initial state", [] {
		const size_t dim = 8;
		std::vector<std::complex<double>> state(dim, 0.0);
		state[0] = 5.0;
		std::vector<std::unique_ptr<krylovBasicObservable>> observables;
		krylovTimeEvolver evolver(1.0, state.data(), 0.5, std::move(observables), hoppingChain(dim));
	});

	//A tolerance this small can never be met, so the step size reduction runs out
	//of attempts. This is a legitimate numerical outcome, not a misuse of the API.
	expectThrow("unreachable tolerance", [] {
		const size_t dim = 50;
		std::vector<std::complex<double>> state(dim, 0.0);
		state[0] = 1.0;
		std::vector<std::unique_ptr<krylovBasicObservable>> observables;
		krylovTimeEvolver evolver(5.0, state.data(), 1.0, std::move(observables), hoppingChain(dim),
			1.0, 1e-30, 40, false, false);
		delete evolver.timeEvolve();
	});

	//Reaching this line at all is the main assertion of this test: none of the
	//failures above terminated the process.
	if (failures != 0)
	{
		std::cerr << failures << " error handling check(s) failed" << std::endl;
		return 1;
	}

	std::cout << "All error handling checks passed." << std::endl;
	return 0;
}
