#include <gtest/gtest.h>
#include <complex>
#include <cstddef>

#include "matrixDataTypes.h"  



TEST(smatrixEqual, trivialMatricsTests)
{
    const size_t problemSize = 20;

    std::complex<double> *vals = new std::complex<double>[problemSize];
    size_t *col = new size_t[problemSize];
    size_t *row = new size_t[problemSize];

    for (int i = 0; i != problemSize; i++)
    {
        col[i] = row[i] = i;
        vals[i] = std::complex<double> (i, 0.0);
    }


    TE::smatrix A(vals, col, row, problemSize, 10, 10);
    TE::smatrix A2(vals, col, row, problemSize, 10, 10);
    TE::smatrix B(vals, col, row, problemSize, 15, 10);

    EXPECT_FALSE(A.approxEqual(B));
    EXPECT_TRUE(A.approxEqual(A2));
    
    vals[0] = std::complex<double>(1.0,1.0);
    TE::smatrix C(vals, col, row, problemSize, 10, 10);
    EXPECT_FALSE(A.approxEqual(C));

    delete[] row;
    delete[] col;
    delete[] vals;
}

