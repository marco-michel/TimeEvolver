#pragma once

#include <complex>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

#ifdef USE_HDF
#include <hdf5.h>
#endif

#include "mathHeader.h"
#include "krylovExceptions.h"


    //Define namespace for matrices and vector classes
namespace TE {


    /**
    * Dense matrix class providing function to store the matrix as well as HDF5 file output on supported platforms
    */
    class matrix
    {
    public:
        size_t n, m;
        std::complex<double>* values;
        size_t numValues;
        matrix(size_t nn, size_t mm);
        matrix(size_t nn, size_t mm, std::complex<double>* vals);
        ~matrix();
    };


    /**
    * Class to store a dense vector
    */
    class vector
    {
    public:
        std::complex<double>* values;
        size_t length;

        vector(unsigned int n) {
            if (n == 0) {
                throw krylovInvalidArgument("Empty vectors are not supported.");
            }
            length = n;
            values = new std::complex<double>[length];
        }

        ~vector() {
            delete[] values;
        }
        //implicit conversion operator to pointer
        operator std::complex<double>* () const { return NULL; }
    };


    /**
    * Internal class for sparse matrix representation. Provides a wrapper for sparse matrix - dense vector multiplications, which can take advantage of optimized libraries.
    */
    class smatrix
    {
    public:
        std::complex<double>* values;
        size_t* columns;
        size_t* rowIndex;
        size_t numValues;
        size_t n, m;
        bool sym, hermitian;
        bool upperTri;
        bool initialized;

        double norm1();
        double normInf();

        smatrix();
        smatrix(std::complex<double>* val, size_t* col, size_t* row, size_t nbV, unsigned int nn, unsigned int mm);
        smatrix(const smatrix& old_obj);
        ~smatrix();

        int spMV(std::complex<double> alpha, std::complex<double>* in, std::complex<double>* out) const;
        int initialize();

        //methods to store/load matrix to/from file
#ifdef USE_HDF
        void saveHDF5(const std::string& filename) const;
        void loadHDF5(const std::string& filename);
#endif

        static constexpr std::complex<double> one = std::complex<double>(1.0, 0.0);
        static constexpr std::complex<double> zero = std::complex<double>(0.0, 0.0);

#ifdef USE_MKL
    private:
        //variables for mkl-library
        sparse_matrix_t* MKLSparseMatrix;
        matrix_descr descriptor;
#endif

#ifdef USE_ARMADILLO
        //variables for armadillo-library
        arma::sp_cx_mat* ArmadilloSparseMatrix;
        arma::umat ArmadillorowIndex;
        arma::umat ArmadillocolIndex;
        arma::umat ArmadilloindexMatrix;
        arma::cx_vec ArmadillovalueVector;
#endif
    };



    /**
    * Wrapper for element-wise vector operations: exp
    */
    inline void expV(size_t len, std::complex<double>* x, std::complex<double>* y)
    {
#ifdef USE_MKL
        vzExp(len, x, y);
#else
        for (size_t i = 0; i < len; i++) {
            y[i] = std::exp(x[i]);
        }
#endif
    }

    /**
    * Wrapper for element-wise vector operations: multiplication
    */
    inline void mulV(size_t len, std::complex<double>* a, std::complex<double>* b, std::complex<double>* y) {
#ifdef USE_MKL
        vzMul(len, a, b, y);
#else
        for (size_t i = 0; i < len; i++) {
            y[i] = b[i] * a[i];
        }
#endif
    }

    /**
    * Wrapper for LAPACK zhseqr
    */
    inline size_t TE_zhseqr(size_t  m,
        std::complex<double>* h,
        std::complex<double>* w,
        std::complex<double>* z
    )
    {
#ifdef ON_APPLE
        long info;
        long onez = 1;
        long workspaceSize = 11 * m;
        const char job = 'S';
        const char COMPZ = 'I';
        long mReplace = (long)m;
        std::complex<double>* workspace = new std::complex<double>[workspaceSize];
        zhseqr_(&job, &COMPZ, &mReplace, &onez, &mReplace, h, &mReplace, w, z, &mReplace, workspace, &workspaceSize, &info);
        delete[] workspace;
        return (size_t)info;
#else
        return LAPACKE_zhseqr(LAPACK_COL_MAJOR, 'S', 'I', m, (size_t)1, m,
            h, m, w, z, m);
#endif
    }


#ifdef USE_HDF

    /**
    * Class to store complex numbers in legacy HDF5 1.14 format
    */
    struct hdf5_complex_t {
        double real;
        double imag;
    };


    // helper to create the compound datatype for complex numbers
    inline hid_t createComplexType()
    {
        hid_t complexType = H5Tcreate(H5T_COMPOUND, sizeof(hdf5_complex_t));
        H5Tinsert(complexType, "r", HOFFSET(hdf5_complex_t, real), H5T_NATIVE_DOUBLE);
        H5Tinsert(complexType, "i", HOFFSET(hdf5_complex_t, imag), H5T_NATIVE_DOUBLE);
        return complexType;
    }

    // small helper for scalar attributes
    template<typename T>
    inline void writeScalarAttribute(hid_t obj, const char* name, hid_t h5Type, const T& value)
    {
        hid_t space = H5Screate(H5S_SCALAR);
        hid_t attr = H5Acreate2(obj, name, h5Type, space, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, h5Type, &value);
        H5Aclose(attr);
        H5Sclose(space);
    }

    template<typename T>
    void readScalarAttribute(hid_t obj, const char* name, hid_t h5Type, T& out)
    {
        hid_t attr = H5Aopen(obj, name, H5P_DEFAULT);
        H5Aread(attr, h5Type, &out);
        H5Aclose(attr);
    }

#endif

}
