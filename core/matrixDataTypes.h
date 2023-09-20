#pragma once

#include <complex>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

#include "mathHeader.h"


#ifdef USE_HDF
    #include <H5Cpp.h>
    using namespace H5;
#endif






    //Define namespace for matrices and vector classes
    namespace TE {

        //Matrix class
        class matrix
        {
        public:
            size_t n, m;
            std::complex<double>* values;
            size_t numValues;
            matrix(size_t nn, size_t mm);
            matrix(size_t nn, size_t mm, std::complex<double>* vals);
            ~matrix();
#ifdef USE_HDF 
            int dumpHDF5(std::string filename);
#endif    
        };


        //Vector class
        class vector
        {
        public:
            std::complex<double>* values;
            size_t length;

            vector(unsigned int n) {
                if (n == 0) {
                    std::cerr << "Empty vectors are not supported." << std::endl;
                    exit(1);
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


        //Sparse matrix class
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

            int spMV(std::complex<double> alpha, std::complex<double>* in, std::complex<double>* out);
            int initialize();

            static constexpr std::complex<double> one = std::complex<double>(1.0, 0.0);
            static constexpr std::complex<double> zero = std::complex<double>(0.0, 0.0);

#ifdef USE_MKL
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
    
#ifdef USE_HDF
            int dumpHDF5(std::string fileName);
#endif
        };


    }
        
    //Wrapper for element-wise vector operations: exp
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

    //Wrapper for element-wise vector operations: multiplication
    inline void mulV(size_t len, std::complex<double>* a, std::complex<double>* b, std::complex<double>* y) {
#ifdef USE_MKL
        vzMul(len, a, b, y);
#else
        for (size_t i = 0; i < len; i++) {
            y[i] = b[i] * a[i];
        }
#endif
    }

    //Wrapper for LAPACK zhseqr

    inline int TE_zhseqr(size_t  m,
		std::complex<double> *  	h,
		std::complex<double> *  	w,
		std::complex<double> *  	z
	) 	
    {
#ifdef ON_APPLE
    long info;
    long onez = 1;
    long workspaceSize = 11*m;
    const char job = 'S';
    const char COMPZ = 'I';
    long mReplace = (long) m;
    std::complex<double>* workspace = new std::complex<double>[workspaceSize];
    zhseqr_(&job, &COMPZ, &mReplace, &onez, &mReplace, h, &mReplace, w, z, &mReplace, workspace,  &workspaceSize, &info);
    std::cerr<< "apple lapack: " << std::endl;
    std::cerr << info << std::endl;
    delete[] workspace;
    return (int) info;
#else
    return LAPACKE_zhseqr(LAPACK_COL_MAJOR, 'S', 'I', m, (size_t) 1, m,
				h, m, w, z, m);
#endif
    }