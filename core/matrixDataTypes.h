#pragma once

#include <complex>
#include <string>
#include <vector>
#include <algorithm>


#ifdef USE_HDF
    #include <H5Cpp.h>
    using namespace H5;
#endif

#ifdef USE_MKL
#define MKL_Complex16 std::complex<double>
#define MKL_INT size_t
#include <mkl.h>
#include <mkl_spblas.h>

#elif defined USE_OPENBLAS
#include <cblas.h>
//#include <lapacke.h>
#endif


    //Define namespace for matrices and vector classes
    namespace TE {

        //Matrix class
        class matrix
        {
        public:
            std::complex<double>* values;
            size_t n, m;
            size_t numValues;

            matrix(size_t nn, size_t mm);
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

            vector(unsigned int n){
                length = n;
                values = new std::complex<double>[length];
            }

            ~vector(){
                delete[] values;
            }
            //implicit conversion operator to pointer
            operator std::complex<double>*() const { return NULL; }
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

            double norm1();
            double normInf();

            smatrix();
            smatrix(std::complex<double>* val, size_t* col, size_t* row, size_t nbV, unsigned int nn, unsigned int mm);
            ~smatrix();

            //inline int spMV(std::complex<double> alpha, vector* in, vector* out); //future implementation for GPU
            int spMV(std::complex<double> alpha, std::complex<double>* in, std::complex<double>* out);

            int initialize();



            static constexpr std::complex<double> one = std::complex<double>(1.0, 0.0);
            static constexpr std::complex<double> zero = std::complex<double>(0.0, 0.0);
            
#ifdef USE_MKL
            //variables for mkl-library
            sparse_matrix_t* MKLSparseMatrix;
            matrix_descr descriptor;
#endif


#ifdef USE_HDF
            int dumpHDF5(std::string fileName);
#endif


        };

    }