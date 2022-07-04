#pragma once

#include <complex>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

#define MKL_Complex16 std::complex<double>
#define MKL_INT size_t

#ifdef USE_MKL
    #include <mkl.h>
    #include <mkl_spblas.h>
#endif
#ifdef USE_HDF
    #include <H5Cpp.h>
    using namespace H5;
#endif
#ifdef USE_CUDA
    #include <cuda_runtime_api.h>
    #include <cusparse.h>
#endif


#define CHECK_CUDA(func)                                                       \
{                                                                              \
    cudaError_t status = (func);                                               \
    if (status != cudaSuccess) {                                               \
        printf("CUDA API failed at line %d with error: %s (%d)\n",             \
               __LINE__, cudaGetErrorString(status), status);                  \
        return EXIT_FAILURE;                                                   \
    }                                                                          \
}

#define CHECK_CUSPARSE(func)                                                   \
{                                                                              \
    cusparseStatus_t status = (func);                                          \
    if (status != CUSPARSE_STATUS_SUCCESS) {                                   \
        printf("CUSPARSE API failed at line %d with error: %s (%d)\n",         \
               __LINE__, cusparseGetErrorString(status), status);              \
        return EXIT_FAILURE;                                                   \
    }                                                                          \
}


//Matrix class
class matrix
{
public:
	alignas(16) std::complex<double>* values;
	size_t n, m;
	size_t numValues;

    matrix(size_t nn, size_t mm)
	{
        if(nn*mm >0)
        {
            n = nn; m = mm;
            numValues = n * m;
#ifdef USE_CUDA
            cudaMallocHost((void**) &values, n * m * sizeof(std::complex<double>));
#else
            values = new std::complex<double>[n * m];
#endif
        }
	}

	~matrix()
	{
        if(n*m >0)
        {
#ifdef USE_CUDA
            cudaFreeHost(values);
#else
            delete[] values;
#endif
        }
	}
#ifdef USE_HDF 
    int dumpHDF5(std::string filename)
    {
        double* realPart = new double[numValues];
        double* imagPart = new double[numValues];
        
        for (unsigned int i = 0; i != numValues; i++)
        {
            realPart[i] = values[i].real();
            imagPart[i] = values[i].imag();
        }
        
        std::string fileNameH5 = filename;
        H5File fileHh(fileNameH5, H5F_ACC_TRUNC);
        int NX = numValues;
        const int RANK = 1;
        hsize_t dimsf[RANK];
        dimsf[0] = NX;
        DataSpace dataspace( RANK, dimsf );
        FloatType datatype( PredType::NATIVE_DOUBLE );
        datatype.setOrder( H5T_ORDER_LE );
        DataSet dataset1 = fileHh.createDataSet("valuesRealPart", datatype,
                                                dataspace );
        dataset1.write(realPart, PredType::NATIVE_DOUBLE );
        DataSet dataset2 = fileHh.createDataSet("valuesImagPart", datatype,
                                                dataspace );
        dataset2.write(imagPart, PredType::NATIVE_DOUBLE );
        
        delete[] realPart;
        delete[] imagPart;
        
        return 0;
    }
#endif    
};

//Vector class
class vector
{
public:
	std::complex<double>* values;
	size_t length;

	vector(int n)
	{
		length = n;
		values = new std::complex<double>[length];
	}

	~vector()
	{
		delete[] values;
	}
};

enum libraryType { MKL, CUDA };
//Sparse matrix class
class smatrix
{
public:
	alignas(16) std::complex<double>* values;
	size_t* columns;
	size_t* rowIndex;
	size_t numValues;
	int n, m;
	bool sym, hermitian;
    bool upperTri;

    static constexpr std::complex<double> one = std::complex<double>(1.0, 0.0);
    static constexpr std::complex<double> zero = std::complex<double>(0.0, 0.0);

    libraryType lType = MKL;

#ifdef USE_MKL
    sparse_matrix_t* A;
    matrix_descr descriptor;
#endif

#ifdef USE_CUDA

    alignas(8) long long *rowIndexCUDA, *colIndexCUDA;
    cusparseHandle_t     handle;
    cusparseSpMatDescr_t matA;
    long long *dA_rows, *dA_columns;
    alignas(16) std::complex<double> *dA_values;
    alignas(16) std::complex<double> *dX, *dY;
    cusparseDnVecDescr_t vecX, vecY;
    void* dBuffer = NULL;
    size_t bufferSize = 0;

    cudaError_t CUDAstatus;
    cusparseStatus_t CUDASpstatus;

#endif

	smatrix()
	{
		sym = hermitian = false;
		rowIndex = columns = nullptr;
		values = nullptr;
		m = n = 0;
		numValues = 0;
	}

    double norm1()
    {
        std::vector<double> colVal(n);
        std::vector<double>::iterator result;

        for(unsigned int i = 0; i != numValues; i++)
        {
            colVal[columns[i]] += std::abs(values[i]);
        }
        
        result = std::max_element(colVal.begin(),colVal.end());
        return *result;
    }

    void createLibraryType(libraryType a)
    {
        lType = a;

#ifdef USE_MKL
        if (a == MKL)
        {
            sparse_status_t mklStatus;
            matrix_descr type; type.type = SPARSE_MATRIX_TYPE_GENERAL;
            descriptor.type = SPARSE_MATRIX_TYPE_GENERAL;

            A = new sparse_matrix_t;
            if (numValues != 0)
                mklStatus = mkl_sparse_z_create_coo(A, SPARSE_INDEX_BASE_ZERO, m, n, numValues, rowIndex, columns, values);

            mklStatus = mkl_sparse_convert_csr(*A, SPARSE_OPERATION_NON_TRANSPOSE, A);
            mklStatus = mkl_sparse_order(*A);
        }
#endif
#ifdef USE_CUDA
        if (a == CUDA)
        {

            rowIndexCUDA = new long long[numValues];
            colIndexCUDA = new long long[numValues];

            for (size_t i = 0; i != numValues; i++)
            {
                rowIndexCUDA[i] = rowIndex[i];
                colIndexCUDA[i] = columns[i];
            }

            
            CUDAstatus = cudaMalloc((void**)&dA_rows, numValues * sizeof(long long));
            CUDAstatus = cudaMalloc((void**)&dA_columns, numValues * sizeof(long long));
            CUDAstatus = cudaMalloc((void**)&dA_values, numValues * sizeof(std::complex<double>));

            CUDAstatus = cudaMalloc((void**)&dX, n * sizeof(std::complex<double>));
            CUDAstatus = cudaMalloc((void**)&dY, n * sizeof(std::complex<double>));

            CUDAstatus = cudaMemcpy(dA_values, values, numValues * sizeof(std::complex<double>), cudaMemcpyHostToDevice);

            std::complex<double>* tests = new std::complex<double>[numValues];
            CUDAstatus = cudaMemcpy(tests, dA_values, numValues * sizeof(std::complex<double>), cudaMemcpyDeviceToHost);

            CUDAstatus = cudaMemcpy(dA_rows, rowIndex, numValues * sizeof(long long), cudaMemcpyHostToDevice);
            CUDAstatus = cudaMemcpy(dA_columns, columns, numValues * sizeof(long long), cudaMemcpyHostToDevice);

            CUDASpstatus = cusparseCreate(&handle);

            CUDASpstatus = cusparseCreateCoo(&matA, n, m, numValues, dA_rows, dA_columns, dA_values, CUSPARSE_INDEX_64I, CUSPARSE_INDEX_BASE_ZERO, CUDA_C_64F);
            CUDASpstatus = cusparseCreateDnVec(&vecX, n, dX, CUDA_C_64F);
            CUDASpstatus = cusparseCreateDnVec(&vecY, n, dY, CUDA_C_64F);
        }
#endif
    }

    inline int sparseSpMV(std::complex<double> alpha, std::complex<double>* X, std::complex<double>* Y)
    {
#ifdef USE_MKL
        if (lType == MKL)
        {
            std::complex<double>* result = 0;
            sparse_status_t mklStatus = mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, *A,
                descriptor, X, zero, Y);
        }

#endif

#ifdef USE_CUDA
        if (lType == CUDA)
        {

            CUDAstatus = cudaMemcpy(dX, X, n * sizeof(std::complex<double>), cudaMemcpyHostToDevice);


            alignas(16) std::complex<double> zeroCUDA;
            alignas(16) std::complex<double> oneCUDA; 
            oneCUDA.real(1);
            alignas(16) std::complex<double> alphaCUDA = alpha;


            CUDASpstatus = cusparseSpMV_bufferSize(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alphaCUDA, matA, vecX, &zeroCUDA, vecY, CUDA_C_64F, CUSPARSE_MV_ALG_DEFAULT, &bufferSize);
            CUDAstatus = cudaMalloc(&dBuffer, bufferSize);



            CUDASpstatus = cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alphaCUDA, matA, vecX, &zeroCUDA, vecY, CUDA_C_64F, CUSPARSE_MV_ALG_DEFAULT, dBuffer);

            CUDAstatus = cudaMemcpy(Y, dY, n * sizeof(std::complex<double>), cudaMemcpyDeviceToHost);



        }
#endif
        return 0;
    }

    double normInf()
    {
        std::vector<double> rowVal(m);
        std::vector<double>::iterator result;

        for(unsigned int i = 0; i != numValues; i++)
        {
            rowVal[rowIndex[i]] += std::abs(values[i]);
        }

        result = std::max_element(rowVal.begin(),rowVal.end());
        return *result;
    }


	smatrix(std::complex<double>* val, size_t* col, size_t* row, size_t nbV, int nn, int mm)
	{
		numValues = nbV; n = nn; m = mm;
		sym = hermitian = false;
		columns = new size_t[nbV];
		rowIndex = new size_t[nbV];
		values = new std::complex<double>[nbV];
        lType = MKL;

		for (unsigned int i = 0; i != nbV; i++)
		{
			values[i] = val[i];
			columns[i] = col[i];
			rowIndex[i] = row[i];
		}
	}

    smatrix(std::complex<double>* val, size_t* col, size_t* row, size_t nbV, int nn, int mm, libraryType type)
    {
        numValues = nbV; n = nn; m = mm;
        sym = hermitian = false;
        columns = new size_t[nbV];
        rowIndex = new size_t[nbV];
        values = new std::complex<double>[nbV];
        lType = type;

        for (unsigned int i = 0; i != nbV; i++)
        {
            values[i] = val[i];
            columns[i] = col[i];
            rowIndex[i] = row[i];
        }
    }

#ifdef USE_HDF
    int dumpHDF5(std::string fileName)
    {
        
        double* realPart = new double[numValues];
        double* imagPart = new double[numValues];
        
        for (unsigned int i = 0; i != numValues; i++)
        {
            realPart[i] = values[i].real();
            imagPart[i] = values[i].imag();
        }
        
        size_t mExport = (size_t)m;
        std::string fileNameH5 = fileName;
        H5File fileHh(fileNameH5, H5F_ACC_TRUNC);
        int NX = numValues;
        const int RANK = 1;
        hsize_t dimsf[RANK];
        hsize_t dimsatt[RANK];
        dimsf[0] = NX;
        dimsatt[0] = 1;
        DataSpace dataspace( RANK, dimsf );
        FloatType datatype( PredType::NATIVE_DOUBLE );

        DataSpace dataspace2(RANK, dimsatt);

        datatype.setOrder( H5T_ORDER_LE );
        DataSet dataset1 = fileHh.createDataSet("valuesRealPart", datatype,
                                               dataspace );
        dataset1.write(realPart, PredType::NATIVE_DOUBLE );
        DataSet dataset2 = fileHh.createDataSet("valuesImagPart", datatype,
                                               dataspace );
        dataset2.write(imagPart, PredType::NATIVE_DOUBLE );
        
        IntType datatypeInt( PredType::NATIVE_HSIZE);
        datatypeInt.setOrder(H5T_ORDER_LE );
        DataSet dataset3 = fileHh.createDataSet("columIndex", datatypeInt,
                                               dataspace );
        dataset3.write(columns, PredType::NATIVE_HSIZE );
        DataSet dataset4 = fileHh.createDataSet("rowIndex", datatypeInt,
                                                dataspace );
        dataset4.write(rowIndex, PredType::NATIVE_HSIZE );

        DataSet dataset5 = fileHh.createDataSet("dimension", datatypeInt, dataspace2);
        dataset5.write(&mExport, PredType::NATIVE_HSIZE);
        
        delete[] realPart;
        delete[] imagPart;
        
        return 0;
    }
#endif
	~smatrix()
	{
		if (n > 1 || m > 1)
		{
			delete[] rowIndex;
			delete[] values;
            delete[] columns;
			
		 
		}
	}

};
