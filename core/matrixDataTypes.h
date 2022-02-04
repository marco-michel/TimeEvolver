#pragma once

#include <complex>
#include <string>
#include <vector>
#include <algorithm>

#define MKL_Complex16 std::complex<double>
#define MKL_INT size_t

#include <mkl.h>
#include <mkl_spblas.h>

#include <H5Cpp.h>

using namespace H5;


//Matrix class
class matrix
{
public:
	std::complex<double>* values;
	size_t n, m;
	size_t numValues;

    matrix(size_t nn, size_t mm)
	{
        if(nn*mm >0)
        {
            n = nn; m = mm;
            values = new std::complex<double>[n*m];
            numValues = n*m;
        }
	}

	~matrix()
	{
        if(n*m >0)
        {
            delete[] values;
        }
	}
    
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


//Sparse matrix class
class smatrix
{
public:
	std::complex<double>* values;
	size_t* columns;
	size_t* rowIndex;
	size_t numValues;
	int n, m;
	bool sym, hermitian;
    bool upperTri;

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

		for (unsigned int i = 0; i != nbV; i++)
		{
			values[i] = val[i];
			columns[i] = col[i];
			rowIndex[i] = row[i];
		}
	}

    
    int dumpHDF5(std::string fileName)
    {
        
        double* realPart = new double[numValues];
        double* imagPart = new double[numValues];
        
        for (unsigned int i = 0; i != numValues; i++)
        {
            realPart[i] = values[i].real();
            imagPart[i] = values[i].imag();
        }
        
        std::string fileNameH5 = fileName;
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
        
        IntType datatypeInt( PredType::NATIVE_HSIZE);
        datatypeInt.setOrder(H5T_ORDER_LE );
        DataSet dataset3 = fileHh.createDataSet("columIndex", datatypeInt,
                                               dataspace );
        dataset3.write(columns, PredType::NATIVE_HSIZE );
        DataSet dataset4 = fileHh.createDataSet("rowIndex", datatypeInt,
                                                dataspace );
        dataset4.write(rowIndex, PredType::NATIVE_HSIZE );
        
        delete[] realPart;
        delete[] imagPart;
        
        return 0;
    }

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
