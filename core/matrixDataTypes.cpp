#include "matrixDataTypes.h"

using namespace TE;

/*
* Constructor for empty matrix 
* @param nn column dimension
* @param mm row dimension
*/
matrix::matrix(size_t nn, size_t mm)
{
    n = nn; m = mm;
    numValues = n * m;
    if (nn * mm > 0)
        values = new std::complex<double>[n * m];
    else
        values = nullptr;
}

/*
* Constructor for matrix 
* @param nn column dimension
* @param mm row dimension
* @param vals values to initialize the matrix with. Note there are no checks for array size. nn*mm is assumed.
*/
TE::matrix::matrix(size_t nn, size_t mm, std::complex<double>* vals) : n(nn), m(mm)
{
    numValues = n * m;
    values = new std::complex<double>[numValues];
    cblas_zcopy(numValues, vals, 1, values, 1);
}

/*
* Deconstructor for matrix
*/
matrix:: ~matrix()
{
    if (n * m > 0)
    {
        delete[] values;
    }
}


#ifdef USE_HDF 
/*
* Save matrix to HDF5 file. 
* @param filename Filename 
* @return The return value has no meaning. It always return 0.
*/
int matrix::dumpHDF5(std::string filename)
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
    size_t NX = numValues;
    const int RANK = 1;
    hsize_t dimsf[RANK];
    dimsf[0] = NX;
    DataSpace dataspace(RANK, dimsf);
    FloatType datatype(PredType::NATIVE_DOUBLE);
    datatype.setOrder(H5T_ORDER_LE);
    DataSet dataset1 = fileHh.createDataSet("valuesRealPart", datatype,
        dataspace);
    dataset1.write(realPart, PredType::NATIVE_DOUBLE);
    DataSet dataset2 = fileHh.createDataSet("valuesImagPart", datatype,
        dataspace);
    dataset2.write(imagPart, PredType::NATIVE_DOUBLE);

    delete[] realPart;
    delete[] imagPart;

    return 0;
}
#endif    


/*
* Default constructor for sparse matrix. Initializes also library specific variables.
*/
smatrix::smatrix()
{
    sym = hermitian = upperTri = initialized = false;
    rowIndex = columns = nullptr;
    values = nullptr;
    m = n = 0;
    numValues = 0;
#ifdef USE_MKL
    //variables for mkl-library
    MKLSparseMatrix = new sparse_matrix_t;
    descriptor.type = SPARSE_MATRIX_TYPE_GENERAL;
    descriptor.diag = SPARSE_DIAG_NON_UNIT;
    initialize();
#endif
}

/**
* Compute 1-norm: maximum absolute column sum
* @return 1-norm of the matrix
*/
double smatrix::norm1()
{
    std::vector<double> colVal(n);
    std::vector<double>::iterator result;

    for (unsigned int i = 0; i != numValues; i++)
    {
        colVal[columns[i]] += std::abs(values[i]);
    }

    result = std::max_element(colVal.begin(), colVal.end());
    return *result;
}

/*
* Compute infinity-norm: maximum absolute row sum
* @return infinity-norm of the matrix
*/
double smatrix::normInf()
{
    std::vector<double> rowVal(m);
    std::vector<double>::iterator result;

    for (unsigned int i = 0; i != numValues; i++)
    {
        rowVal[rowIndex[i]] += std::abs(values[i]);
    }

    result = std::max_element(rowVal.begin(), rowVal.end());
    return *result;
}

/*
* Constructor for sparse matrix with initializing values
* @param val Values to initialize the (sparse) matrix with
* @param col Column indices for non-zero values 
* @param row Row indices for non-zero values
* @param nbV Number of non-zero values
* @param nn Column dimension
* @param mm Row dimension
*/
smatrix::smatrix(std::complex<double>* val, size_t* col, size_t* row, size_t nbV, unsigned int nn, unsigned int mm)
{
    if (nn == 0 || mm == 0) {
        std::cerr << "Empty matrices are not supported." << std::endl;
        exit(1);
    }
    numValues = nbV; n = nn; m = mm;
    sym = hermitian = upperTri = false;
    columns = new size_t[nbV];
    rowIndex = new size_t[nbV];
    values = new std::complex<double>[nbV];
    initialized = false;

    for (unsigned int i = 0; i != nbV; i++)
    {
        values[i] = val[i];
        columns[i] = col[i];
        rowIndex[i] = row[i];
    }
#ifdef USE_MKL
    //variables for mkl-library
    MKLSparseMatrix = new sparse_matrix_t;
    descriptor.type = SPARSE_MATRIX_TYPE_GENERAL;
    descriptor.diag = SPARSE_DIAG_NON_UNIT;
    initialize();
#endif
}

/*
* Copy constructor
* @param old_obj Reference object to construct a copy from
*/
smatrix::smatrix(const smatrix& old_obj) {
    numValues = old_obj.numValues; n = old_obj.n; m = old_obj.m;
    sym = old_obj.sym; hermitian = old_obj.hermitian; upperTri = old_obj.upperTri;
    columns = new size_t[numValues];
    rowIndex = new size_t[numValues];
    values = new std::complex<double>[numValues];
    initialized = false;

    for (unsigned int i = 0; i != numValues; i++) {
        values[i] = old_obj.values[i];
        columns[i] = old_obj.columns[i];
        rowIndex[i] = old_obj.rowIndex[i];
    }

#ifdef USE_MKL
    //variables for mkl-library
    MKLSparseMatrix = new sparse_matrix_t;
    descriptor = old_obj.descriptor;
    initialize();
#endif
}

/*
* Sparse Matrix dense vector multiplication: out = alpha * this * in. Uses MKL for optimized routines but also has an alternative branch in case there is no MKL installed. 
* @param alpha Scalar factor (usually set to one)
* @param in (Dense) vector multiplying (this) matrix
* @param out Result sparse matrix
* @return Status indicating success or failure of the operation
*/
int smatrix::spMV(std::complex<double> alpha, std::complex<double>* in, std::complex<double> *out) {

#if defined USE_MKL
    sparse_status_t mklStatus = mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, *MKLSparseMatrix,
        descriptor, in, zero, out);
    return (int) mklStatus;
#elif defined USE_ARMADILLO
    arma::cx_vec ArmadilloVector(in, m, false, true);
    auto res = alpha * (*ArmadilloSparseMatrix) * ArmadilloVector;
    for (unsigned int i = 0; i != m; i++)
        out[i] = res(i);
    return 0;

#else  //not recommended, very slow,  please install an optimized sparse BLAS library for reasonable performance
    for (size_t i = 0; i != n; i++)
        out[i] = 0;
    for (size_t i = 0; i < numValues; i++) {
        out[rowIndex[i]] += alpha * values[i] * in[columns[i]];
    }
    return 0;
#endif
}

/*
* Setting up library variables to speed up spMV. Or do nothing if none is installed.
*/
int smatrix::initialize() {

    if (initialized == true)
        return -1;

#ifdef USE_MKL
    if (numValues == 0)
        return 1;

    sparse_status_t mklStatus;
    matrix_descr type; type.type = SPARSE_MATRIX_TYPE_GENERAL; type.diag = SPARSE_DIAG_NON_UNIT; type.mode = SPARSE_FILL_MODE_FULL;

    descriptor.type = SPARSE_MATRIX_TYPE_GENERAL;
    descriptor.diag = SPARSE_DIAG_NON_UNIT;

    MKLSparseMatrix = new sparse_matrix_t;

    mklStatus = mkl_sparse_z_create_coo(MKLSparseMatrix, SPARSE_INDEX_BASE_ZERO, m, n, numValues, rowIndex, columns, values);

    if (mklStatus != SPARSE_STATUS_SUCCESS) {
        std::cerr << "Problem with MKL sparse matrix creation" << std::endl;
        return -1;
    }
    
    mkl_sparse_convert_csr(*MKLSparseMatrix, SPARSE_OPERATION_NON_TRANSPOSE, MKLSparseMatrix);
    mkl_sparse_order(*MKLSparseMatrix);
    mkl_sparse_set_mv_hint(*MKLSparseMatrix, SPARSE_OPERATION_NON_TRANSPOSE, type, 20000);
    mkl_sparse_set_memory_hint(*MKLSparseMatrix, SPARSE_MEMORY_AGGRESSIVE);
    mkl_sparse_optimize(*MKLSparseMatrix);
    
 #endif
    
#ifdef USE_ARMADILLO
    ArmadillorowIndex = arma::umat((unsigned long long *) rowIndex, 1, numValues, false, true);
    ArmadillocolIndex = arma::umat((unsigned long long *) columns, 1, numValues, false, true);
    ArmadillovalueVector = arma::cx_vec(values, numValues, false, true);
    ArmadilloindexMatrix = arma::join_cols(ArmadillorowIndex, ArmadillocolIndex);

    ArmadilloSparseMatrix = new arma::sp_cx_mat(ArmadilloindexMatrix, ArmadillovalueVector, m, n, true, true);
#endif

    initialized = true;

    return 0;
}


#ifdef USE_HDF
/*
* Save matrix to HDF5 file.
* @param filename Filename
* @return The return value has no meaning. It always return 0.
*/
int smatrix::dumpHDF5(std::string fileName)
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
    size_t NX =  numValues;
    const int RANK = 1;
    hsize_t dimsf[RANK] = {};
    hsize_t dimsatt[RANK] = {};
    dimsf[0] = NX;
    dimsatt[0] = 1;
    DataSpace dataspace(RANK, dimsf);
    FloatType datatype(PredType::NATIVE_DOUBLE);

    DataSpace dataspace2(RANK, dimsatt);

    datatype.setOrder(H5T_ORDER_LE);
    DataSet dataset1 = fileHh.createDataSet("valuesRealPart", datatype,
        dataspace);
    dataset1.write(realPart, PredType::NATIVE_DOUBLE);
    DataSet dataset2 = fileHh.createDataSet("valuesImagPart", datatype,
        dataspace);
    dataset2.write(imagPart, PredType::NATIVE_DOUBLE);

    IntType datatypeInt(PredType::NATIVE_HSIZE);
    datatypeInt.setOrder(H5T_ORDER_LE);
    DataSet dataset3 = fileHh.createDataSet("columIndex", datatypeInt,
        dataspace);
    dataset3.write(columns, PredType::NATIVE_HSIZE);
    DataSet dataset4 = fileHh.createDataSet("rowIndex", datatypeInt,
        dataspace);
    dataset4.write(rowIndex, PredType::NATIVE_HSIZE);

    DataSet dataset5 = fileHh.createDataSet("dimension", datatypeInt, dataspace2);
    dataset5.write(&mExport, PredType::NATIVE_HSIZE);

    delete[] realPart;
    delete[] imagPart;

    return 0;
}
#endif

/*
* Default deconstructor
*/
smatrix::~smatrix()
{
    if (n > 1 || m > 1)
    {
        delete[] rowIndex;
        delete[] values;
        delete[] columns;
    }

#ifdef USE_MKL
    if (MKLSparseMatrix != nullptr) {
        mkl_sparse_destroy(*MKLSparseMatrix);
        delete MKLSparseMatrix;
    }
#endif
}


