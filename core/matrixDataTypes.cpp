#include "matrixDataTypes.h"

using namespace TE;

/**
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

/**
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

/**
* Deconstructor for matrix
*/
matrix:: ~matrix()
{
    if (n * m > 0)
    {
        delete[] values;
    }
}

/**
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
    //variables for mkl-library. The handle is allocated by initialize(), which
    //is the only place that can fill it.
    MKLSparseMatrix = nullptr;
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

/**
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

/**
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
        throw krylovInvalidArgument("Empty matrices are not supported.");
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
    //variables for mkl-library. The handle is allocated by initialize(), which
    //is the only place that can fill it.
    MKLSparseMatrix = nullptr;
    descriptor.type = SPARSE_MATRIX_TYPE_GENERAL;
    descriptor.diag = SPARSE_DIAG_NON_UNIT;
    initialize();
#endif
}

/**
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
    //see the remark on ownership of MKLSparseMatrix in the other constructors
    MKLSparseMatrix = nullptr;
    descriptor = old_obj.descriptor;
    initialize();
#endif
}

/**
* Sparse Matrix dense vector multiplication: out = alpha * this * in. Uses MKL for optimized routines but also has an alternative branch in case there is no MKL installed. 
* @param alpha Scalar factor (usually set to one)
* @param in (Dense) vector multiplying (this) matrix
* @param out Result sparse matrix
* @return Status indicating success or failure of the operation
*/
int smatrix::spMV(std::complex<double> alpha, std::complex<double>* in, std::complex<double> *out) const {

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

/**
* Setting up library variables to speed up spMV. Or do nothing if none is installed.
* @return Status if the operation was successful (0) or not (error code != 0). 
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


/**
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



#ifdef USE_HDF

/**
* Save sparse matrix to a HDF5 file. The matrix is stored in COO form: one value,
* one column index and one row index per non zero entry. Everything describing
* the matrix itself is attached to the file as scalar attributes.
* @param filename Name of the output file
*/
void smatrix::saveHDF5(const std::string& filename) const
{
TE_HDF5_TRY
    H5::H5File file(filename, H5F_ACC_TRUNC);

    hdf5::writeVersion(file);
    hdf5::writeAttribute(file, "n", this->n);
    hdf5::writeAttribute(file, "m", this->m);
    hdf5::writeAttribute(file, "numValues", this->numValues);
    hdf5::writeAttribute(file, "sym", this->sym);
    hdf5::writeAttribute(file, "hermitian", this->hermitian);
    hdf5::writeAttribute(file, "upperTri", this->upperTri);

    //Stored as a compound type rather than as separate real and imaginary
    //datasets, so that a reader sees a single complex array.
    std::vector<hdf5::complexType> valueBuffer(this->numValues);
    for (size_t i = 0; i != this->numValues; i++)
    {
        valueBuffer[i].r = this->values[i].real();
        valueBuffer[i].i = this->values[i].imag();
    }

    std::vector<hsize_t> columnBuffer(this->numValues);
    std::vector<hsize_t> rowBuffer(this->numValues);
    for (size_t i = 0; i != this->numValues; i++)
    {
        columnBuffer[i] = static_cast<hsize_t>(this->columns[i]);
        rowBuffer[i] = static_cast<hsize_t>(this->rowIndex[i]);
    }

    H5::CompType complexType = hdf5::complexDataType();
    hdf5::writeDataset(file, "values", complexType, complexType,
        valueBuffer.data(), this->numValues, sizeof(hdf5::complexType));
    hdf5::writeDataset(file, "columns", H5::PredType::STD_U64LE, H5::PredType::NATIVE_HSIZE,
        columnBuffer.data(), this->numValues, sizeof(hsize_t));
    hdf5::writeDataset(file, "rowIndex", H5::PredType::STD_U64LE, H5::PredType::NATIVE_HSIZE,
        rowBuffer.data(), this->numValues, sizeof(hsize_t));
TE_HDF5_CATCH("Could not write sparse matrix to " + filename)
}

/**
* Load a sparse matrix from a HDF5 file written by saveHDF5. Any data held by
* this matrix is replaced.
* @param filename Name of the input file
*/
void smatrix::loadHDF5(const std::string& filename)
{
TE_HDF5_TRY
    H5::H5File file(filename, H5F_ACC_RDONLY);

    size_t newN = hdf5::readSizeAttribute(file, "n");
    size_t newM = hdf5::readSizeAttribute(file, "m");
    size_t newNumValues = hdf5::readSizeAttribute(file, "numValues");

    if (newN == 0 || newM == 0)
        throw krylovInvalidArgument("Empty matrices are not supported.");

    std::vector<hdf5::complexType> valueBuffer(newNumValues);
    std::vector<hsize_t> columnBuffer(newNumValues);
    std::vector<hsize_t> rowBuffer(newNumValues);

    if (newNumValues != 0)
    {
        H5::CompType complexType = hdf5::complexDataType();
        file.openDataSet("values").read(valueBuffer.data(), complexType);
        file.openDataSet("columns").read(columnBuffer.data(), H5::PredType::NATIVE_HSIZE);
        file.openDataSet("rowIndex").read(rowBuffer.data(), H5::PredType::NATIVE_HSIZE);
    }

    //Only replace the current contents once everything has been read, so that a
    //failed load leaves the matrix as it was.
    delete[] values;
    delete[] columns;
    delete[] rowIndex;

    numValues = newNumValues;
    n = newN;
    m = newM;
    sym = hdf5::readBoolAttribute(file, "sym");
    hermitian = hdf5::readBoolAttribute(file, "hermitian");
    upperTri = hdf5::readBoolAttribute(file, "upperTri");

    values = new std::complex<double>[numValues];
    columns = new size_t[numValues];
    rowIndex = new size_t[numValues];

    for (size_t i = 0; i != numValues; i++)
    {
        values[i] = std::complex<double>(valueBuffer[i].r, valueBuffer[i].i);
        columns[i] = static_cast<size_t>(columnBuffer[i]);
        rowIndex[i] = static_cast<size_t>(rowBuffer[i]);
    }

    initialized = false;
    initialize();
TE_HDF5_CATCH("Could not read sparse matrix from " + filename)
}

#endif
