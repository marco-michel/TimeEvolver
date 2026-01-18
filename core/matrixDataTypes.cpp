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
* Compares two matrices modulo floating point noise. Flags are not checked. Differently sorted indices are not accounted for so far.
* @return true if matrices are approximitaley (up to fp noise) equal.
*/
bool smatrix::approxEqual(const smatrix& other, double absTol)
{
    //first check basic info
    if (other.m != this->m || other.n != this->n || other.numValues != this->numValues)
        return false;
    
    for(size_t i = 0; i != this->numValues; i++)
    {
        if(other.columns[i] != this->columns[i] || other.rowIndex[i] != this->rowIndex[i] || std::abs(other.values[i] - this->values[i]) > absTol)
            return false;        
    }
    //all checks passed therefore matrices are equal
    return true;
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
    //variables for mkl-library
    MKLSparseMatrix = new sparse_matrix_t;
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

#ifdef USE_HDF

/**
* Save sparse matrix to a HDF5 file
* @param M sparse matrix to be stored on file
* @param filename of the output file
*/
void smatrix::saveHDF5(const std::string& filename) const
{
    // --- convert complex values to hdf5_complex_t buffer ---
    std::vector<hdf5_complex_t> valBuf(this->numValues);
    for (size_t i = 0; i < this->numValues; ++i) {
        valBuf[i].real = this->values[i].real();
        valBuf[i].imag = this->values[i].imag();
    }

    // --- convert indices to hsize_t ---
    std::vector<hsize_t> colBuf(this->numValues);
    for (size_t i = 0; i < this->numValues; ++i) {
        colBuf[i] = static_cast<hsize_t>(this->columns[i]);
    }

    std::vector<hsize_t> rowBuf(this->numValues);
    for (size_t i = 0; i < this->numValues; ++i) {
        rowBuf[i] = static_cast<hsize_t>(this->rowIndex[i]);
    }

    // --- create file (overwrite if existing) ---
    hid_t file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC,
        H5P_DEFAULT, H5P_DEFAULT);

    // --- write metadata as attributes on the file ---
    hsize_t n_attr = static_cast<hsize_t>(this->n);
    hsize_t m_attr = static_cast<hsize_t>(this->m);
    hsize_t nnz_attr = static_cast<hsize_t>(this->numValues);
    int     sym_attr = this->sym ? 1 : 0;
    int     herm_attr = this->hermitian ? 1 : 0;
    int     upper_attr = this->upperTri ? 1 : 0;

    writeScalarAttribute(file, "n", H5T_NATIVE_HSIZE, n_attr);
    writeScalarAttribute(file, "m", H5T_NATIVE_HSIZE, m_attr);
    writeScalarAttribute(file, "numValues", H5T_NATIVE_HSIZE, nnz_attr);
    writeScalarAttribute(file, "sym", H5T_NATIVE_INT, sym_attr);
    writeScalarAttribute(file, "hermitian", H5T_NATIVE_INT, herm_attr);
    writeScalarAttribute(file, "upperTri", H5T_NATIVE_INT, upper_attr);

    // --- create datatype for complex values ---
    hid_t complexType = createComplexType();

    // ---------- dataset: values ----------
    {
        hsize_t dims[1] = { static_cast<hsize_t>(this->numValues) };
        hid_t space = H5Screate_simple(1, dims, nullptr);
        hid_t dset = H5Dcreate2(file, "values", complexType, space,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite(dset, complexType, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, valBuf.data());

        H5Dclose(dset);
        H5Sclose(space);
    }

    // ---------- dataset: columns ----------
    {
        hsize_t dims[1] = { static_cast<hsize_t>(this->numValues) };
        hid_t space = H5Screate_simple(1, dims, nullptr);
        hid_t dset = H5Dcreate2(file, "columns", H5T_NATIVE_HSIZE, space,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite(dset, H5T_NATIVE_HSIZE, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, colBuf.data());

        H5Dclose(dset);
        H5Sclose(space);
    }

    // ---------- dataset: rowIndex ----------
    {
        hsize_t dims[1] = { static_cast<hsize_t>(this->numValues) };
        hid_t space = H5Screate_simple(1, dims, nullptr);
        hid_t dset = H5Dcreate2(file, "rowIndex", H5T_NATIVE_HSIZE, space,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite(dset, H5T_NATIVE_HSIZE, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, rowBuf.data());

        H5Dclose(dset);
        H5Sclose(space);
    }

    H5Tclose(complexType);
    H5Fclose(file);
}

void smatrix::loadHDF5(const std::string& filename)
{
    // ---- open file ----
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        std::cerr << "Error: cannot open HDF5 file " << filename << std::endl;
        return;
    }

    // ---- read attributes into temporaries ----
    hsize_t n_attr, m_attr, nnz_attr;
    int sym_attr, herm_attr, upper_attr;

    readScalarAttribute(file, "n", H5T_NATIVE_HSIZE, n_attr);
    readScalarAttribute(file, "m", H5T_NATIVE_HSIZE, m_attr);
    readScalarAttribute(file, "numValues", H5T_NATIVE_HSIZE, nnz_attr);
    readScalarAttribute(file, "sym", H5T_NATIVE_INT, sym_attr);
    readScalarAttribute(file, "hermitian", H5T_NATIVE_INT, herm_attr);
    readScalarAttribute(file, "upperTri", H5T_NATIVE_INT, upper_attr);

    // ---- free old contents if this smatrix already had something ----
    if (this->values)   delete[] this->values;
    if (this->columns)  delete[] this->columns;
    if (this->rowIndex) delete[] this->rowIndex;

    // ---- assign metadata ----
    this->n = static_cast<size_t>(n_attr);
    this->m = static_cast<size_t>(m_attr);
    this->numValues = static_cast<size_t>(nnz_attr);
    this->sym = (sym_attr != 0);
    this->hermitian = (herm_attr != 0);
    this->upperTri = (upper_attr != 0);

    if (this->initialized) {
        if (MKLSparseMatrix != nullptr) {
            mkl_sparse_destroy(*MKLSparseMatrix);
            delete MKLSparseMatrix;
        }
        this->initialized = false;
    }


    // ---- allocate new arrays ----
    this->values = new std::complex<double>[this->numValues];
    this->columns = new size_t[this->numValues];
    this->rowIndex = new size_t[this->numValues];

    // ---- create complex datatype ----
    hid_t complexType = createComplexType();

    // ---------- read: values ----------
    {
        hid_t dset = H5Dopen2(file, "values", H5P_DEFAULT);
        hid_t space = H5Dget_space(dset);

        std::vector<hdf5_complex_t> tmp(this->numValues);

        H5Dread(dset, complexType, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, tmp.data());

        H5Dclose(dset);
        H5Sclose(space);

        // convert back to std::complex<double>
        for (size_t i = 0; i < this->numValues; ++i) {
            this->values[i] = std::complex<double>(tmp[i].real, tmp[i].imag);
        }
    }

    // ---------- read: columns ----------
    {
        hid_t dset = H5Dopen2(file, "columns", H5P_DEFAULT);
        hid_t space = H5Dget_space(dset);

        std::vector<hsize_t> tmp(this->numValues);

        H5Dread(dset, H5T_NATIVE_HSIZE, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, tmp.data());

        H5Dclose(dset);
        H5Sclose(space);

        for (size_t i = 0; i < this->numValues; ++i)
            this->columns[i] = static_cast<size_t>(tmp[i]);
    }

    // ---------- read: rowIndex ----------
    {
        hid_t dset = H5Dopen2(file, "rowIndex", H5P_DEFAULT);
        hid_t space = H5Dget_space(dset);

        std::vector<hsize_t> tmp(this->numValues);

        H5Dread(dset, H5T_NATIVE_HSIZE, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, tmp.data());

        H5Dclose(dset);
        H5Sclose(space);

        for (size_t i = 0; i < this->numValues; ++i)
            this->rowIndex[i] = static_cast<size_t>(tmp[i]);
    }

    H5Tclose(complexType);
    H5Fclose(file);
    initialize();
}

#endif