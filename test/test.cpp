#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>

using namespace std;

#ifdef USE_HDF
#include <H5Cpp.h>
using namespace H5;
#endif



int main()
{
    std::system("./Example/main");

#ifdef USE_HDF
    std::string reference = "../output/ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0.h5";
    std::string testData = "ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0.h5";

    H5File fileRef (reference.c_str(), H5F_ACC_RDONLY);
    H5File fileTest (testData.c_str(), H5F_ACC_RDONLY);
    double* data0 = new double[1001];
    double* data1 = new double[1001];

    for(int i = 0; i != 1; i++)
    {
        DataSet dsRef = fileRef.openDataSet("mode"+std::to_string(i));
        dsRef.read(data0, PredType::NATIVE_DOUBLE);

        DataSet dsTest = fileTest.openDataSet("mode"+std::to_string(i));
        dsTest.read(data1, PredType::NATIVE_DOUBLE);

        for(int j = 0; j != 1001; j++)
        {
            double diff = std::abs(data1[j]-data0[j]);
            std::cout << diff << std::endl;
            if (diff > 1e-8)
                return 1; //return != 0 indicated test failure 
        }
    }
    delete[] data0;
    delete[] data1;
#endif

    std::system("rm ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0.h5");
    return 0;
}