#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <vector>
#include <sstream>

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

    for(int i = 0; i != 10; i++)
    {
        DataSet dsRef = fileRef.openDataSet("mode"+std::to_string(i));
        dsRef.read(data0, PredType::NATIVE_DOUBLE);

        DataSet dsTest = fileTest.openDataSet("mode"+std::to_string(i));
        dsTest.read(data1, PredType::NATIVE_DOUBLE);

        for(int j = 0; j != 1001; j++)
        {
            double diff = std::abs(data1[j]-data0[j]);
            if (diff > 1e-8)
                return 1; //return != 0 indicated test failure 
        }
    }
    delete[] data0;
    delete[] data1;
    std::system("rm ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0.h5");

#else


    const int numFiles = 10;

    std::string reference[numFiles] = {"../output/ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode0.csv", "../output/ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode1.csv",
    "../output/ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode2.csv","../output/ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode3.csv",
    "../output/ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode4.csv","../output/ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode5.csv",
    "../output/ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode6.csv","../output/ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode7.csv",
    "../output/ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode8.csv","../output/ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode9.csv"};
    std::string testData[numFiles] = {"ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode0.csv", "ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode1.csv",
    "ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode2.csv","ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode3.csv",
    "ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode4.csv","ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode5.csv",
    "ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode6.csv","ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode7.csv",
    "ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode8.csv","ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode9.csv"};

    std::ifstream finRef;
    std::ifstream finOut;

    std::vector<double> valRef;
    std::vector<double> valOut;


    for(int i = 0; i != numFiles; i++)
    {
        finRef.open(reference[i], std::ifstream::ate);
        finOut.open(testData[i], std::ifstream::ate);

        if(finRef.fail() || finOut.fail())
        {
            return 1;
        }

        if (finRef.tellg() != finOut.tellg()) 
        {
            return 1; 
        }

        finRef.seekg(0, std::ifstream::beg);
        finOut.seekg(0, std::ifstream::beg);

        std::string line;

        while(getline(finRef, line))
        {
            std::stringstream sep(line);
            std::string va;

            while (std::getline(sep, va, ','))
            {
                valRef.push_back(std::stod(va));
            }
        }

        while(getline(finOut, line))
        {
            std::stringstream sep(line);
            std::string va;

            while (std::getline(sep, va, ','))
            {
                valOut.push_back(std::stod(va));
            }
        }

        finOut.close();
        finRef.close();

        for(int i = 0; i != valOut.size(); i++)
        {
            double diff = std::abs(valOut[i]-valRef[i]);
            if (diff > 1e-8)
                return 1; //return != 0 indicated test failure 
        }
    }

    for (int i = 0; i != numFiles; i++)
    {          
        std::system(("rm " + testData[i]).c_str());
    }

#endif

    return 0;
}