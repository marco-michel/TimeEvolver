#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "TestUtils.h"

TEST(BlackHoleRegression, OutputMatchesReference)
{
    test_utils::run_command("./Example/main --N0 20 --Nm 2 --K 4 --C0 1 --Cm 1 --maxT 10 --samplingStep 0.01 --tol 1e-08 --m 40 --DeltaN 12 --capacity 1 --fastIntegration 0");

#ifdef USE_HDF
    const std::string reference_path = "../output/ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0.h5";
    const std::string output_path = "ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0.h5";
    test_utils::cleanup_guard cleanup({ output_path });

    for (int mode = 0; mode < 10; ++mode) {
        const std::string dataset_name = "mode" + std::to_string(mode);
        SCOPED_TRACE(dataset_name);
        test_utils::expect_hdf_dataset_near(reference_path, output_path, dataset_name, 1.0e-8);
    }
#else
    const std::vector<std::string> reference_files = {
        "../output/ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode0.csv",
        "../output/ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode1.csv",
        "../output/ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode2.csv",
        "../output/ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode3.csv",
        "../output/ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode4.csv",
        "../output/ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode5.csv",
        "../output/ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode6.csv",
        "../output/ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode7.csv",
        "../output/ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode8.csv",
        "../output/ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode9.csv"
    };
    const std::vector<std::string> output_files = {
        "ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode0.csv",
        "ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode1.csv",
        "ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode2.csv",
        "ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode3.csv",
        "ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode4.csv",
        "ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode5.csv",
        "ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode6.csv",
        "ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode7.csv",
        "ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode8.csv",
        "ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0mode9.csv"
    };
    test_utils::cleanup_guard cleanup(output_files);

    ASSERT_EQ(reference_files.size(), output_files.size());
    for (std::size_t i = 0; i < reference_files.size(); ++i) {
        SCOPED_TRACE(output_files[i]);
        test_utils::expect_csv_file_near(reference_files[i], output_files[i], 1.0e-8);
    }
#endif
}
