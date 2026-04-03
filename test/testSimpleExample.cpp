#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "TestUtils.h"

TEST(SimpleExampleRegression, OutputMatchesReference)
{
    test_utils::run_command("./Example/simpleExample");

    const std::vector<std::string> reference_files = {
        "../output/SimpleExampleOutputOccupationNumber0.csv",
        "../output/SimpleExampleOutputOccupationNumber1.csv"
    };
    const std::vector<std::string> output_files = {
        "SimpleExampleOutputOccupationNumber0.csv",
        "SimpleExampleOutputOccupationNumber1.csv"
    };

    test_utils::cleanup_guard cleanup(output_files);

    ASSERT_EQ(reference_files.size(), output_files.size());
    for (std::size_t i = 0; i < reference_files.size(); ++i) {
        SCOPED_TRACE(output_files[i]);
        test_utils::expect_csv_file_near(reference_files[i], output_files[i], 1.0e-8);
    }
}
