#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>

int main()
{

    std::system("./Example/simpleExample");

    const int numFiles = 2;

    std::string reference[numFiles] = {"../output/SimpleExampleOutputOccupationNumber0.csv", "../output/SimpleExampleOutputOccupationNumber1.csv"};
    std::string testData[numFiles] = {"SimpleExampleOutputOccupationNumber0.csv", "SimpleExampleOutputOccupationNumber1.csv"};


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
        std::system(("rm " + testData[i]).c_str());
    return 0;

}

