#include <vector>
#include <memory>
#include <iomanip>
#include <type_traits>
#include <fstream>
#include <sstream>
#include <iostream>

#include "krylovHelper.h"


#ifdef USE_HDF
#include <H5Cpp.h>
using namespace H5;
#endif

//--------------------------------------------------------------------------
//Class Parameter


std::string parameter::getName()
{
    return p_name;
}


//--------------------------------------------------------------------------
//Class outputHelper

void outputHelper::saveResult()
{

    std::string outputFileName = fileName;

    parameter_list::iterator paraIter;
    observable_list::iterator obsIter;

    //Create string with parameter and its value
    for (paraIter = para_list.begin(); paraIter != para_list.end(); paraIter++)
        outputFileName += "_" + (*paraIter)->getName() + (*paraIter)->getData();

    //Sort data from time ordering to observable ordering
    double** nicelySorted = new double* [nbObservables];
    for (int i = 0; i < nbObservables; ++i)
        nicelySorted[i] = new double[results->nSamples];

    for (int j = 0; j < nbObservables; j++)
    {
        for (unsigned int i = 0; i < results->nSamples; i++)
        {
            nicelySorted[j][i] = (results->sampling->values + nbObservables * i + j)->real();
        }
    }

    obsIter = obs_list.begin();

    //Use HDF output if HDF5 libraries are discovered during compiling
#ifdef USE_HDF

    outputFileName += ".h5";
    H5File fileHh(outputFileName.c_str(), H5F_ACC_TRUNC);


    for (int j = 0; j < nbObservables && obsIter != obs_list.end(); j++)
    {
        int NX = results->nSamples;
        const int RANK = 1;
        hsize_t dimsf[RANK];
        dimsf[0] = NX;
        DataSpace dataspace(RANK, dimsf);
        FloatType datatype(PredType::NATIVE_DOUBLE);
        datatype.setOrder(H5T_ORDER_LE);
        std::string observableName = (*obsIter)->getName();
        dataset = fileHh.createDataSet(observableName.c_str(), datatype, dataspace);
        dataset.write(nicelySorted[j], PredType::NATIVE_DOUBLE);

        writeAttributes();

        dataset.close();
        obsIter++;

    }
    //If not write data to simple csv files
#else

    for (int j = 0; j != nbObservables && obsIter != obs_list.end(); j++)
    {
        std::string fileNameCSV = outputFileName + (*obsIter)->getName() + ".csv";
        std::ofstream outputfile;
        outputfile.open(fileNameCSV);
        for (int i = 0; i != (results->nSamples) - 1; i++)
            outputfile << nicelySorted[j][i] << ", ";
        outputfile << nicelySorted[j][(results->nSamples) - 1];
        outputfile.close();
    }

#endif
    //end write obserables to file


    //clean up
    for (int i = 0; i < nbObservables; i++)
    {
        delete[] nicelySorted[i];
    }
    delete[] nicelySorted; 


}

#ifdef USE_HDF

void outputHelper::writeAttributes() 
{
    //write all parameters as attributes to file
    hsize_t dims[1] = { 1 };

    DataSpace** attr_dataspace = new DataSpace * [nbOutputParameters];
    Attribute** attributes = new Attribute * [nbOutputParameters];

    for (int i = 0; i < nbOutputParameters; ++i)
    {
        attr_dataspace[i] = new DataSpace(1, dims);
        attributes[i] = new Attribute();
    }

    parameter_list::iterator paraIter;
    int counter = 0;
    for (paraIter = para_list.begin(); paraIter != para_list.end(); paraIter++, counter++)
    {
        if ((*paraIter)->isDouble())
        {
    	     double paraValue = dynamic_cast<typedParameter<double>&>(*(*paraIter)).getValue();
            *attributes[counter] = dataset.createAttribute((*paraIter)->getName(), PredType::NATIVE_DOUBLE, *attr_dataspace[counter]); attributes[counter]->write(PredType::NATIVE_DOUBLE, &paraValue);
        }
        else if ((*paraIter)->isInt())
        {
            int paraValue = dynamic_cast<typedParameter<int>&>(*(*paraIter)).getValue();
            *attributes[counter] = dataset.createAttribute((*paraIter)->getName(), PredType::NATIVE_INT, *attr_dataspace[counter]); attributes[counter]->write(PredType::NATIVE_INT, &paraValue);
        } else if ((*paraIter)->isBool()) 
        {
            int paraValue = dynamic_cast<typedParameter<bool>&>(*(*paraIter)).getValue();
            *attributes[counter] = dataset.createAttribute((*paraIter)->getName(), PredType::NATIVE_INT, *attr_dataspace[counter]); attributes[counter]->write(PredType::NATIVE_INT, &paraValue);
        }
    }

    for (int i = 0; i < nbOutputParameters; ++i)
    {
        delete attr_dataspace[i];
        delete attributes[i];
    }

    delete[] attr_dataspace;
    delete[] attributes;
}

#endif
