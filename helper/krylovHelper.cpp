#include "krylovHelper.h"


/**
 * Writes sampled expectation values of a list of observables to file
 * @param obs_list List of observables with sampled expectation values
 * @param para List of parameters that are included as attributes in the HDF5 file as well as filename (if requested)
 * @param name Requested filename
 */
void saveResult(const std::vector<std::unique_ptr<krylovBasicObservable>>& obs_list, parameter_list& para, const std::string& name)
{
    std::string outputFileName = name;

    parameter_list::iterator paraIter;
    std::vector<std::unique_ptr<krylovBasicObservable>>::const_iterator obsIter;

    //Create string with parameter and its value
    for (paraIter = para.begin(); paraIter != para.end(); paraIter++) {
        if ((*paraIter)->getPrintFilename() == true)
            outputFileName += "_" + (*paraIter)->getName() + (*paraIter)->getData();
    }



    obsIter = obs_list.begin();
    size_t nbOutputParameters = para.size();


#ifdef USE_HDF
    outputFileName += ".h5";
    H5File fileHh(outputFileName.c_str(), H5F_ACC_TRUNC);
    DataSet dataset;

    for (unsigned int j = 0; j < obs_list.size() && obsIter != obs_list.end(); j++)
    {
        int NX = (int)(*obsIter)->retNumSamples();
        const int RANK = 1;
        hsize_t dimsf[RANK] = {};
        dimsf[0] = NX;
        DataSpace dataspace(RANK, dimsf);
        FloatType datatype(PredType::NATIVE_DOUBLE);
        datatype.setOrder(H5T_ORDER_LE);
        std::string observableName = (*obsIter)->retName();
        dataset = fileHh.createDataSet(observableName.c_str(), datatype, dataspace);
        dataset.write((*obsIter)->retExpectationValues(), PredType::NATIVE_DOUBLE);

        //write  parameters as attributes to file
        hsize_t dims[1] = { 1 };

        DataSpace** attr_dataspace = new DataSpace * [nbOutputParameters];
        Attribute** attributes = new Attribute * [nbOutputParameters];

        for (unsigned int i = 0; i < nbOutputParameters; ++i)
        {
            attr_dataspace[i] = new DataSpace(1, dims);
            attributes[i] = new Attribute();
        }

        int counter = 0;
        for (paraIter = para.begin(); paraIter != para.end(); paraIter++, counter++)
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
            }
            else if ((*paraIter)->isBool())
            {
                int paraValue = dynamic_cast<typedParameter<bool>&>(*(*paraIter)).getValue();
                *attributes[counter] = dataset.createAttribute((*paraIter)->getName(), PredType::NATIVE_INT, *attr_dataspace[counter]); attributes[counter]->write(PredType::NATIVE_INT, &paraValue);
            }
        }

        for (unsigned int i = 0; i < nbOutputParameters; ++i)
        {
            delete attr_dataspace[i];
            delete attributes[i];
        }

        delete[] attr_dataspace;
        delete[] attributes;

        dataset.close();
        obsIter++;

    }
    //If not write data to simple csv files
#else

    for (; obsIter != obs_list.end(); obsIter++)
    {
        std::string fileNameCSV = outputFileName + (*obsIter)->retName() + ".csv";
        std::ofstream outputfile;
        outputfile.open(fileNameCSV);


        for (int i = 0; i != (*obsIter)->retNumSamples() - 1; i++)
            outputfile << (*obsIter)->retExpectationValues()[i] << ", ";
        outputfile << (*obsIter)->retExpectationValues()[((*obsIter)->retNumSamples()) - 1];
        outputfile.close();
    }

#endif

}