#pragma once

#ifdef USE_HDF

#include <algorithm>
#include <complex>
#include <string>
#include <vector>

#include <H5Cpp.h>

#include "krylovExceptions.h"
#include "version.h"

/**
* Shared helpers for the HDF5 output of the library. Everything HDF5 specific
* is collected here so that the rest of the code only deals with plain types.
*/
namespace TE {
    namespace hdf5 {

        /**
        * Layout used to store a complex number. The member names "r" and "i" are
        * the convention h5py uses, so a dataset written with this type is read
        * back as a native complex array without any work on the reader side.
        */
        struct complexType {
            double r;
            double i;
        };

        /**
        * Compound datatype matching complexType.
        */
        inline H5::CompType complexDataType()
        {
            H5::CompType type(sizeof(complexType));
            type.insertMember("r", HOFFSET(complexType, r), H5::PredType::IEEE_F64LE);
            type.insertMember("i", HOFFSET(complexType, i), H5::PredType::IEEE_F64LE);
            return type;
        }

        /**
        * HDF5 has no native boolean, so a two valued enumeration is used. This is
        * the representation h5py writes and reads as a real Python bool.
        */
        inline H5::EnumType boolDataType()
        {
            H5::EnumType type(sizeof(unsigned char));
            unsigned char value = 0;
            type.insert("FALSE", &value);
            value = 1;
            type.insert("TRUE", &value);
            return type;
        }

        /**
        * Datasets below this many bytes are written uncompressed. Compressing
        * small datasets makes files larger rather than smaller, because the
        * chunk index and filter overhead outweigh anything deflate can save.
        */
        constexpr hsize_t compressionThreshold = 1024 * 1024;

        /**
        * Deflate level. Higher levels cost noticeably more time without
        * compressing numeric data much further.
        */
        constexpr int compressionLevel = 4;

        /**
        * Creation properties for a one dimensional dataset. Sufficiently large
        * datasets are chunked and compressed; the shuffle filter is applied
        * first because it substantially improves deflate on numeric data.
        */
        inline H5::DSetCreatPropList datasetProperties(hsize_t numElements, size_t elementSize)
        {
            H5::DSetCreatPropList properties;

            if (numElements * elementSize < compressionThreshold || numElements == 0)
                return properties;

            //Aim for chunks of about 1 MB, but never more than the dataset itself.
            hsize_t chunk = compressionThreshold / elementSize;
            if (chunk > numElements)
                chunk = numElements;
            if (chunk == 0)
                chunk = 1;

            properties.setChunk(1, &chunk);
            properties.setShuffle();
            properties.setDeflate(compressionLevel);
            return properties;
        }

        /**
        * Creation properties for a two dimensional dataset that is filled one
        * row at a time and grows as it goes.
        *
        * A dataset can only be extended if it is chunked. A chunk spans whole
        * rows while a row is small and a limited number of columns once it is
        * not, because a chunk that does not fit into the cache is decompressed
        * and rewritten on every row that touches it.
        *
        * @param expectedRows Number of rows the caller intends to write, used for sizing only
        * @param columns Fixed length of a row
        * @param elementSize Size of a single value in the file
        */
        inline H5::DSetCreatPropList extendibleDatasetProperties(hsize_t expectedRows, hsize_t columns,
            size_t elementSize)
        {
            H5::DSetCreatPropList properties;

            hsize_t chunk[2];
            chunk[1] = std::min<hsize_t>(columns, compressionThreshold / elementSize);
            if (chunk[1] == 0)
                chunk[1] = 1;
            chunk[0] = compressionThreshold / (chunk[1] * elementSize);
            if (chunk[0] == 0)
                chunk[0] = 1;
            if (expectedRows != 0 && chunk[0] > expectedRows)
                chunk[0] = expectedRows;

            properties.setChunk(2, chunk);

            if (expectedRows * columns * elementSize >= compressionThreshold)
            {
                properties.setShuffle();
                properties.setDeflate(compressionLevel);
            }

            return properties;
        }

        /**
        * Access properties to go with extendibleDatasetProperties. The default
        * chunk cache holds a single chunk, which makes a partly filled chunk be
        * evicted and recompressed on nearly every row; a few chunks of room
        * turn that back into one compression per chunk.
        */
        inline H5::DSetAccPropList extendibleDatasetAccess()
        {
            H5::DSetAccPropList access;
            access.setChunkCache(521, 8 * compressionThreshold, 0.75);
            return access;
        }

        /**
        * Write a one dimensional dataset, compressed if it is large enough.
        * @param location File or group to write into
        * @param name Name of the new dataset
        * @param fileType Type used to store the data in the file
        * @param memType Type describing the data in memory
        * @param data Pointer to numElements values
        * @param numElements Number of values
        * @param elementSize Size of a single value in the file
        */
        inline void writeDataset(const H5::H5Location& location, const std::string& name,
            const H5::DataType& fileType, const H5::DataType& memType,
            const void* data, hsize_t numElements, size_t elementSize)
        {
            H5::DataSpace space(1, &numElements);
            H5::DataSet dataset = location.createDataSet(name, fileType, space,
                datasetProperties(numElements, elementSize));
            if (numElements != 0)
                dataset.write(data, memType);
        }

        /**
        * Attach a scalar attribute. Scalar rather than an array of length one, so
        * that a reader gets a plain number instead of a one element array.
        */
        template <typename T>
        inline void writeAttribute(const H5::H5Object& object, const std::string& name,
            const H5::DataType& fileType, const H5::DataType& memType, const T& value)
        {
            H5::DataSpace scalar(H5S_SCALAR);
            H5::Attribute attribute = object.createAttribute(name, fileType, scalar);
            attribute.write(memType, &value);
        }

        inline void writeAttribute(const H5::H5Object& object, const std::string& name, double value)
        {
            writeAttribute(object, name, H5::PredType::IEEE_F64LE, H5::PredType::NATIVE_DOUBLE, value);
        }

        inline void writeAttribute(const H5::H5Object& object, const std::string& name, int value)
        {
            writeAttribute(object, name, H5::PredType::STD_I32LE, H5::PredType::NATIVE_INT, value);
        }

        inline void writeAttribute(const H5::H5Object& object, const std::string& name, size_t value)
        {
            hsize_t stored = static_cast<hsize_t>(value);
            writeAttribute(object, name, H5::PredType::STD_U64LE, H5::PredType::NATIVE_HSIZE, stored);
        }

        inline void writeAttribute(const H5::H5Object& object, const std::string& name, bool value)
        {
            H5::EnumType type = boolDataType();
            unsigned char stored = value ? 1 : 0;
            writeAttribute(object, name, type, type, stored);
        }

        inline void writeAttribute(const H5::H5Object& object, const std::string& name, const std::string& value)
        {
            H5::StrType type(H5::PredType::C_S1, H5T_VARIABLE);
            H5::DataSpace scalar(H5S_SCALAR);
            H5::Attribute attribute = object.createAttribute(name, type, scalar);
            attribute.write(type, value);
        }

        /**
        * Attach a scalar attribute, replacing one of the same name if it is
        * already there.
        */
        template <typename T>
        inline void writeOrReplaceAttribute(const H5::H5Object& object, const std::string& name, const T& value)
        {
            if (object.attrExists(name))
                object.removeAttr(name);
            writeAttribute(object, name, value);
        }

        /**
        * Read a scalar attribute.
        */
        template <typename T>
        inline void readAttribute(const H5::H5Object& object, const std::string& name,
            const H5::DataType& memType, T& out)
        {
            H5::Attribute attribute = object.openAttribute(name);
            attribute.read(memType, &out);
        }

        inline size_t readSizeAttribute(const H5::H5Object& object, const std::string& name)
        {
            hsize_t value = 0;
            readAttribute(object, name, H5::PredType::NATIVE_HSIZE, value);
            return static_cast<size_t>(value);
        }

        inline bool readBoolAttribute(const H5::H5Object& object, const std::string& name)
        {
            unsigned char value = 0;
            H5::EnumType type = boolDataType();
            readAttribute(object, name, type, value);
            return value != 0;
        }

        /**
        * Record which version of the library produced the file, so that its
        * structure can be looked up later.
        */
        inline void writeVersion(const H5::H5Object& object)
        {
            writeAttribute(object, "timeEvolverVersion", std::string(TIMEEVOLVER_VERSION_STRING));
            writeAttribute(object, "timeEvolverVersionNumber", (int)TIMEEVOLVER_VERSION);
        }

    }
}

/**
* H5::Exception does not derive from std::exception, so it would slip past every
* handler in the library and in user code. Every entry point that touches HDF5
* wraps its body in this macro to translate it.
*/
#define TE_HDF5_TRY H5::Exception::dontPrint(); try {

#define TE_HDF5_CATCH(context)                                                    \
    } catch (const H5::Exception& e) {                                            \
        throw TE::krylovIOError(std::string(context) + ": " + e.getDetailMsg());  \
    }

#endif
