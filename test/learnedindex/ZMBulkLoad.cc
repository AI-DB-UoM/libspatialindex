/******************************************************************************
 * Project:  libspatialindex - A C++ library for spatial indexing
 * Author:   Guanli Liu, liuguanli22@gmail.com
 ******************************************************************************
 * Copyright (c) 2024, Guanli Liu
 *
 * All rights reserved.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
******************************************************************************/

// NOTE: Please read README.txt before browsing this code.
#include <cstring>

// include library header file.
#include <spatialindex/SpatialIndex.h>

using namespace SpatialIndex;

#define INSERT 1
#define DELETE 0
#define QUERY 2

class MyDataStream : public IDataStream
{
public:
	MyDataStream(std::string inputFile) : m_pNext(nullptr)
	{
		m_fin.open(inputFile.c_str());

		if (! m_fin)
			throw Tools::IllegalArgumentException("Input file not found.");

		readNextEntry();
	}

	~MyDataStream() override
	{
		if (m_pNext != nullptr) delete m_pNext;
	}

	IData* getNext() override
	{
		if (m_pNext == nullptr) return nullptr;

		LearnedIndex::Data* ret = m_pNext;
		m_pNext = nullptr;
		readNextEntry();
		return ret;
	}

	bool hasNext() override
	{
		return (m_pNext != nullptr);
	}

	uint32_t size() override
	{
		throw Tools::NotSupportedException("Operation not supported.");
	}

	void rewind() override
	{
		if (m_pNext != nullptr)
		{
			delete m_pNext;
			m_pNext = nullptr;
		}

		m_fin.seekg(0, std::ios::beg);
		readNextEntry();
	}

	void readNextEntry()
	{
		id_type id;
		uint32_t op;
		uint64_t z_value;
		double low[2], high[2];

		m_fin >> op >> id >> low[0] >> low[1] >> high[0] >> high[1] >> z_value;

		if (m_fin.good())
		{
			if (op != INSERT && op != QUERY)
				throw Tools::IllegalArgumentException(
					"The data input should contain insertions only."
				);

			Region r(low, high, 2);

			uint8_t z_value_bytes[sizeof(uint64_t)];
        	memcpy(z_value_bytes, &z_value, sizeof(uint64_t));

			m_pNext = new LearnedIndex::Data(sizeof(uint64_t), z_value_bytes, r, id);
				// Associate a bogus data array with every entry for testing purposes.
				// Once the data array is given to KDTree:Data a local copy will be created.
				// Hence, the input data array can be deleted after this operation if not
				// needed anymore.
		}
	}

	std::ifstream m_fin;
	LearnedIndex::Data* m_pNext;
};
(z_order_output_default):

            logger.info(f"{z_order_output_default} NOT exists")
            transform_command = f"python tools/rank_space_z.py {ablosute_data_file_name} {RANK_SPACE_Z_ORDER_OUTPUT} {bit_num}"
            elapsed_time_ns_order = execute_command(transform_command)

            query_adapt_command = f"python tools/libspatialindex_zm_query_adapter.py --bits {bit_num} --data {ablosute_data_file_name} --query_list {query_list_str}"
            execute_command(query_adapt_command)

            format_data_command = f"python tools/libspatialindex_data_adapter.py --type data --input {RANK_SPACE_Z_ORDER_OUTPUT} --output {data_file}"
            execute_command(format_data_command)

            copy_and_rename(data_file, z_order_output_default)

        else:
            copy_and_rename(z_order_output_default, data_file)

            elapsed_time_ns_order = 0

        logger.info("build zm")


        build_output_path = ZM_BUILD_OUTPUT_PATH.format(
                        data_file_prefix=data_file_prefix,
                        bit_num=bit_num,
                    )

        command = f"test-learnedindex-ZMBulkLoad {data_file} {INDEX_PATH}/zm {page_size} {fill_factor} {BLOCK_SIZE} {BUFFER}"
        result, elapsed_time_ns_build = execute_command_with_err(command)

    
        os.makedirs(os.path.dirname(build_output_path), exist_ok=True)

        with open(build_output_path, "w") as f:
            if result:
                f.write(result.stderr)
            f.write(f"Elapsed Learn Time: {elapsed_time_ns_order}\n")
            f.write(f"Elapsed Build Time: {elapsed_time_ns_build}\n")
            f.write(f"Tree.dat Size: {os.path.getsize(f'{INDEX_PATH}/zm.dat')}\n")
            f.write(f"Tree.idx Size: {os.pat
		}

		std::string baseName = argv[2];
		double utilization = 1.0;
		utilization = atof(argv[4]);
		int pagesize = atoi(argv[5]);
		int mainMemoryRandomBuffer = atoi(argv[6]);

		// if (strcmp(argv[1], "qdtree") == 0)
		// {
		// 	myVariant = SpatialIndex::KDTree::QD_NORMAL;
		// 	loadMethod = SpatialIndex::KDTree::LOAD_QD;
		// }

		IStorageManager* diskfile = StorageManager::createNewDiskStorageManager(baseName, pagesize);
			// Create a new storage manager with the provided base name and a 4K page size.

		StorageManager::IBuffer* file = StorageManager::createNewRandomEvictionsBuffer(*diskfile, mainMemoryRandomBuffer, true);
			// applies a main memory random buffer on top of the persistent storage manager
			// (LRU buffer, etc can be created the same way).

		MyDataStream stream(argv[1]);

		id_type indexIdentifier;
		ISpatialIndex* tree = LearnedIndex::createAndBulkLoadNewLearnedIndex(
			SpatialIndex::LearnedIndex::LOAD_ZM, stream, *file, utilization, atoi(argv[3]), atoi(argv[3]), 2, 
			SpatialIndex::LearnedIndex::LI_ZM, indexIdentifier);
	
		
		std::cerr << *tree;
		std::cerr << "Buffer hits: " << file->getHits() << std::endl;
		std::cerr << "Index ID: " << indexIdentifier << std::endl;

		bool ret = tree->isIndexValid();
		if (ret == false) std::cerr << "ERROR: Structure is invalid!" << std::endl;
		else std::cerr << "The stucture seems O.K." << std::endl;

		delete tree;
		delete file;
		delete diskfile;
			// delete the buffer first, then the storage manager
			// (otherwise the the buffer will fail trying to write the dirty entries).
	}
	catch (Tools::Exception& e)
	{
		std::cerr << "******ERROR******" << std::endl;
		std::string s = e.what();
		std::cerr << s << std::endl;
		return -1;
	}

	return 0;
}
