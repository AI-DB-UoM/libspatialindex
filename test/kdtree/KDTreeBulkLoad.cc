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

		KDTree::Data* ret = m_pNext;
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
		double low[2], high[2];

		m_fin >> op >> id >> low[0] >> low[1] >> high[0] >> high[1];

		if (m_fin.good())
		{
			if (op != INSERT && op != QUERY)
				throw Tools::IllegalArgumentException(
					"The data input should contain insertions only."
				);

			Region r(low, high, 2);
			m_pNext = new KDTree::Data(sizeof(double), reinterpret_cast<uint8_t*>(low), r, id);
				// Associate a bogus data array with every entry for testing purposes.
				// Once the data array is given to KDTree:Data a local copy will be created.
				// Hence, the input data array can be deleted after this operation if not
				// needed anymore.
		}
	}

	std::ifstream m_fin;
	KDTree::Data* m_pNext;
};

int main(int argc, char** argv)
{
	try
	{
		if (argc != 7)
		{
			std::cerr << "Usage: " << argv[0] << "kd_type input_file query_file tree_file capacity utilization." << std::endl;
			return -1;
		}

		std::string baseName = argv[4];
		double utilization = 1.0;
		utilization = atof(argv[6]);
		SpatialIndex::KDTree::KDTreeVariant myVariant = SpatialIndex::KDTree::KD_NORMAL;
		SpatialIndex::KDTree::LoadMethod loadMethod = SpatialIndex::KDTree::LOAD_KD;

		if (strcmp(argv[1], "kdtree") == 0)
		{
			myVariant = SpatialIndex::KDTree::KD_NORMAL;
			loadMethod = SpatialIndex::KDTree::LOAD_KD;
		}
		if (strcmp(argv[1], "greedy_kdtree") == 0)
		{
			myVariant = SpatialIndex::KDTree::KD_GREEDY;
			loadMethod = SpatialIndex::KDTree::LOAD_KD_GREEDY;
		}
		if (strcmp(argv[1], "qdtree") == 0)
		{
			myVariant = SpatialIndex::KDTree::QD_NORMAL;
			loadMethod = SpatialIndex::KDTree::LOAD_QD;
		}

		IStorageManager* diskfile = StorageManager::createNewDiskStorageManager(baseName, 4096);
			// Create a new storage manager with the provided base name and a 4K page size.

		StorageManager::IBuffer* file = StorageManager::createNewRandomEvictionsBuffer(*diskfile, 10, false);
			// applies a main memory random buffer on top of the persistent storage manager
			// (LRU buffer, etc can be created the same way).

		MyDataStream stream(argv[2]);
		id_type indexIdentifier;
		ISpatialIndex* tree = nullptr;  
		// Create and bulk load a new KDTree with dimensionality 2, using "file" as
		// the StorageManager.
		if (myVariant == SpatialIndex::KDTree::KD_NORMAL)
		{
			std::cerr << "myVariant: " << myVariant << std::endl;
			tree = KDTree::createNewKDTree(
				loadMethod, stream, *file, utilization, 2, atoi(argv[5]), 2, myVariant, indexIdentifier);
		}
		else if (myVariant == SpatialIndex::KDTree::KD_GREEDY || myVariant == SpatialIndex::KDTree::QD_NORMAL)
		{
			MyDataStream queryStream(argv[3]);
			tree = KDTree::createNewKDTree(
				loadMethod, stream, *file, utilization, 2, atoi(argv[5]), 2, myVariant, indexIdentifier, queryStream);
		}
		
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
