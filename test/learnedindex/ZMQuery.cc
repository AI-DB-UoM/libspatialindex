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

#include <chrono>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <cmath> // For calculating variance

// include library header file.
#include <spatialindex/SpatialIndex.h>

#include <iomanip> // Include for setprecision

using namespace SpatialIndex;
using namespace std;

#define INSERT 1
#define DELETE 0
#define QUERY 2

// example of a Visitor pattern.
// findes the index and leaf IO for answering the query and prints
// the resulting data IDs to stdout.
class MyVisitor : public IVisitor
{
public:
	size_t m_indexIO{0};
	size_t m_leafIO{0};
	size_t m_dataTime{0};
	size_t m_dataCount{0};
	size_t m_indexTime{0};
	size_t m_leafTime{0};

public:
    MyVisitor() = default;

	void visitNode(const INode& n) override
	{
		if (n.isLeaf()) m_leafIO++;
		else m_indexIO++;
	}

	void visitData(const IData& d) override
	{
		auto start = chrono::high_resolution_clock::now();
		IShape* pS;
		d.getShape(&pS);
		// do something.

		// Region* pr = dynamic_cast<Region*>(pS);

		// cerr << pr->m_pLow[0] << " " << pr->m_pLow[1] << endl;
		// cerr << pr->m_pHigh[0] << " " << pr->m_pHigh[1] << endl << endl;

		delete pS;

		// data should be an array of characters representing a Region as a string.
		uint8_t* pData = nullptr;
		uint32_t cLen = 0;
		d.getData(cLen, &pData);
		// do something.
		//string s = reinterpret_cast<char*>(pData);
		//cout << s << endl;
		delete[] pData;

		// TODO Uncomment the following code if you want to see the query result
		// cout << d.getIdentifier() << endl;
			// the ID of this data entry is an answer to the query. I will just print it to stdout.
		auto end = chrono::high_resolution_clock::now(); 
		auto duration = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
		m_dataTime += duration;
		m_dataCount++;
	}

	void visitData(std::vector<const IData*>& v) override
	{
		cout << v[0]->getIdentifier() << " " << v[1]->getIdentifier() << endl;
	}

	void visitNodeCost(const INode& n, size_t time_cost) override 
	{
		if (n.isLeaf()) m_leafTime += time_cost;
		else m_indexTime += time_cost;
	}
};

// example of a Strategy pattern.
// traverses the tree by level.
class MyQueryStrategy : public SpatialIndex::IQueryStrategy
{
private:
	queue<id_type> ids;

public:
	void getNextEntry(const IEntry& entry, id_type& nextEntry, bool& hasNext) override
	{
		IShape* ps;
		entry.getShape(&ps);
		Region* pr = dynamic_cast<Region*>(ps);

		cout << pr->m_pLow[0] << " " << pr->m_pLow[1] << endl;
		cout << pr->m_pHigh[0] << " " << pr->m_pLow[1] << endl;
		cout << pr->m_pHigh[0] << " " << pr->m_pHigh[1] << endl;
		cout << pr->m_pLow[0] << " " << pr->m_pHigh[1] << endl;
		cout << pr->m_pLow[0] << " " << pr->m_pLow[1] << endl << endl << endl;
			// print node MBRs gnuplot style!

		delete ps;

		const INode* n = dynamic_cast<const INode*>(&entry);

		// traverse only index nodes at levels 2 and higher.
		if (n != nullptr && n->getLevel() > 1)
		{
			for (uint32_t cChild = 0; cChild < n->getChildrenCount(); cChild++)
			{
				ids.push(n->getChildIdentifier(cChild));
			}
		}

		if (! ids.empty())
		{
			nextEntry = ids.front(); ids.pop();
			hasNext = true;
		}
		else
		{
			hasNext = false;
		}
	}
};

// example of a Strategy pattern.
// find the total indexed space managed by the index (the MBR of the root).
class MyQueryStrategy2 : public IQueryStrategy
{
public:
	Region m_indexedSpace;

public:
	void getNextEntry(const IEntry& entry, id_type& /* nextEntry */, bool& hasNext) override
	{
		// the first time we are called, entry points to the root.

		// stop after the root.
		hasNext = false;

		IShape* ps;
		entry.getShape(&ps);
		ps->getMBR(m_indexedSpace);
		delete ps;
	}
};

int main(int argc, char** argv)
{
	try
	{
		if (argc != 5)
		{
			cerr << "Usage: " << argv[0] << " query_file tree_file query_type [intersection | kNN | selfjoin] buffer." << endl;
			return -1;
		}

		uint32_t queryType = 0;
		int k = 0; // To store the value of k for kNN queries
		int len = strlen(argv[3]);
		int mainMemoryRandomBuffer = atof(argv[4]);

		if (strcmp(argv[3], "intersection") == 0) queryType = 0;
		else if (len >= 3 && strcmp(&argv[3][len - 2], "NN") == 0) 
        {
            queryType = 1; // Assuming kNN query type is represented by 1
            char kStr[10];
			strncpy(kStr, argv[3], len - 2);
			kStr[len - 2] = '\0'; 
			k = atoi(kStr); 
			cerr << "knn query: k =" << k << endl;

        }
		else if (strcmp(argv[3], "selfjoin") == 0) queryType = 2;
		else
		{
			cerr << "Unknown query type." << endl;
			return -1;
		}

		ifstream fin(argv[1]);
		if (! fin)
		{
			cerr << "Cannot open query file " << argv[1] << "." << endl;
			return -1;
		}

		string baseName = argv[2];
		IStorageManager* diskfile = StorageManager::loadDiskStorageManager(baseName);
			// this will try to locate and open an already existing storage manager.

		StorageManager::IBuffer* file = StorageManager::createNewRandomEvictionsBuffer(*diskfile, mainMemoryRandomBuffer, true);
			// applies a main memory random buffer on top of the persistent storage manager
			// (LRU buffer, etc can be created the same way).

		// If we need to open an existing tree stored in the storage manager, we only
		// have to specify the index identifier as follows
		ISpatialIndex* tree = LearnedIndex::loadLearnedIndex(*file, 1);

		size_t count = 0;
		size_t indexIO = 0;
		size_t leafIO = 0;
		size_t dataTime = 0;
		size_t indexTime = 0;
		size_t leafTime = 0;
		id_type id;
		uint32_t op;
		double x1, x2, y1, y2;
		double plow[2], phigh[2];
		uint64_t key_low;
		uint64_t key_high;

		vector<double> queryTimes; 
		vector<double> insertTimes; 
		vector<size_t> dataCounts; 

		cerr << fixed << setprecision(10);

		while (fin)
		{
			fin >> op >> id >> x1 >> y1 >> x2 >> y2 >> key_low >> key_high;
			if (! fin.good()) continue; // skip newlines, etc.

			if (op == QUERY)
			{
				auto start = chrono::high_resolution_clock::now();
				plow[0] = x1; plow[1] = y1;
				phigh[0] = x2; phigh[1] = y2;

				MyVisitor vis;

				if (queryType == 0)
				{
					Region r = Region(plow, phigh, 2);
					// tree->intersectsWithQuery(r, vis);
					tree->intersectsWithQueryLearnedIndex(r, vis, key_low, key_high);

						// this will find all data that intersect with the query range.
				}
				else if (queryType == 1)
				{
					Point p = Point(plow, 2);
					tree->nearestNeighborQuery(k, p, vis);
						// this will find the 10 nearest neighbors.
				}
				else
				{
					Region r = Region(plow, phigh, 2);
					tree->selfJoinQuery(r, vis);
				}

				indexIO += vis.m_indexIO;
				leafIO += vis.m_leafIO;
				dataTime += vis.m_dataTime;
				leafTime += vis.m_leafTime;
				indexTime += vis.m_indexTime;
				auto end = chrono::high_resolution_clock::now(); 
				auto duration = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
				queryTimes.push_back(duration); 
				dataCounts.push_back(vis.m_dataCount);
			}
			else if (op == INSERT)
			{
				auto start = chrono::high_resolution_clock::now();
				plow[0] = x1; plow[1] = y1;
				phigh[0] = x2; phigh[1] = y2;
				Region r = Region(plow, phigh, 2);

				std::ostringstream os;
				os << r;
				std::string data = os.str();

				tree->insertData((uint32_t)(data.size() + 1), reinterpret_cast<const uint8_t*>(data.c_str()), r, id);
				//tree->insertData(0, 0, r, id);
					// example of passing zero size and a null pointer as the associated data.
				auto end = chrono::high_resolution_clock::now(); 
				auto duration = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
				insertTimes.push_back(duration); 
			}
			else
			{
				cerr << "This is not a query operation." << endl;
			}

			count++;
		}

		

		MyQueryStrategy2 qs;
		tree->queryStrategy(qs);

		cerr << "Indexed space: " << qs.m_indexedSpace << endl;
		cerr << "Operations: " << count << endl;
		cerr << *tree;
		cerr << "Index I/O: " << indexIO << endl;
		cerr << "Leaf I/O: " << leafIO << endl;
		cerr << "Index Time: " << indexTime << endl;
		cerr << "Leaf Time: " << leafTime << endl;
		cerr << "Data Time: " << dataTime << endl;
		cerr << "Buffer hits: " << file->getHits() << endl;

		size_t totalSum = std::accumulate(dataCounts.begin(), dataCounts.end(), 0ULL);
		cerr << "Total sum of dataCounts: " << totalSum << endl;

		if (queryTimes.size() > 0)
		{
			// Query mean
			double mean = accumulate(queryTimes.begin(), queryTimes.end(), 0.0) / queryTimes.size();

			// Query stdDev
			double variance = 0.0;
			for (auto time : queryTimes) {
				variance += (time - mean) * (time - mean);
			}
			if (queryTimes.size() > 1)
				variance /= (queryTimes.size() - 1);
			double stdDev = sqrt(variance);
	
			cerr << "Query num: " << queryTimes.size() << endl;
			cerr << "Query mean: " << mean << endl;
			cerr << "Query variance: " << variance << endl;
			cerr << "Query stdDev: " << stdDev << endl;
			// cerr << "Query p50: " << p50 << endl;
			// cerr << "Query p99: " << p99 << endl;
			// Query P1 to P100
			size_t n = queryTimes.size();
			for (int p = 0; p < 100; ++p) {
				nth_element(queryTimes.begin(), queryTimes.begin() + n * p / 100, queryTimes.end());
				size_t percentile = queryTimes[n * p / 100];
				size_t dataCount = dataCounts[n * p / 100];
				cerr << "Query P" << p << ": " << percentile << endl;
				cerr << "Query P" << p << " Count: " << dataCount << endl;
			}
		}

		if (insertTimes.size() > 0)
		{
			// Insert mean
			double mean = accumulate(insertTimes.begin(), insertTimes.end(), 0.0) / insertTimes.size();

			// Insert stdDev
			double variance = 0.0;
			for (auto time : insertTimes) {
				variance += (time - mean) * (time - mean);
			}
			if (insertTimes.size() > 1)
				variance /= (insertTimes.size() - 1);
			double stdDev = sqrt(variance);
	
			// Insert P50
			nth_element(insertTimes.begin(), insertTimes.begin() + insertTimes.size() / 2, insertTimes.end());
			double p50 = insertTimes[insertTimes.size() / 2];

			// Insert P99
			size_t n = insertTimes.size();
			nth_element(insertTimes.begin(), insertTimes.begin() + n * 99 / 100, insertTimes.end());
			double p99 = insertTimes[n * 99 / 100];

			cerr << "Insert num: " << insertTimes.size() << endl;
			cerr << "Insert mean: " << mean << endl;
			cerr << "Insert variance: " << variance << endl;
			cerr << "Insert stdDev: " << stdDev << endl;
			cerr << "Insert p50: " << p50 << endl;
			cerr << "Insert p99: " << p99 << endl;
		}		

		delete tree;
		delete file;
		delete diskfile;
			// delete the buffer first, then the storage manager
			// (otherwise the the buffer will fail writting the dirty entries).
	}
	catch (Tools::Exception& e)
	{
		cerr << "******ERROR******" << endl;
		std::string s = e.what();
		cerr << s << endl;
		return -1;
	}
	catch (...)
	{
		cerr << "******ERROR******" << endl;
		cerr << "other exception" << endl;
		return -1;
	}

	return 0;
}
