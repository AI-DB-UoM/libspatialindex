/******************************************************************************
 * Project:  libspatialindex - A C++ library for spatial indexing
 * Author:   Marios Hadjieleftheriou, mhadji@gmail.com
 ******************************************************************************
 * Copyright (c) 2002, Marios Hadjieleftheriou
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

#include <cstring>
#include <cstdio>
#include <cmath>

#ifndef _MSC_VER
#include <unistd.h>
#endif

#include <spatialindex/SpatialIndex.h>

#include "LearnedIndex.h"
#include "Leaf.h"
#include "Index.h"
#include "BulkLoader.h"

using namespace SpatialIndex;
using namespace SpatialIndex::LearnedIndex;

//
// ExternalSorter::Record
//
ExternalSorter::Record::Record()
= default;

ExternalSorter::Record::Record(const Region& r, id_type id, uint32_t len, uint8_t* pData, uint32_t s)
: m_r(r), m_id(id), m_len(len), m_pData(pData), m_s(s)
{
}

ExternalSorter::Record::~Record()
{
	delete[] m_pData;
}

bool ExternalSorter::Record::operator<(const Record& r) const
{
	int result = std::memcmp(m_pData, r.m_pData, r.m_len);
	return result < 0;
	// if (m_s != r.m_s)
	// 	throw Tools::IllegalStateException("ExternalSorter::Record::operator<: Incompatible sorting dimensions.");

	// if (m_r.m_pHigh[m_s] + m_r.m_pLow[m_s] < r.m_r.m_pHigh[m_s] + r.m_r.m_pLow[m_s])
	// 	return true;
	// else
	// 	return false;
}

void ExternalSorter::Record::storeToFile(Tools::TemporaryFile& f)
{
	f.write(static_cast<uint64_t>(m_id));
	f.write(m_r.m_dimension);
	f.write(m_s);

	for (uint32_t i = 0; i < m_r.m_dimension; ++i)
	{
		f.write(m_r.m_pLow[i]);
		f.write(m_r.m_pHigh[i]);
	}

	f.write(m_len);
	if (m_len > 0) f.write(m_len, m_pData);
}

void ExternalSorter::Record::loadFromFile(Tools::TemporaryFile& f)
{
	m_id = static_cast<id_type>(f.readUInt64());
	uint32_t dim = f.readUInt32();
	m_s = f.readUInt32();

	if (dim != m_r.m_dimension)
	{
		delete[] m_r.m_pLow;
		delete[] m_r.m_pHigh;
		m_r.m_dimension = dim;
		m_r.m_pLow = new double[dim];
		m_r.m_pHigh = new double[dim];
	}

	for (uint32_t i = 0; i < m_r.m_dimension; ++i)
	{
		m_r.m_pLow[i] = f.readDouble();
		m_r.m_pHigh[i] = f.readDouble();
	}

	m_len = f.readUInt32();
	delete[] m_pData; m_pData = nullptr;
	if (m_len > 0) f.readBytes(m_len, &m_pData);
}

//
// ExternalSorter
//
ExternalSorter::ExternalSorter(uint32_t u32PageSize, uint32_t u32BufferPages)
: m_bInsertionPhase(true), m_u32PageSize(u32PageSize),
  m_u32BufferPages(u32BufferPages), m_u64TotalEntries(0), m_stI(0)
{
}

ExternalSorter::~ExternalSorter()
{
	for (m_stI = 0; m_stI < m_buffer.size(); ++m_stI) delete m_buffer[m_stI];
}

void ExternalSorter::insert(Record* r)
{
	if (m_bInsertionPhase == false)
		throw Tools::IllegalStateException("ExternalSorter::insert: Input has already been sorted.");

	m_buffer.push_back(r);
	++m_u64TotalEntries;

	// this will create the initial, sorted buckets before the
	// external merge sort.
	if (m_buffer.size() >= m_u32PageSize * m_u32BufferPages)
	{
		std::sort(m_buffer.begin(), m_buffer.end(), Record::SortAscending());
		Tools::TemporaryFile* tf = new Tools::TemporaryFile();
		for (size_t j = 0; j < m_buffer.size(); ++j)
		{
			m_buffer[j]->storeToFile(*tf);
			delete m_buffer[j];
		}
		m_buffer.clear();
		tf->rewindForReading();
		m_runs.push_back(std::shared_ptr<Tools::TemporaryFile>(tf));
	}
}

void ExternalSorter::sort()
{
	if (m_bInsertionPhase == false)
		throw Tools::IllegalStateException("ExternalSorter::sort: Input has already been sorted.");

	if (m_runs.empty())
	{
		// The data fits in main memory. No need to store to disk.
		std::sort(m_buffer.begin(), m_buffer.end(), Record::SortAscending());
		m_bInsertionPhase = false;
		return;
	}

	if (m_buffer.size() > 0)
	{
		// Whatever remained in the buffer (if not filled) needs to be stored
		// as the final bucket.
		std::sort(m_buffer.begin(), m_buffer.end(), Record::SortAscending());
		Tools::TemporaryFile* tf = new Tools::TemporaryFile();
		for (size_t j = 0; j < m_buffer.size(); ++j)
		{
			m_buffer[j]->storeToFile(*tf);
			delete m_buffer[j];
		}
		m_buffer.clear();
		tf->rewindForReading();
		m_runs.push_back(std::shared_ptr<Tools::TemporaryFile>(tf));
	}

	if (m_runs.size() == 1)
	{
		m_sortedFile = m_runs.front();
	}
	else
	{
		Record* r = nullptr;

		while (m_runs.size() > 1)
		{
            std::shared_ptr<Tools::TemporaryFile> tf(new Tools::TemporaryFile());
			std::vector<std::shared_ptr<Tools::TemporaryFile> > buckets;
			std::vector<std::queue<Record*> > buffers;
			std::priority_queue<PQEntry, std::vector<PQEntry>, PQEntry::SortAscending> pq;

			// initialize buffers and priority queue.
			std::list<std::shared_ptr<Tools::TemporaryFile> >::iterator it = m_runs.begin();
			for (uint32_t i = 0; i < (std::min)(static_cast<uint32_t>(m_runs.size()), m_u32BufferPages); ++i)
			{
				buckets.push_back(*it);
				buffers.emplace_back();

				r = new Record();
				r->loadFromFile(**it);
					// a run cannot be empty initially, so this should never fail.
				pq.push(PQEntry(r, i));

				for (uint32_t j = 0; j < m_u32PageSize - 1; ++j)
				{
					// fill the buffer with the rest of the page of records.
					try
					{
						r = new Record();
						r->loadFromFile(**it);
						buffers.back().push(r);
					}
					catch (Tools::EndOfStreamException&)
					{
						delete r;
						break;
					}
				}
				++it;
			}

			// exhaust buckets, buffers, and priority queue.
			while (! pq.empty())
			{
				PQEntry e = pq.top(); pq.pop();
				e.m_r->storeToFile(*tf);
				delete e.m_r;

				if (! buckets[e.m_u32Index]->eof() && buffers[e.m_u32Index].empty())
				{
					for (uint32_t j = 0; j < m_u32PageSize; ++j)
					{
						try
						{
							r = new Record();
							r->loadFromFile(*buckets[e.m_u32Index]);
							buffers[e.m_u32Index].push(r);
						}
						catch (Tools::EndOfStreamException&)
						{
							delete r;
							break;
						}
					}
				}

				if (! buffers[e.m_u32Index].empty())
				{
					e.m_r = buffers[e.m_u32Index].front();
					buffers[e.m_u32Index].pop();
					pq.push(e);
				}
			}

			tf->rewindForReading();

			// check if another pass is needed.
			uint32_t u32Count = std::min(static_cast<uint32_t>(m_runs.size()), m_u32BufferPages);
			for (uint32_t i = 0; i < u32Count; ++i)
			{
				m_runs.pop_front();
			}

			if (m_runs.size() == 0)
			{
				m_sortedFile = tf;
				break;
			}
			else
			{
				m_runs.push_back(tf);
			}
		}
	}

	m_bInsertionPhase = false;
}

void ExternalSorter::finishLoad()
{
	m_bInsertionPhase = false;
}

ExternalSorter::Record* ExternalSorter::getNextRecord()
{
	// if (m_bInsertionPhase == true)
	// 	throw Tools::IllegalStateException("ExternalSorter::getNextRecord: Input has not been sorted yet.");

	Record* ret;

	if (m_sortedFile.get() == nullptr)
	{
		if (m_stI < m_buffer.size())
		{
			ret = m_buffer[m_stI];
			m_buffer[m_stI] = nullptr;
			++m_stI;
		}
		else
			throw Tools::EndOfStreamException("");
	}
	else
	{
		ret = new Record();
		ret->loadFromFile(*m_sortedFile);
	}

	return ret;
}

inline uint64_t ExternalSorter::getTotalEntries() const
{
	return m_u64TotalEntries;
}

//
// BulkLoader
//
void BulkLoader::bulkLoadUsingZM(
	SpatialIndex::LearnedIndex::LearnedIndex* pTree,
	IDataStream& stream,
	uint32_t bindex,
	uint32_t bleaf,
	uint32_t pageSize,
	uint32_t numberOfPages
) {
	if (! stream.hasNext())
		throw Tools::IllegalArgumentException(
			"LearnedIndex::BulkLoader::bulkLoadUsingZM: Empty data stream given."
		);

	std::cerr << "bulkLoadUsingZM." << std::endl;

	NodePtr n = pTree->readNode(pTree->m_rootID);
	pTree->deleteNode(n.get());

	#ifndef NDEBUG
	std::cerr << "LearnedIndex::BulkLoader: Sorting data." << std::endl;
	#endif

	std::vector<ExternalSorter::Record*> es;

	while (stream.hasNext())
	{
		Data* d = reinterpret_cast<Data*>(stream.getNext());
		if (d == nullptr)
			throw Tools::IllegalArgumentException(
				"bulkLoadUsingSFC: LearnedIndex bulk load expects SpatialIndex::LearnedIndex::Data entries."
			);

		// uint64_t z_value = 0; // Initialize the uint64_t variable
		// memcpy(&z_value, d->m_pData, sizeof(uint64_t)); // Copy the bytes from m_pData into z_value

		// std::cerr << "d->m_pData (as uint64_t): " << z_value << std::endl;

		es.push_back(new ExternalSorter::Record(d->m_region, d->m_id, d->m_dataLength, d->m_pData, 0));
		d->m_pData = nullptr;
		delete d;
	}
	// es->sort();
	// es->finishLoad();
	// std::sort(es.begin(), es.end(), ExternalSorter::Record::SortAscending());

	pTree->m_stats.m_u64Data = es.size();

	std::cerr << "es.size(): " << es.size() << std::endl;

	pTree->m_stats.m_nodesInLevel.push_back(0);

	// create index levels.
	// uint32_t level = 0;

	std::vector<ExternalSorter::Record *> es2;
	createLevelZM(pTree, es, pTree->m_dimension, bleaf, bindex, 1, es2, pageSize, numberOfPages);

	pTree->storeHeader();
}


// std::vector<double> BulkLoader::calculateLinearRegression(std::vector<ExternalSorter::Record*>& es) {
//     if (es.empty()) {
//         return {0.0, 0.0};  // No data
//     }

//     uint64_t min_x = std::numeric_limits<uint64_t>::max();
//     uint64_t max_x = std::numeric_limits<uint64_t>::min();
//     double min_y = 0.0, max_y = 0.0;

//     for (int i = 0; i < es.size(); ++i) {
//         const auto& record = es[i];
//         uint64_t x;
//         std::memcpy(&x, record->m_pData, sizeof(uint64_t));

//         if (x < min_x) {
//             min_x = x;
//             min_y = static_cast<double>(i);
//         }
//         if (x > max_x) {
//             max_x = x;
//             max_y = static_cast<double>(i);
//         }
//     }

//     if (min_x == max_x) {
//         return {min_y, 0.0};  // All x values are the same, slope is zero
//     }

//     double slope = (max_y - min_y) / static_cast<double>(max_x - min_x);
//     double intercept = min_y - slope * static_cast<double>(min_x);

//     return {slope, intercept};
// }


std::vector<double> BulkLoader::calculateLinearRegression(std::vector<ExternalSorter::Record*>& es, uint32_t bindex) {
    double mean_x = 0.0, mean_y = 0.0, covariance = 0.0, variance = 0.0;
    uint64_t n = 0;  // Count of elements
    int data_size = es.size();
	double scale = 1.0 * bindex / data_size;

    for (int i = 0; i < data_size; ++i) {
        const auto& record = es[i];
        
        uint64_t x;
        std::memcpy(&x, record->m_pData, sizeof(uint64_t));
        double y = static_cast<double>(i * scale);

        ++n;

        // Update means
        double dx = x - mean_x;
        mean_x += dx / static_cast<double>(n);
        mean_y += (y - mean_y) / static_cast<double>(n);

        // Update covariance and variance
        covariance += dx * (y - mean_y);
        variance += dx * (x - mean_x);
    }

    // Handle edge cases
    if (data_size == 0) {
        return {0.0, 0.0};  // No data
    }

    if (data_size == 1) {
        return {mean_y, 0.0};  // Only one data point
    }

    variance /= static_cast<double>(n - 1);
    covariance /= static_cast<double>(n - 1);

    if (variance == 0.0) {
        return {mean_y, 0.0};  // No variation in x
    }

    // Calculate regression coefficients
    double slope = covariance / variance;  // Slope
    double intercept = mean_y - slope * mean_x;  // Intercept

    return {slope, intercept};
}


void BulkLoader::createLevelZM(
	SpatialIndex::LearnedIndex::LearnedIndex* pTree,
	std::vector<ExternalSorter::Record *> es,
	uint32_t dimension,
	uint32_t bleaf,
	uint32_t bindex,
	uint32_t level,
	std::vector<ExternalSorter::Record *>& es2,
	uint32_t pageSize,
	uint32_t numberOfPages
) {
	uint64_t total_entries = es.size();
    uint64_t entries_per_child = total_entries / bindex;
    uint64_t remainder = total_entries % bindex;

	// std::cerr << "createLevelZM->  total_entries:" << total_entries << " entries_per_child:" << entries_per_child << std::endl;

    if (pTree->m_stats.m_nodesInLevel.size() < level) {
        pTree->m_stats.m_nodesInLevel.push_back(0);
    }

	if (total_entries <= bleaf) {
		Node* n = createNode(pTree, es, 0);
        es.clear();
        pTree->writeNode(n);
        es2.push_back(new ExternalSorter::Record(n->m_nodeMBR, n->m_identifier, 0, nullptr, 0));
        pTree->m_rootID = n->m_identifier;

        pTree->m_stats.m_u32TreeHeight = std::max(pTree->m_stats.m_u32TreeHeight, static_cast<uint32_t>(level));

        delete n;
	} else {
		std::vector<std::vector<ExternalSorter::Record *>> children_es(bindex);
        ExternalSorter::Record* pR;


		uint64_t i = 0;
		// entries_per_child = entries_per_child + (remainder > 0 ? 1 : 0);
		std::vector<double> regressionParams = calculateLinearRegression(es, bindex - 1);

		while (i< total_entries) {
			try { pR = es[i++]; }
            catch (Tools::EndOfStreamException&) { break; }
			uint64_t x;
    		std::memcpy(&x, pR->m_pData, sizeof(uint64_t)); 
			double a = regressionParams[0];
    		double b = regressionParams[1];
    		double y = a * x + b;

			uint32_t predict_index = static_cast<uint32_t>(y);
			// std::cerr << "prediction a:" << a << " b:" << b << " x:" << x << " y:" << y << " predict_index:" << predict_index << std::endl;

			predict_index = std::max(static_cast<uint32_t>(0), std::min(predict_index, static_cast<uint32_t>(bindex - 1)));
			// uint64_t predict_index = i / entries_per_child;
			children_es[predict_index].push_back(pR);
		}

		// recursively prediction
        std::vector<ExternalSorter::Record*> child_records;
		/* TODO do not record child_records, only record model, but we need 
		 to execute knn, so we use this to calculate the MBR*/
        for (uint32_t i = 0; i < bindex; ++i) {
            createLevelZM(pTree, children_es[i], dimension, bleaf, bindex, level + 1, es2, pageSize, numberOfPages);
            child_records.push_back(es2.back());
            es2.pop_back();
        }

		Node* n_parent = createNode(pTree, child_records, pTree->m_stats.m_u32TreeHeight - level);
		n_parent->insertModel(regressionParams);

        pTree->writeNode(n_parent);
        es2.push_back(new ExternalSorter::Record(n_parent->m_nodeMBR, n_parent->m_identifier, 0, nullptr, 0));
        pTree->m_rootID = n_parent->m_identifier;
        delete n_parent;
	}
}

Node* BulkLoader::createNode(SpatialIndex::LearnedIndex::LearnedIndex* li, std::vector<ExternalSorter::Record*>& e, uint32_t level)
{
	Node* n;

	if (level == 0) n = new Leaf(li, -1);
	else n = new Index(li, -1, level);

	for (size_t cChild = 0; cChild < e.size(); ++cChild)
	{
		n->insertEntry(e[cChild]->m_len, e[cChild]->m_pData, e[cChild]->m_r, e[cChild]->m_id);
		e[cChild]->m_pData = nullptr;
		delete e[cChild];
	}

	return n;
}