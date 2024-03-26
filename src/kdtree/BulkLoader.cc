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

#include <cstring>
#include <cstdio>
#include <cmath>

#ifndef _MSC_VER
#include <unistd.h>
#endif

#include <spatialindex/SpatialIndex.h>

#include "KDTree.h"
#include "Leaf.h"
#include "Index.h"
#include "BulkLoader.h"

using namespace SpatialIndex;
using namespace SpatialIndex::KDTree;

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
	if (m_s != r.m_s)
		throw Tools::IllegalStateException("ExternalSorter::Record::operator<: Incompatible sorting dimensions.");

	if (m_r.m_pHigh[m_s] + m_r.m_pLow[m_s] < r.m_r.m_pHigh[m_s] + r.m_r.m_pLow[m_s])
		return true;
	else
		return false;
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
	// if (m_bInsertionPhase == false)
	// 	throw Tools::IllegalStateException("ExternalSorter::insert: Input has already been sorted.");

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
	// if (m_bInsertionPhase == false)
	// 	throw Tools::IllegalStateException("ExternalSorter::sort: Input has already been sorted.");

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
void BulkLoader::topDownPartitioning(
	SpatialIndex::KDTree::KDTree* pTree,
	IDataStream& stream,
	uint32_t bindex,
	uint32_t bleaf,
	uint32_t pageSize,
	uint32_t numberOfPages
) {
	if (! stream.hasNext())
		throw Tools::IllegalArgumentException(
			"KDTree::BulkLoader::topDownPartitioning: Empty data stream given."
		);

	NodePtr n = pTree->readNode(pTree->m_rootID);
	pTree->deleteNode(n.get());

	#ifndef NDEBUG
	std::cerr << "KDTree::BulkLoader: Sorting data." << std::endl;
	#endif

    std::shared_ptr<ExternalSorter> es = std::shared_ptr<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));

	while (stream.hasNext())
	{
		Data* d = reinterpret_cast<Data*>(stream.getNext());
		if (d == nullptr)
			throw Tools::IllegalArgumentException(
				"topDownPartitioning: KDTree bulk load expects SpatialIndex::KDTree::Data entries."
			);

		es->insert(new ExternalSorter::Record(d->m_region, d->m_id, d->m_dataLength, d->m_pData, 0));
		d->m_pData = nullptr;
		delete d;
	}
	es->sort();

	pTree->m_stats.m_u64Data = es->getTotalEntries();

	pTree->m_stats.m_nodesInLevel.push_back(0);

	std::shared_ptr<ExternalSorter> es2 = std::shared_ptr<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));
	partition(pTree, es, pTree->m_dimension, bleaf, bindex, 1, es2, pageSize, numberOfPages);

	pTree->storeHeader();
}


void BulkLoader::partition(
	SpatialIndex::KDTree::KDTree* pTree,
	std::shared_ptr<ExternalSorter> es,
	uint32_t dimension,
	uint32_t bleaf,
	uint32_t bindex,
	uint32_t level,
	std::shared_ptr<ExternalSorter> es2,
	uint32_t pageSize,
	uint32_t numberOfPages
) {

	uint64_t total_entries = es->getTotalEntries();
	uint64_t left_node_en = static_cast<uint64_t>(total_entries / 2);
	uint64_t right_node_en = static_cast<uint64_t>(total_entries - left_node_en);

	if (pTree->m_stats.m_nodesInLevel.size() < level)
	{
		pTree->m_stats.m_nodesInLevel.push_back(0);
	}

	if (total_entries <= bleaf)
	{
		std::vector<ExternalSorter::Record*> node;
		ExternalSorter::Record* r;
		uint32_t i = 0;

		while (i < total_entries)
		{
			try { r = es->getNextRecord(); } catch (Tools::EndOfStreamException&) { break; }
			node.push_back(r);
			i++;
		}

		Node* n = createNode(pTree, node, 0);
		node.clear();
		pTree->writeNode(n);
		es2->insert(new ExternalSorter::Record(n->m_nodeMBR, n->m_identifier, 0, nullptr, 0));
		pTree->m_rootID = n->m_identifier;

		pTree->m_stats.m_u32TreeHeight = std::max(pTree->m_stats.m_u32TreeHeight, static_cast<uint32_t>(level));

		delete n;
	}
	else
	{
        std::shared_ptr<ExternalSorter> left_node_es = std::shared_ptr<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));
        std::shared_ptr<ExternalSorter> right_node_es = std::shared_ptr<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));
		ExternalSorter::Record* pR_left;
		ExternalSorter::Record* pR_right;

		uint32_t sort_dim_index = level % dimension;

		for (uint64_t i = 0; i < left_node_en; ++i)
		{
			try { pR_left = es->getNextRecord(); }
			catch (Tools::EndOfStreamException&) { break; }
			pR_left->m_s = sort_dim_index;
			left_node_es->insert(pR_left);
		}
		left_node_es->sort();
		
		partition(pTree, left_node_es, dimension, bleaf, bindex, level+1, es2, pageSize, numberOfPages);
		ExternalSorter::Record* r_left = es2->getNextRecord();


		for (uint64_t i = 0; i < right_node_en; ++i)
		{
			try { pR_right = es->getNextRecord(); }
			catch (Tools::EndOfStreamException&) { break; }
			pR_right->m_s = sort_dim_index;
			right_node_es->insert(pR_right);
		}
		right_node_es->sort();

		partition(pTree, right_node_es, dimension, bleaf, bindex, level+1, es2, pageSize, numberOfPages);
		ExternalSorter::Record* r_right = es2->getNextRecord();

		std::vector<ExternalSorter::Record*> parent;
		parent.push_back(r_left);
		parent.push_back(r_right);

		Node* n_parent = createNode(pTree, parent, pTree->m_stats.m_u32TreeHeight - level);

		parent.clear();
		pTree->writeNode(n_parent);

		es2->insert(new ExternalSorter::Record(n_parent->m_nodeMBR, n_parent->m_identifier, 0, nullptr, 0));
		pTree->m_rootID = n_parent->m_identifier;

		delete n_parent;
	}

}

Node* BulkLoader::createNode(SpatialIndex::KDTree::KDTree* pTree, std::vector<ExternalSorter::Record*>& e, uint32_t level)
{
	Node* n;

	if (level == 0) n = new Leaf(pTree, -1);
	else n = new Index(pTree, -1, level);

	for (size_t cChild = 0; cChild < e.size(); ++cChild)
	{
		n->insertEntry(e[cChild]->m_len, e[cChild]->m_pData, e[cChild]->m_r, e[cChild]->m_id);
		e[cChild]->m_pData = nullptr;
		delete e[cChild];
	}

	return n;
}
