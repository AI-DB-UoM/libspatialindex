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
#include <queue>
#include <iostream>
#include <fstream>

#ifndef _MSC_VER
#include <unistd.h>
#endif

#include <spatialindex/SpatialIndex.h>

#include "RTree.h"
#include "Leaf.h"
#include "Index.h"
#include "BulkLoader.h"

using namespace SpatialIndex;
using namespace SpatialIndex::RTree;

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
	if (m_bInsertionPhase == true)
		throw Tools::IllegalStateException("ExternalSorter::getNextRecord: Input has not been sorted yet.");

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
void BulkLoader::bulkLoadUsingSTR(
	SpatialIndex::RTree::RTree* pTree,
	IDataStream& stream,
	uint32_t bindex,
	uint32_t bleaf,
	uint32_t pageSize,
	uint32_t numberOfPages
) {
	if (! stream.hasNext())
		throw Tools::IllegalArgumentException(
			"RTree::BulkLoader::bulkLoadUsingSTR: Empty data stream given."
		);

	NodePtr n = pTree->readNode(pTree->m_rootID);
	pTree->deleteNode(n.get());

	#ifndef NDEBUG
	std::cerr << "RTree::BulkLoader: Sorting data." << std::endl;
	#endif

    std::shared_ptr<ExternalSorter> es = std::shared_ptr<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));

	while (stream.hasNext())
	{
		Data* d = reinterpret_cast<Data*>(stream.getNext());
		if (d == nullptr)
			throw Tools::IllegalArgumentException(
				"bulkLoadUsingSTR: RTree bulk load expects SpatialIndex::RTree::Data entries."
			);

		es->insert(new ExternalSorter::Record(d->m_region, d->m_id, d->m_dataLength, d->m_pData, 0));
		d->m_pData = nullptr;
		delete d;
	}
	es->sort();

	pTree->m_stats.m_u64Data = es->getTotalEntries();

	// create index levels.
	uint32_t level = 0;

	while (true)
	{
		#ifndef NDEBUG
		std::cerr << "RTree::BulkLoader: Building level " << level << std::endl;
		#endif

		pTree->m_stats.m_nodesInLevel.push_back(0);

        std::shared_ptr<ExternalSorter> es2 = std::shared_ptr<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));
		createLevel(pTree, es, 0, bleaf, bindex, level++, es2, pageSize, numberOfPages);
		es = es2;

		if (es->getTotalEntries() == 1) break;
		es->sort();
	}

	pTree->m_stats.m_u32TreeHeight = level;
	pTree->storeHeader();
}

//
// BulkLoader
//
void BulkLoader::bulkLoadUsingSFC(
	SpatialIndex::RTree::RTree* pTree,
	IDataStream& stream,
	uint32_t bindex,
	uint32_t bleaf,
	uint32_t pageSize,
	uint32_t numberOfPages
) {
	if (! stream.hasNext())
		throw Tools::IllegalArgumentException(
			"RTree::BulkLoader::bulkLoadUsingSFC: Empty data stream given."
		);

	NodePtr n = pTree->readNode(pTree->m_rootID);
	pTree->deleteNode(n.get());

	#ifndef NDEBUG
	std::cerr << "RTree::BulkLoader: Sorting data." << std::endl;
	#endif

    std::shared_ptr<ExternalSorter> es = std::shared_ptr<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));

	while (stream.hasNext())
	{
		Data* d = reinterpret_cast<Data*>(stream.getNext());
		if (d == nullptr)
			throw Tools::IllegalArgumentException(
				"bulkLoadUsingSFC: RTree bulk load expects SpatialIndex::RTree::Data entries."
			);

		es->insert(new ExternalSorter::Record(d->m_region, d->m_id, d->m_dataLength, d->m_pData, 0));
		d->m_pData = nullptr;
		delete d;
	}
	// es->sort();
	es->finishLoad();

	pTree->m_stats.m_u64Data = es->getTotalEntries();

	// create index levels.
	uint32_t level = 0;

	while (true)
	{
		#ifndef NDEBUG
		std::cerr << "RTree::BulkLoader: Building level " << level << std::endl;
		#endif

		pTree->m_stats.m_nodesInLevel.push_back(0);

        std::shared_ptr<ExternalSorter> es2 = std::shared_ptr<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));
		createLevelSFC(pTree, es, 0, bleaf, bindex, level++, es2, pageSize, numberOfPages);
		es = es2;

		if (es->getTotalEntries() == 1) break;
		// es->sort();
		es->finishLoad();
	}

	pTree->m_stats.m_u32TreeHeight = level;
	pTree->storeHeader();
}

void BulkLoader::createLevel(
	SpatialIndex::RTree::RTree* pTree,
	std::shared_ptr<ExternalSorter> es,
	uint32_t dimension,
	uint32_t bleaf,
	uint32_t bindex,
	uint32_t level,
	std::shared_ptr<ExternalSorter> es2,
	uint32_t pageSize,
	uint32_t numberOfPages
) {
	uint64_t b = (level == 0) ? bleaf : bindex;
	uint64_t P = static_cast<uint64_t>(std::ceil(static_cast<double>(es->getTotalEntries()) / static_cast<double>(b)));
	uint64_t S = static_cast<uint64_t>(std::ceil(std::sqrt(static_cast<double>(P))));

	if (S == 1 || dimension == pTree->m_dimension - 1 || S * b == es->getTotalEntries())
	{
		std::vector<ExternalSorter::Record*> node;
		ExternalSorter::Record* r;

		while (true)
		{
			try { r = es->getNextRecord(); } catch (Tools::EndOfStreamException&) { break; }
			node.push_back(r);

			if (node.size() == b)
			{
				Node* n = createNode(pTree, node, level);
				node.clear();
				pTree->writeNode(n);
				es2->insert(new ExternalSorter::Record(n->m_nodeMBR, n->m_identifier, 0, nullptr, 0));
				pTree->m_rootID = n->m_identifier;
					// special case when the root has exactly bindex entries.
				delete n;
			}
		}

		if (! node.empty())
		{
			Node* n = createNode(pTree, node, level);
			pTree->writeNode(n);
			es2->insert(new ExternalSorter::Record(n->m_nodeMBR, n->m_identifier, 0, nullptr, 0));
			pTree->m_rootID = n->m_identifier;
			delete n;
		}
	}
	else
	{
		bool bMore = true;

		while (bMore)
		{
			ExternalSorter::Record* pR;
            std::shared_ptr<ExternalSorter> es3 = std::shared_ptr<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));

			for (uint64_t i = 0; i < S * b; ++i)
			{
				try { pR = es->getNextRecord(); }
				catch (Tools::EndOfStreamException&) { bMore = false; break; }
				pR->m_s = dimension + 1;
				es3->insert(pR);
			}
			es3->sort();
			createLevel(pTree, es3, dimension + 1, bleaf, bindex, level, es2, pageSize, numberOfPages);
		}
	}
}


void BulkLoader::createLevelTGS(
	SpatialIndex::RTree::RTree* pTree,
	std::shared_ptr<std::vector<ExternalSorter::Record*>> orderedRecords,
	uint32_t dimension,
	uint32_t bleaf,
	uint32_t bindex,
	uint32_t level,
	std::shared_ptr<std::vector<ExternalSorter::Record*>> orderedRecords2,
	uint32_t pageSize,
	uint32_t numberOfPages
) {
	uint64_t b = (level == 0) ? bleaf : bindex;
	uint64_t P = static_cast<uint64_t>(std::ceil(static_cast<double>(orderedRecords->size() / static_cast<double>(b))));

	std::vector<ExternalSorter::Record*> node;
	ExternalSorter::Record* r;

	for (int i = 0; i < orderedRecords->size(); i++)
	{
		node.push_back((*orderedRecords)[i]);

		if (node.size() == b)
		{
			Node* n = createNode(pTree, node, level);
			node.clear();
			pTree->writeNode(n);
			orderedRecords2->push_back(new ExternalSorter::Record(n->m_nodeMBR, n->m_identifier, 0, nullptr, 0));
			pTree->m_rootID = n->m_identifier;
				// special case when the root has exactly bindex entries.
			delete n;
		}
	}

	if (! node.empty())
	{
		std::cout << node.size() << " records remains, put remaining records in one block\n";
		Node* n = createNode(pTree, node, level);
		pTree->writeNode(n);
		orderedRecords2->push_back(new ExternalSorter::Record(n->m_nodeMBR, n->m_identifier, 0, nullptr, 0));
		pTree->m_rootID = n->m_identifier;
		delete n;
	}
}


void BulkLoader::bulkLoadUsingTGS(
	SpatialIndex::RTree::RTree* pTree,
	IDataStream& stream,
	uint32_t bindex,
	uint32_t bleaf,
	uint32_t pageSize,
	uint32_t numberOfPages,
	std::string cutListFile
) {
	std::vector<std::pair<int, int>> cutList;
	int dim, pos;
	std::ifstream readFile(cutListFile);
	while (readFile >> dim >> pos) {
		cutList.push_back(std::make_pair(dim, pos));
	}


	if (! stream.hasNext())
		throw Tools::IllegalArgumentException(
			"RTree::BulkLoader::bulkLoadUsingTGS: Empty data stream given."
		);

	NodePtr n = pTree->readNode(pTree->m_rootID);
	pTree->deleteNode(n.get());

	#ifndef NDEBUG
	std::cerr << "RTree::BulkLoader: Sorting data." << std::endl;
	#endif

    // std::shared_ptr<ExternalSorter> es = std::shared_ptr<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));
	std::shared_ptr<std::vector<ExternalSorter::Record*>> es = std::shared_ptr<std::vector<ExternalSorter::Record*>>(new std::vector<ExternalSorter::Record*>);

	std::cout << "Start loading data...\n";

	while (stream.hasNext())
	{
		Data* d = reinterpret_cast<Data*>(stream.getNext());
		if (d == nullptr)
			throw Tools::IllegalArgumentException(
				"bulkLoadUsingTGS: RTree bulk load expects SpatialIndex::RTree::Data entries."
			);

		// es->insert(new ExternalSorter::Record(d->m_region, d->m_id, d->m_dataLength, d->m_pData, 0));
		es->push_back(new ExternalSorter::Record(d->m_region, d->m_id, d->m_dataLength, d->m_pData, 0));
		d->m_pData = nullptr;
		delete d;
	}
	// es->sort();
	std::sort(es->begin(), es->end(), ExternalSorter::Record::SortAscending());

	// pTree->m_stats.m_u64Data = es->getTotalEntries();
	pTree->m_stats.m_u64Data = es->size();

	// std::cout << es->getTotalEntries() << " Data loaded\n";
	std::cout << es->size() << " Data loaded\n";

	std::queue<PartitionState*> q;

	// uint32_t totalLevel = std::ceil(log(static_cast<double>(es->getTotalEntries()) / static_cast<double>(bleaf)) / log(static_cast<double>(bindex)));
	uint32_t totalLevel = std::ceil(log(static_cast<double>(es->size()) / static_cast<double>(bleaf)) / log(static_cast<double>(bindex)));

	std::shared_ptr<std::vector<ExternalSorter::Record*>> orderedRecords = std::shared_ptr<std::vector<ExternalSorter::Record*>>(new std::vector<ExternalSorter::Record*>);
	std::shared_ptr<std::vector<ExternalSorter::Record*>> orderedRecords2 = std::shared_ptr<std::vector<ExternalSorter::Record*>>(new std::vector<ExternalSorter::Record*>);
	
	int cutCount = 0;
	q.push(new PartitionState(es, totalLevel));
	while (!q.empty()) {
		PartitionState* currentState = q.front();
		q.pop();
		// std::vector<std::shared_ptr<ExternalSorter>> partitionList;
		std::vector<std::shared_ptr<std::vector<ExternalSorter::Record*>>> partitionList;
		partitionList.push_back(currentState->es);
		uint32_t currentLevel = currentState->level;
		uint64_t childSize = bleaf * pow(bindex, currentLevel - 1);
		// uint32_t nChilds = currentState->es->getTotalEntries() / childSize;
		uint32_t nChilds = ceil(static_cast<double>(currentState->es->size()) / childSize);
		std::cout << "currentLevel: " << currentLevel << ", childSize: " << childSize << ", nChilds: " << nChilds << "\n";
		for (int cCut = 0; cCut < nChilds - 1; cCut++){
			int currentIdx = 0;
			for (int i = 0; i < partitionList.size(); i++){
				// if (partitionList[i]->getTotalEntries() > childSize) {
				if (partitionList[i]->size() > childSize) {
					currentIdx = i;
					break;
				}
			}
			// std::cout << "partition index: " << currentIdx << ", size: " << partitionList[currentIdx]->getTotalEntries() << "\n";
			std::cout << "partition index: " << currentIdx << ", size: " << partitionList[currentIdx]->size() << "\n";
			// std::shared_ptr<ExternalSorter> currentEs = partitionList[currentIdx];
			std::shared_ptr<std::vector<ExternalSorter::Record*>> currentEs = partitionList[currentIdx];
			// uint64_t currentTotalEntries = currentEs->getTotalEntries();
			uint64_t currentTotalEntries = currentEs->size();
			std::cout << "currentTotalEntries: " << currentTotalEntries << "\n";
			std::cout << "cutCount: " << cutCount << "\n";

			ExternalSorter::Record* r;
			std::pair<int, int> currentCut = cutList[cutCount];
			
			cutCount++;
			int dim = currentCut.first;
			int pos = currentCut.second;
			std::cout << "Start sorting by cut dimension: " << dim << " " << pos << "\n";
			// std::shared_ptr<ExternalSorter> intermediateEs = std::shared_ptr<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));
			std::shared_ptr<std::vector<ExternalSorter::Record*>> intermediateEs = std::shared_ptr<std::vector<ExternalSorter::Record*>>(new std::vector<ExternalSorter::Record*>);
			// for (int i = 0; i < currentTotalEntries; i++){
			// 	try { r = currentEs->getNextRecord(); } catch (Tools::EndOfStreamException&) { break; }
			// 	r->m_s = dim;
			// 	intermediateEs->insert(r);
			// }
			for (int i = 0; i < currentTotalEntries; i++) {
				r = (*currentEs)[i];
				r->m_s = dim;
				intermediateEs->push_back(r);
			}
			std::sort(intermediateEs->begin(), intermediateEs->end(), ExternalSorter::Record::SortAscending());
			currentEs = intermediateEs;
			std::cout << "Sorted by cut dimension: " << dim << " " << pos << "\n";
			// std::shared_ptr<ExternalSorter> es1 = std::shared_ptr<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));
			// std::shared_ptr<ExternalSorter> es2 = std::shared_ptr<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));
			std::shared_ptr<std::vector<ExternalSorter::Record*>> es1 = std::shared_ptr<std::vector<ExternalSorter::Record*>>(new std::vector<ExternalSorter::Record*>);
			std::shared_ptr<std::vector<ExternalSorter::Record*>> es2 = std::shared_ptr<std::vector<ExternalSorter::Record*>>(new std::vector<ExternalSorter::Record*>);
			// for (int i = 0; i < pos * childSize; i++) {
			// 	try { r = currentEs->getNextRecord(); } catch (Tools::EndOfStreamException&) { break; }
			// 	es1->insert(r);
			// }
			// es1->sort();
			// for (int i = 0; i < currentEs->getTotalEntries() - pos * childSize; i++) {
			// 	try { r = currentEs->getNextRecord(); } catch (Tools::EndOfStreamException&) { break; }
			// 	es2->insert(r);
			// }
			// es2->sort();
			// std::cout << "currentEs->size(): " << currentEs->size() << " \n";

			for (int i = 0; i < pos * childSize; i++) {
				r = (*currentEs)[i];
				es1->push_back(r);
			}

			for (int i = pos * childSize; i < currentEs->size(); i++) {
				r = (*currentEs)[i];
				es2->push_back(r);
			}

			std::vector<std::shared_ptr<std::vector<ExternalSorter::Record*>>>::iterator it;
			it = partitionList.begin();
			partitionList.erase(it+currentIdx);
			partitionList.insert(it+currentIdx, es1);
			partitionList.insert(it+currentIdx+1, es2);
			// std::cout << "partitioned into size: " << es1->getTotalEntries() << " " << es2->getTotalEntries() << " \n";
			std::cout << "partitioned into size: " << es1->size() << " " << es2->size() << " \n";
		}
		if (currentState->level == 1) {
			ExternalSorter::Record* r;
			for (int i = 0; i < partitionList.size(); i++){
				// while (true) {
				// 	try { r = partitionList[i]->getNextRecord(); } catch (Tools::EndOfStreamException&) { break; }
				// 	orderedRecords->push_back(r);
				// }
				for (int j = 0; j < partitionList[i]->size(); j++) {
					r = (*partitionList[i])[j];
					orderedRecords->push_back(r);
				}
			}
		}
		else {
			for (int i = 0; i < partitionList.size(); i++){
				q.push(new PartitionState(partitionList[i], currentLevel-1));
			}
		}
	}

	std::cout << "Ordered records in TGS order\n";

	

	// create index levels.
	uint32_t level = 0;

	while (true)
	{
		#ifndef NDEBUG
		std::cerr << "RTree::BulkLoader: Building level " << level << std::endl;
		#endif
		std::cout << "RTree::BulkLoader: Building level " << level << std::endl;
		std::cout << "Number of records: " << orderedRecords->size() << std::endl;
		if (level == 2 || level == 3) {
			for (int i = 0; i < orderedRecords->size(); i++) {
				std::cout << (*orderedRecords)[i]->m_r << std::endl;
			}
		}
		std::shared_ptr<std::vector<ExternalSorter::Record*>> orderedRecords2 = std::shared_ptr<std::vector<ExternalSorter::Record*>>(new std::vector<ExternalSorter::Record*>);

		pTree->m_stats.m_nodesInLevel.push_back(0);

        std::shared_ptr<ExternalSorter> es2 = std::shared_ptr<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));
		createLevelTGS(pTree, orderedRecords, 0, bleaf, bindex, level++, orderedRecords2, pageSize, numberOfPages);
		orderedRecords = orderedRecords2;

		if (orderedRecords->size() == 1) break;
		// es->sort();
	}
	std::cout << "Finished building all levels\n";

	pTree->m_stats.m_u32TreeHeight = level;
	pTree->storeHeader();
}

void BulkLoader::createLevelSFC(
	SpatialIndex::RTree::RTree* pTree,
	std::shared_ptr<ExternalSorter> es,
	uint32_t dimension,
	uint32_t bleaf,
	uint32_t bindex,
	uint32_t level,
	std::shared_ptr<ExternalSorter> es2,
	uint32_t pageSize,
	uint32_t numberOfPages
) {
	uint64_t b = (level == 0) ? bleaf : bindex;
	uint64_t P = static_cast<uint64_t>(std::ceil(static_cast<double>(es->getTotalEntries()) / static_cast<double>(b)));
	uint64_t S = static_cast<uint64_t>(std::ceil(std::sqrt(static_cast<double>(P))));

	if (S == 1 || dimension == pTree->m_dimension - 1 || S * b == es->getTotalEntries())
	{
		std::vector<ExternalSorter::Record*> node;
		ExternalSorter::Record* r;

		while (true)
		{
			try { r = es->getNextRecord(); } catch (Tools::EndOfStreamException&) { break; }
			node.push_back(r);

			if (node.size() == b)
			{
				Node* n = createNode(pTree, node, level);
				node.clear();
				pTree->writeNode(n);
				es2->insert(new ExternalSorter::Record(n->m_nodeMBR, n->m_identifier, 0, nullptr, 0));
				pTree->m_rootID = n->m_identifier;
					// special case when the root has exactly bindex entries.
				delete n;
			}
		}

		if (! node.empty())
		{
			Node* n = createNode(pTree, node, level);
			pTree->writeNode(n);
			es2->insert(new ExternalSorter::Record(n->m_nodeMBR, n->m_identifier, 0, nullptr, 0));
			pTree->m_rootID = n->m_identifier;
			delete n;
		}
	}
	else
	{
		bool bMore = true;

		while (bMore)
		{
			ExternalSorter::Record* pR;
            std::shared_ptr<ExternalSorter> es3 = std::shared_ptr<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));

			for (uint64_t i = 0; i < S * b; ++i)
			{
				try { pR = es->getNextRecord(); }
				catch (Tools::EndOfStreamException&) { bMore = false; break; }
				pR->m_s = dimension + 1;
				es3->insert(pR);
			}
			es3->finishLoad();
			createLevelSFC(pTree, es3, dimension + 1, bleaf, bindex, level, es2, pageSize, numberOfPages);
		}
	}
}

Node* BulkLoader::createNode(SpatialIndex::RTree::RTree* pTree, std::vector<ExternalSorter::Record*>& e, uint32_t level)
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

//
// PartitionState
//
PartitionState::PartitionState(std::shared_ptr<std::vector<ExternalSorter::Record*>> es, uint32_t level)
: es(es), level(level)
{
}

PartitionState::~PartitionState()
{
	es = nullptr;
}