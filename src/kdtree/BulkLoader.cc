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
#include <random>
#include <ctime>

#include <algorithm>

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

    // std::shared_ptr<ExternalSorter> es = std::shared_ptr<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));
	std::vector<ExternalSorter::Record*> es;

	while (stream.hasNext())
	{
		Data* d = reinterpret_cast<Data*>(stream.getNext());
		if (d == nullptr)
			throw Tools::IllegalArgumentException(
				"topDownPartitioning: KDTree bulk load expects SpatialIndex::KDTree::Data entries."
			);

		// es->insert(new ExternalSorter::Record(d->m_region, d->m_id, d->m_dataLength, d->m_pData, 0));
		es.push_back(new ExternalSorter::Record(d->m_region, d->m_id, d->m_dataLength, d->m_pData, 0));
		d->m_pData = nullptr;
		delete d;
	}
	// es->sort();
	std::sort(es.begin(), es.end(), ExternalSorter::Record::SortAscending());

	// pTree->m_stats.m_u64Data = es->getTotalEntries();
	pTree->m_stats.m_u64Data = es.size();

	pTree->m_stats.m_nodesInLevel.push_back(0);

	// std::shared_ptr<ExternalSorter> es2 = std::shared_ptr<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));
	std::vector<ExternalSorter::Record *> es2;
	partition(pTree, es, pTree->m_dimension, bleaf, bindex, 1, es2, pageSize, numberOfPages);

	pTree->storeHeader();
}


void BulkLoader::partition(
    SpatialIndex::KDTree::KDTree* pTree,
    std::vector<ExternalSorter::Record *> es,
    uint32_t dimension,
    uint32_t bleaf,
    uint32_t bindex,
    uint32_t level,
	std::vector<ExternalSorter::Record *>& es2,
    uint32_t pageSize,
    uint32_t numberOfPages
) {
    // uint64_t total_entries = es->getTotalEntries();
	uint64_t total_entries = es.size();
    uint64_t left_node_en = static_cast<uint64_t>(total_entries / 2);
    uint64_t right_node_en = static_cast<uint64_t>(total_entries - left_node_en);

    if (pTree->m_stats.m_nodesInLevel.size() < level) {
        pTree->m_stats.m_nodesInLevel.push_back(0);
    }

    if (total_entries <= bleaf) {
        // std::vector<ExternalSorter::Record*> node;
        // ExternalSorter::Record* r;
        // uint32_t i = 0;

        // while (i < total_entries) {
        //     try { r = es->getNextRecord(); } catch (Tools::EndOfStreamException&) { break; }
        //     node.push_back(r);
        //     i++;
        // }

        Node* n = createNode(pTree, es, 0);
        es.clear();
        pTree->writeNode(n);
        es2.push_back(new ExternalSorter::Record(n->m_nodeMBR, n->m_identifier, 0, nullptr, 0));
        pTree->m_rootID = n->m_identifier;

        pTree->m_stats.m_u32TreeHeight = std::max(pTree->m_stats.m_u32TreeHeight, static_cast<uint32_t>(level));

        delete n;
    } else {
        // auto left_node_es = std::make_shared<ExternalSorter>(pageSize, numberOfPages);
        // auto right_node_es = std::make_shared<ExternalSorter>(pageSize, numberOfPages);
		std::vector<ExternalSorter::Record *> left_node_es;
		std::vector<ExternalSorter::Record *> right_node_es;
        ExternalSorter::Record* pR_left;
        ExternalSorter::Record* pR_right;

        uint32_t sort_dim_index = level % dimension;

        for (uint64_t i = 0; i < left_node_en; ++i) {
            try { pR_left = es[i]; }
            catch (Tools::EndOfStreamException&) { break; }
            pR_left->m_s = sort_dim_index;
			left_node_es.push_back(pR_left);
            // left_node_es->insert(pR_left);
        }
        // left_node_es->sort();
		std::sort(left_node_es.begin(), left_node_es.end(), ExternalSorter::Record::SortAscending());

		for (uint64_t i = left_node_en; i < total_entries; ++i) {
            try { pR_right = es[i]; }
            catch (Tools::EndOfStreamException&) { break; }
            pR_right->m_s = sort_dim_index;
			right_node_es.push_back(pR_right);
            // right_node_es->insert(pR_right);
        }
        // right_node_es->sort();
		std::sort(right_node_es.begin(), right_node_es.end(), ExternalSorter::Record::SortAscending());

        partition(pTree, left_node_es, dimension, bleaf, bindex, level + 1, es2, pageSize, numberOfPages);
        ExternalSorter::Record* r_left = es2.back();
		es2.pop_back();

        partition(pTree, right_node_es, dimension, bleaf, bindex, level + 1, es2, pageSize, numberOfPages);
        ExternalSorter::Record* r_right = es2.back();
		es2.pop_back();

        std::vector<ExternalSorter::Record*> parent;
        parent.push_back(r_left);
        parent.push_back(r_right);

        Node* n_parent = createNode(pTree, parent, pTree->m_stats.m_u32TreeHeight - level);

        parent.clear();
        pTree->writeNode(n_parent);

        es2.push_back(new ExternalSorter::Record(n_parent->m_nodeMBR, n_parent->m_identifier, 0, nullptr, 0));
        pTree->m_rootID = n_parent->m_identifier;

        delete n_parent;
    }
}


std::vector<std::pair<uint32_t, double>> BulkLoader::generageCandidateCutPos(std::vector<Region>& regions)
{
	std::vector<std::pair<uint32_t, double>> candidateCutPos;
	for (size_t i = 0; i < regions.size(); ++i)
	{
		for (size_t j = 0; j < regions[i].m_dimension; ++j)
		{
			candidateCutPos.push_back(std::make_pair(j, regions[i].m_pLow[j]));
			candidateCutPos.push_back(std::make_pair(j, regions[i].m_pHigh[j]));
		}
	}
	return candidateCutPos;
}

std::vector<std::pair<uint32_t, double>> BulkLoader::generageModelCandidateCutPos(int dimension, std::vector<Region>& regions, int sample_size)
{
	assert(!regions.empty());
    std::vector<std::vector<std::pair<uint32_t, double>>> actions(dimension);
	for (int i = 0; i < dimension; ++i) {
		std::vector<std::pair<uint32_t, double>> dim_actions;
		for (const auto& region : regions) {
			dim_actions.emplace_back(i, region.m_pLow[i]);
			dim_actions.emplace_back(i, region.m_pHigh[i]);
		}
		// Sort actions based on the cut positions
		std::sort(dim_actions.begin(), dim_actions.end(), [](const auto& a, const auto& b) {
			return a.second < b.second;
		});
		actions[i] = std::move(dim_actions);
	}

	int sample_size_per_dim = sample_size / dimension;
	sample_size_per_dim += 1;  // Adjust for the first entry removal

	std::vector<std::vector<std::pair<uint32_t, double>>> sampled_actions(dimension);
	for (int i = 0; i < dimension; ++i) {
		int total_actions = actions[i].size();
		int interval = total_actions / sample_size_per_dim;
		assert(interval > 0);
		for (int j = 0; j < sample_size_per_dim && j * interval < total_actions; ++j) {
			sampled_actions[i].push_back(actions[i][j * interval]);
		}
		
		if (!sampled_actions[i].empty()) {
			sampled_actions[i].erase(sampled_actions[i].begin());  // Remove the first sampled element
		}
	}
	actions = std::move(sampled_actions);

	// Flatten the actions vector
	std::vector<std::pair<uint32_t, double>> candidateCutPos;
	for (const auto& action_list : actions) {
		candidateCutPos.insert(candidateCutPos.end(), action_list.begin(), action_list.end());
	}
	return candidateCutPos;
}

uint64_t BulkLoader::calculateSkip(std::vector<std::pair<Region, uint64_t>>& candidateSplits, std::vector<Region>& regions, uint32_t dimension, double value)
{
	uint64_t skip;
	for (const auto& candidateSplit : candidateSplits) {
		const Region& splitRegion = candidateSplit.first;
		uint64_t splitCount = candidateSplit.second; 

		for (const Region& region : regions) {
			if (!splitRegion.intersectsRegion(region)) {
				skip += splitCount;
			}
		}
	}
	return skip;
}

void BulkLoader::topDownGreedyPartitioning(
	SpatialIndex::KDTree::KDTree* pTree,
	IDataStream& stream,
	IDataStream& queryStream,
	uint32_t bindex,
	uint32_t bleaf,
	uint32_t pageSize,
	uint32_t numberOfPages
) {
	if (! stream.hasNext())
		throw Tools::IllegalArgumentException(
			"KDTree::BulkLoader::topDownGreedyPartitioning: Empty data stream given."
		);

	NodePtr n = pTree->readNode(pTree->m_rootID);
	pTree->deleteNode(n.get());

	#ifndef NDEBUG
	std::cerr << "KDTree::BulkLoader: Sorting data." << std::endl;
	#endif

	Region parentMBR = pTree->m_infiniteRegion;
	std::vector<ExternalSorter::Record*> tupleSet;
	while (stream.hasNext())
	{
		Data* d = reinterpret_cast<Data*>(stream.getNext());
		if (d == nullptr)
			throw Tools::IllegalArgumentException(
				"topDownGreedyPartitioning: KDTree bulk load expects SpatialIndex::KDTree::Data entries."
			);
		// es->insert(new ExternalSorter::Record(d->m_region, d->m_id, d->m_dataLength, d->m_pData, 0));
		tupleSet.push_back(new ExternalSorter::Record(d->m_region, d->m_id, d->m_dataLength, d->m_pData, 0));
		parentMBR.combineRegion(d->m_region);
		d->m_pData = nullptr;
		delete d;
	}
	// std::cerr << " before sort" << tupleSet.size() << std::endl;
	std::sort(tupleSet.begin(), tupleSet.end(), ExternalSorter::Record::SortAscending());
	// std::cerr << " after sort" << std::endl;

	std::vector<Region> allqueryRegions;
	while (queryStream.hasNext())
	{
		Data* d = reinterpret_cast<Data*>(queryStream.getNext());
		allqueryRegions.push_back(d->m_region);
		delete d;
	}

	std::vector<Region> queryRegions;

	std::sample(allqueryRegions.begin(), allqueryRegions.end(), std::back_inserter(queryRegions), 100, std::mt19937{std::random_device{}()});

	// std::vector<Region> queryRegions;
	// while (queryStream.hasNext())
	// {
	// 	Data* d = reinterpret_cast<Data*>(queryStream.getNext());
	// 	queryRegions.push_back(d->m_region);
	// 	delete d;
	// }

	pTree->m_stats.m_u64Data = tupleSet.size();

	pTree->m_stats.m_nodesInLevel.push_back(0);

	std::vector<ExternalSorter::Record *> tupleSet2;
	std::vector<std::pair<uint32_t, double>> candidateCutPos = generageCandidateCutPos(queryRegions);

	greedyPartition(pTree, tupleSet, bleaf, bindex, 1, parentMBR, tupleSet2, queryRegions, candidateCutPos);

	pTree->storeHeader();
}

bool BulkLoader::greedyPartition(
	SpatialIndex::KDTree::KDTree* pTree,
	std::vector<ExternalSorter::Record *> tupleSet,
	uint32_t bleaf,
	uint32_t bindex,
	uint32_t level,
	Region parentMBR,
	std::vector<ExternalSorter::Record *>& tupleSet2,
	std::vector<Region>& queryRegions,
	std::vector<std::pair<uint32_t, double>>& candidateCutPos)
{

	bool can_split = false;
	uint64_t total_entries = tupleSet.size();

	// std::cerr << "KDTree::greedyPartition: total_entries: " << total_entries << std::endl;
	if (pTree->m_stats.m_nodesInLevel.size() < level)
	{
		pTree->m_stats.m_nodesInLevel.push_back(0);
	}

	if (total_entries <= bleaf)
	{
		// std::cerr << "KDTree::greedyPartition: ----leaf----" << " level: " << level << std::endl;
		
		Node* n = createNode(pTree, tupleSet, 0);
		tupleSet.clear();
		pTree->writeNode(n);
		tupleSet2.push_back(new ExternalSorter::Record(n->m_nodeMBR, n->m_identifier, 0, nullptr, 0));
		pTree->m_rootID = n->m_identifier;
		pTree->m_stats.m_u32TreeHeight = std::max(pTree->m_stats.m_u32TreeHeight, static_cast<uint32_t>(level));
		delete n;
		return can_split;
	}
	// std::cerr << "KDTree::greedyPartition: ----node----" << " level: " << level << std::endl;

	uint64_t maxSkip = std::numeric_limits<uint64_t>::min();
	uint32_t maxSplitDimension = 0;
	double maxSplitValue;

	std::vector<ExternalSorter::Record *> bestLeftData;
	std::vector<ExternalSorter::Record *> bestRightData;

	Region bestLeftMBR = pTree->m_infiniteRegion;
	Region bestRightMBR = pTree->m_infiniteRegion;

	for (size_t j = 0; j < candidateCutPos.size(); j++)
	{
		uint32_t splitDimension = candidateCutPos[j].first;
		double splitValue = candidateCutPos[j].second;
		if (parentMBR.getLow(splitDimension) >= splitValue || parentMBR.getHigh(splitDimension) <= splitValue)
			continue;
		// create two sub Regions.
		double *low = new double[pTree->m_dimension];
		double *high = new double[pTree->m_dimension];

		std::copy(parentMBR.m_pLow, parentMBR.m_pLow + pTree->m_dimension, low);
		std::copy(parentMBR.m_pHigh, parentMBR.m_pHigh + pTree->m_dimension, high);
		high[splitDimension] = splitValue;
		Region leftMBR(low, high, pTree->m_dimension);

		std::copy(parentMBR.m_pLow, parentMBR.m_pLow + pTree->m_dimension, low);
		std::copy(parentMBR.m_pHigh, parentMBR.m_pHigh + pTree->m_dimension, high);
		low[splitDimension] = splitValue;
		Region rightMBR(low, high, pTree->m_dimension);

		std::vector<ExternalSorter::Record *> leftData;
		std::vector<ExternalSorter::Record *> rightData;

		for (uint64_t l = 0; l < total_entries; ++l)
		{
			ExternalSorter::Record *r = tupleSet[l];
            r->m_s = splitDimension;
			if (leftMBR.intersectsRegion(r->m_r))
			{
				leftData.push_back(r);
				continue;
			}
			if (rightMBR.intersectsRegion(r->m_r))
			{
				rightData.push_back(r);
				continue;
			}
			throw Tools::IllegalStateException("BulkLoader::topDownGreedyPartitioning: data missing!");
		}

		delete[] low;
		delete[] high;

		//  chech if one of the candidlate split leads to super small split.
		if (leftData.size() <= bleaf || rightData.size() <= bleaf)
			continue;

		std::vector<std::pair<Region, uint64_t>> candidateSplits;
		candidateSplits.push_back(std::make_pair(leftMBR, leftData.size()));
		candidateSplits.push_back(std::make_pair(rightMBR, rightData.size()));

		// calculate the number of split
		uint64_t skip = calculateSkip(candidateSplits, queryRegions, splitDimension, splitValue);

		if (skip > maxSkip)
		{
			maxSkip = skip;
			maxSplitDimension = splitDimension;
			maxSplitValue = splitValue;
			bestLeftData = std::move(leftData);
			bestRightData = std::move(rightData);
			bestLeftMBR = leftMBR; 
        	bestRightMBR = rightMBR; 
		}

		if (maxSkip > 0)
		{
			can_split = true;
		}
	}

	// std::cerr << "KDTree::greedyPartition: ----can_split----" << can_split << std::endl;

	if (!can_split) 
	{
		uint64_t left_node_en = static_cast<uint64_t>(total_entries / 2);
		std::vector<ExternalSorter::Record *> leftData;
		std::vector<ExternalSorter::Record *> rightData;
		Region leftMBR = pTree->m_infiniteRegion;
		Region rightMBR = pTree->m_infiniteRegion;

        // uint32_t sort_dim_index = level % dimension;

		for (uint64_t l = 0; l < left_node_en; ++l)
		{
			ExternalSorter::Record *r = tupleSet[l];
            r->m_s = maxSplitDimension;
			leftData.push_back(r);
			leftMBR.combineRegion(r->m_r);
		}
		for (uint64_t l = left_node_en; l < total_entries; ++l)
		{
			ExternalSorter::Record *r = tupleSet[l];
            r->m_s = maxSplitDimension;
			rightData.push_back(r);
			rightMBR.combineRegion(r->m_r);
		}
		bestLeftData = std::move(leftData);
		bestRightData = std::move(rightData);
		bestLeftMBR = leftMBR;
		bestRightMBR = rightMBR;
	}
	// std::cerr << " before sort left" << bestLeftData.size() << std::endl;
	std::sort(bestLeftData.begin(), bestLeftData.end(), ExternalSorter::Record::SortAscending());
	// std::cerr << " after sort left" << std::endl;

	greedyPartition(pTree, bestLeftData, bleaf, bindex, level+1, bestLeftMBR, tupleSet2, queryRegions, candidateCutPos);
	// ExternalSorter::Record* r_left = tupleSet2
	ExternalSorter::Record* r_left = tupleSet2.back();
	tupleSet2.pop_back();

	// std::cerr << " before sort right" << bestRightData.size() << std::endl;
	std::sort(bestRightData.begin(), bestRightData.end(), ExternalSorter::Record::SortAscending());
	// std::cerr << " after sort right" << std::endl;

	greedyPartition(pTree, bestRightData, bleaf, bindex, level+1, bestRightMBR, tupleSet2, queryRegions, candidateCutPos);

	ExternalSorter::Record* r_right = tupleSet2.back();

	tupleSet2.pop_back();


	std::vector<ExternalSorter::Record*> parent;
	parent.push_back(r_left);
	parent.push_back(r_right);

	Node* n_parent = createNode(pTree, parent, pTree->m_stats.m_u32TreeHeight - level);
	parent.clear();
	pTree->writeNode(n_parent);

	tupleSet2.push_back(new ExternalSorter::Record(n_parent->m_nodeMBR, n_parent->m_identifier, 0, nullptr, 0));
	pTree->m_rootID = n_parent->m_identifier;

	delete n_parent;

	return can_split;
}

void BulkLoader::topDownModelPartitioning(
	SpatialIndex::KDTree::KDTree* pTree,
	IDataStream& stream,
	IDataStream& queryStream,
	uint32_t bindex,
	uint32_t bleaf,
	uint32_t pageSize,
	uint32_t numberOfPages,
	const std::string& modelPath,
	int action_space_sizes
) 
{
	if (! stream.hasNext())
		throw Tools::IllegalArgumentException(
			"KDTree::BulkLoader::topDownGreedyPartitioning: Empty data stream given."
		);

	if (modelPath.empty())
	{
		throw Tools::IllegalArgumentException("Model path is empty, exiting.");
	}
	try
	{
		torch::Device device(torch::kCPU);
		pTree->m_splitModel = torch::jit::load(modelPath + "/qdtree.pth");
		pTree->m_splitModel.to(device);
	}
	catch(const std::exception& e)
	{
		throw Tools::IllegalArgumentException("Split Model path is empty, exiting.");
	}

	NodePtr n = pTree->readNode(pTree->m_rootID);
	pTree->deleteNode(n.get());

	#ifndef NDEBUG
	std::cerr << "KDTree::BulkLoader: Sorting data." << std::endl;
	#endif

	Region parentMBR = pTree->m_infiniteRegion;
	std::vector<ExternalSorter::Record*> tupleSet;
	while (stream.hasNext())
	{
		Data* d = reinterpret_cast<Data*>(stream.getNext());
		if (d == nullptr)
			throw Tools::IllegalArgumentException(
				"topDownGreedyPartitioning: KDTree bulk load expects SpatialIndex::KDTree::Data entries."
			);
		// es->insert(new ExternalSorter::Record(d->m_region, d->m_id, d->m_dataLength, d->m_pData, 0));
		tupleSet.push_back(new ExternalSorter::Record(d->m_region, d->m_id, d->m_dataLength, d->m_pData, 0));
		parentMBR.combineRegion(d->m_region);
		d->m_pData = nullptr;
		delete d;
	}
	std::sort(tupleSet.begin(), tupleSet.end(), ExternalSorter::Record::SortAscending());

	std::vector<Region> queryRegions;
	while (queryStream.hasNext())
	{
		Data* d = reinterpret_cast<Data*>(queryStream.getNext());
		queryRegions.push_back(d->m_region);
		delete d;
	}

	pTree->m_stats.m_u64Data = tupleSet.size();

	pTree->m_stats.m_nodesInLevel.push_back(0);

	std::vector<ExternalSorter::Record *> tupleSet2;
	std::vector<std::pair<uint32_t, double>> candidateCutPos = generageModelCandidateCutPos(pTree->m_dimension, queryRegions, action_space_sizes);

	modelPartition(pTree, tupleSet, bleaf, bindex, 1, parentMBR, tupleSet2, queryRegions, candidateCutPos);

	pTree->storeHeader();
}


std::mt19937 generator(static_cast<unsigned int>(std::time(nullptr)));
std::uniform_real_distribution<double> distribution(0.0, 1.0);

bool BulkLoader::modelPartition(
	SpatialIndex::KDTree::KDTree* pTree,
	std::vector<ExternalSorter::Record *> tupleSet,
	uint32_t bleaf,
	uint32_t bindex,
	uint32_t level,
	Region parentMBR,
	std::vector<ExternalSorter::Record *>& tupleSet2,
	std::vector<Region>& queryRegions,
	std::vector<std::pair<uint32_t, double>>& candidateCutPos)
{

	bool can_split = false;
	uint64_t total_entries = tupleSet.size();

	if (pTree->m_stats.m_nodesInLevel.size() < level)
	{
		pTree->m_stats.m_nodesInLevel.push_back(0);
	}

	if (total_entries <= bleaf)
	{
		Node* n = createNode(pTree, tupleSet, 0);
		tupleSet.clear();
		pTree->writeNode(n);
		tupleSet2.push_back(new ExternalSorter::Record(n->m_nodeMBR, n->m_identifier, 0, nullptr, 0));
		pTree->m_rootID = n->m_identifier;
		pTree->m_stats.m_u32TreeHeight = std::max(pTree->m_stats.m_u32TreeHeight, static_cast<uint32_t>(level));
		delete n;
		return can_split;
	}
	// std::cerr << "KDTree::modelPartition: ----node----" << " level: " << level << std::endl;

	std::vector<ExternalSorter::Record *> bestLeftData;
	std::vector<ExternalSorter::Record *> bestRightData;

	Region bestLeftMBR = pTree->m_infiniteRegion;
	Region bestRightMBR = pTree->m_infiniteRegion;

	// generate states
	std::vector<int> state;

	for (uint32_t i = 0; i < pTree->m_dimension; ++i) {
		std::vector<int> low_bits = float_to_bit_array(parentMBR.getLow(i));
		std::vector<int> high_bits = float_to_bit_array(parentMBR.getHigh(i));
		state.insert(state.end(), low_bits.begin(), low_bits.end());
		state.insert(state.end(), high_bits.begin(), high_bits.end());
	}



	double randomValue = distribution(generator);
	uint32_t action = 0;

	if (randomValue >= 0.5) 
	{
		std::uniform_int_distribution<uint32_t> actionDistribution(0, candidateCutPos.size() - 1);
		action = actionDistribution(generator);
	}
	else
	{
		action = pTree->splitModelForward(state);
	}

	uint32_t splitDimension = candidateCutPos[action].first;
	double splitValue = candidateCutPos[action].second;
	if (parentMBR.getLow(splitDimension) >= splitValue || parentMBR.getHigh(splitDimension) <= splitValue)
		can_split = false;

	// create two sub Regions.
	double *low = new double[pTree->m_dimension];
	double *high = new double[pTree->m_dimension];

	std::copy(parentMBR.m_pLow, parentMBR.m_pLow + pTree->m_dimension, low);
	std::copy(parentMBR.m_pHigh, parentMBR.m_pHigh + pTree->m_dimension, high);
	high[splitDimension] = splitValue;
	Region leftMBR(low, high, pTree->m_dimension);

	std::copy(parentMBR.m_pLow, parentMBR.m_pLow + pTree->m_dimension, low);
	std::copy(parentMBR.m_pHigh, parentMBR.m_pHigh + pTree->m_dimension, high);
	low[splitDimension] = splitValue;
	Region rightMBR(low, high, pTree->m_dimension);

	std::vector<ExternalSorter::Record *> leftData;
	std::vector<ExternalSorter::Record *> rightData;

	for (uint64_t l = 0; l < total_entries; ++l)
	{
		ExternalSorter::Record *r = tupleSet[l];
		r->m_s = splitDimension;
		if (leftMBR.intersectsRegion(r->m_r))
		{
			leftData.push_back(r);
			continue;
		}
		if (rightMBR.intersectsRegion(r->m_r))
		{
			rightData.push_back(r);
			continue;
		}
		throw Tools::IllegalStateException("BulkLoader::topDownModelPartitioning: data missing!");
	}

	delete[] low;
	delete[] high;

	//  chech if one of the candidlate split leads to super small split.
	if (can_split)
	{
		if (leftData.size() <= bleaf || rightData.size() <= bleaf)
		{
			can_split = false;
		}
		else
		{
			bestLeftData = std::move(leftData);
			bestRightData = std::move(rightData);
			bestLeftMBR = leftMBR; 
			bestRightMBR = rightMBR; 
		}
	}

	if (!can_split) 
	{
		uint64_t left_node_en = static_cast<uint64_t>(total_entries / 2);
		std::vector<ExternalSorter::Record *> leftData_;
		std::vector<ExternalSorter::Record *> rightData_;
		Region leftMBR_ = pTree->m_infiniteRegion;
		Region rightMBR_ = pTree->m_infiniteRegion;

		for (uint64_t l = 0; l < left_node_en; ++l)
		{
			ExternalSorter::Record *r = tupleSet[l];
			r->m_s = splitDimension;
			leftData_.push_back(r);
			leftMBR_.combineRegion(r->m_r);
		}
		for (uint64_t l = left_node_en; l < total_entries; ++l)
		{
			ExternalSorter::Record *r = tupleSet[l];
			r->m_s = splitDimension;
			rightData_.push_back(r);
			rightMBR_.combineRegion(r->m_r);
		}
		bestLeftData = std::move(leftData_);
		bestRightData = std::move(rightData_);
		bestLeftMBR = leftMBR_;
		bestRightMBR = rightMBR_;
	}

	// std::cerr << "bestLeftData: " << " bestLeftData[0]->m_s: " << bestLeftData[0]->m_s << std::endl;
	std::sort(bestLeftData.begin(), bestLeftData.end(), ExternalSorter::Record::SortAscending());
	// for (int i = 0; i < bestLeftData.size(); ++i)
	// {
	// 	std::cerr << bestLeftData[i]->m_r.m_pLow[0] << "," << bestLeftData[i]->m_r.m_pHigh[0] << "  " << bestLeftData[i]->m_r.m_pLow[1] << "," << bestLeftData[i]->m_r.m_pHigh[1] << std::endl;
	// }

	modelPartition(pTree, bestLeftData, bleaf, bindex, level+1, bestLeftMBR, tupleSet2, queryRegions, candidateCutPos);
	ExternalSorter::Record* r_left = tupleSet2.back();
	tupleSet2.pop_back();

	std::sort(bestRightData.begin(), bestRightData.end(), ExternalSorter::Record::SortAscending());
	modelPartition(pTree, bestRightData, bleaf, bindex, level+1, bestRightMBR, tupleSet2, queryRegions, candidateCutPos);
	ExternalSorter::Record* r_right = tupleSet2.back();
	tupleSet2.pop_back();

	std::vector<ExternalSorter::Record*> parent;
	parent.push_back(r_left);
	parent.push_back(r_right);

	Node* n_parent = createNode(pTree, parent, pTree->m_stats.m_u32TreeHeight - level);
	parent.clear();
	pTree->writeNode(n_parent);

	tupleSet2.push_back(new ExternalSorter::Record(n_parent->m_nodeMBR, n_parent->m_identifier, 0, nullptr, 0));
	pTree->m_rootID = n_parent->m_identifier;

	delete n_parent;

	return can_split;
}

std::vector<int> BulkLoader::float_to_bit_array(float f) 
{
	// Convert float to a 32-bit binary representation
	uint32_t as_int;
	std::memcpy(&as_int, &f, sizeof(float));  // Use memcpy to avoid breaking strict-aliasing rules

	// Convert the integer to a 32-bit binary string
	std::bitset<32> bits(as_int);

	// Convert the binary string to a vector of integers (0s and 1s)
	std::vector<int> bit_array;
	for (int i = 31; i >= 0; --i) {
		bit_array.push_back(bits[i]);
	}
	return bit_array;
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
