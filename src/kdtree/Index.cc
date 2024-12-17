/******************************************************************************
 * Project:  libspatialindex - A C++ library for spatial indexing
 ******************************************************************************
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

#include <limits>

#include <spatialindex/SpatialIndex.h>
#include "KDTree.h"
#include "Node.h"
#include "Leaf.h"
#include "Index.h"

using namespace SpatialIndex;
using namespace SpatialIndex::KDTree;

Index::~Index()
= default;

Index::Index(SpatialIndex::KDTree::KDTree* pTree, id_type id, uint32_t level) : Node(pTree, id, level, pTree->m_indexCapacity)
{
}

NodePtr Index::chooseSubtree(const Region& mbr, uint32_t insertionLevel, std::stack<id_type>& pathBuffer)
{
	if (m_level == insertionLevel) return NodePtr(this, &(m_pTree->m_indexPool));

	pathBuffer.push(m_identifier);

	uint32_t child = 0;

	switch (m_pTree->m_treeVariant)
	{
		case KD_NORMAL:
		case KD_GREEDY:
		case QD_NORMAL:
			child = findInsertChild(mbr);
			break;
		default:
			throw Tools::NotSupportedException("Index::chooseSubtree: Tree variant not supported.");
	}
	assert(child != std::numeric_limits<uint32_t>::max());

	NodePtr n = m_pTree->readNode(m_pIdentifier[child]);
	NodePtr ret = n->chooseSubtree(mbr, insertionLevel, pathBuffer);
	assert(n.unique());
	if (ret.get() == n.get()) n.relinquish();

	return ret;
}

uint32_t Index::findInsertChild(const Region& r) const
{

	RegionPtr t = m_pTree->m_regionPool.acquire();

	for (uint32_t cChild = 0; cChild < m_children; ++cChild)
	{
		if (m_ptrMBR[cChild]->intersectsRegion(r) || m_ptrMBR[cChild]->touchesRegion(r))
		{
			return cChild;
		}
	}

	return 0;
}

NodePtr Index::findLeaf(const Region& mbr, id_type id, std::stack<id_type>& pathBuffer)
{
	pathBuffer.push(m_identifier);

	for (uint32_t cChild = 0; cChild < m_children; ++cChild)
	{
		if (m_ptrMBR[cChild]->containsRegion(mbr))
		{
			NodePtr n = m_pTree->readNode(m_pIdentifier[cChild]);
			NodePtr l = n->findLeaf(mbr, id, pathBuffer);
			if (n.get() == l.get()) n.relinquish();
			if (l.get() != nullptr) return l;
		}
	}

	pathBuffer.pop();

	return NodePtr();
}

void Index::split(uint32_t dataLength, uint8_t* pData, Region& mbr, id_type id, NodePtr& ptrLeft, NodePtr& ptrRight)
{
	++(m_pTree->m_stats.m_u64Splits);

	// std::vector<ExternalSorter::Record*> es;

	// uint32_t sort_dim_index = m_level % m_pTree->m_dimension;


	// for(uint64_t i = 0; i < m_capacity; ++i)
	// {
	// 	es.push_back(new ExternalSorter::Record(m_ptrMBR[i], m_pIdentifier[i], m_pDataLength[i], m_pData[i], sort_dim_index));
	// }

	// es.push_back(new ExternalSorter::Record(mbr, id, dataLength,pData, sort_dim_index));

	// std::sort(es.begin(), es.end(), ExternalSorter::Record::SortAscending());

	// // sort_dim_index = (m_level + 1) % m_pTree->m_dimension;
	// uint64_t total_entries = es.size();
    // uint64_t left_node_en = static_cast<uint64_t>(total_entries / 2);
    // uint64_t right_node_en = static_cast<uint64_t>(total_entries - left_node_en);
	
	// std::vector<ExternalSorter::Record *> left_node_es;
	// std::vector<ExternalSorter::Record *> right_node_es;
	// ExternalSorter::Record* pR_left;
	// ExternalSorter::Record* pR_right;

	// Node* left = new Leaf(m_pTree, -1);

	// for (uint64_t i = 0; i < left_node_en; ++i) {
	// 	try { pR_left = es[i]; }
	// 	catch (Tools::EndOfStreamException&) { break; }
	// 	// pR_left->m_s = sort_dim_index;
	// 	// left_node_es.push_back(pR_left);
	// 	left->insertEntry(pR_left->m_len, pR_left->m_pData, pR_left->m_r, pR_left->m_id);
    // }
	// left->m_level = m_level + 1;
	// left->m_identifier = -1;
	// m_pTree->writeNode(left);

	// Node* right = new Leaf(m_pTree, -1);
	// for (uint64_t i = left_node_en; i < total_entries; ++i) {
	// 	try { pR_right = es[i]; }
	// 	catch (Tools::EndOfStreamException&) { break; }
	// 	// pR_right->m_s = sort_dim_index;
	// 	// right_node_es.push_back(pR_right);
	// 	// right_node_es->insert(pR_right);
	// 	right->insertEntry(pR_right->m_len, pR_right->m_pData, pR_right->m_r, pR_right->m_id);
	// }
	// right->m_level = m_level + 1;
	// right->m_identifier = -1;
	// m_pTree->writeNode(right);

	// pTree->m_stats.m_u32TreeHeight += 1;

	// std::vector<ExternalSorter::Record*> parent;
	// parent.push_back(r_left);
	// parent.push_back(r_right);
	// if (pathBuffer.empty())
	// {

	// 	NodePtr ptrR = m_pTree->m_indexPool.acquire();
	// 	if (ptrR.get() == nullptr)
	// 	{
	// 		ptrR = NodePtr(new Index(m_pTree, m_pTree->m_rootID, m_level + 1), &(m_pTree->m_indexPool));
	// 	}
	// 	else
	// 	{
	// 		//ptrR->m_pTree = m_pTree;
	// 		ptrR->m_identifier = m_pTree->m_rootID;
	// 		ptrR->m_level = m_level + 1;
	// 		ptrR->m_nodeMBR = m_pTree->m_infiniteRegion;
	// 	}

	// 	m_capacity = 2;
	// 	m_children = 0;

	// 	ptrR->insertEntry(0, nullptr, left->m_nodeMBR, left->m_identifier);
	// 	ptrR->insertEntry(0, nullptr, right->m_nodeMBR, right->m_identifier);

	// 	// m_pTree->writeNode(ptrR.get());

	// 	// m_pTree->m_stats.m_nodesInLevel[m_level] = 2;
	// 	// m_pTree->m_stats.m_nodesInLevel.push_back(1);
	// 	// m_pTree->m_stats.m_u32TreeHeight = m_level + 2;
	// }
	// // else
	// // {


	// 	// id_type cParent = pathBuffer.top(); pathBuffer.pop();
	// 	// NodePtr ptrN = m_pTree->readNode(cParent);
	// 	// Index* p = static_cast<Index*>(ptrN.get());
	// 	// p->adjustTree(n.get(), nn.get(), pathBuffer, overflowTable);
	// // }

}

uint32_t Index::findLeastEnlargement(const Region& r) const
{
	double area = std::numeric_limits<double>::infinity();
	uint32_t best = std::numeric_limits<uint32_t>::max();

	RegionPtr t = m_pTree->m_regionPool.acquire();

	for (uint32_t cChild = 0; cChild < m_children; ++cChild)
	{
		m_ptrMBR[cChild]->getCombinedRegion(*t, r);

		double a = m_ptrMBR[cChild]->getArea();
		double enl = t->getArea() - a;

		if (enl < area)
		{
			area = enl;
			best = cChild;
		}
		else if (enl == area)
		{
			// this will rarely happen, so compute best area on the fly only
			// when necessary.
			if (enl == std::numeric_limits<double>::infinity()
			    || a < m_ptrMBR[best]->getArea()) best = cChild;
		}
	}

	return best;
}

uint32_t Index::findLeastOverlap(const Region& r) const
{
	OverlapEntry** entries = new OverlapEntry*[m_children];

	double leastOverlap = std::numeric_limits<double>::max();
	double me = std::numeric_limits<double>::max();
	OverlapEntry* best = nullptr;

	// find combined region and enlargement of every entry and store it.
	for (uint32_t cChild = 0; cChild < m_children; ++cChild)
	{
		try
		{
			entries[cChild] = new OverlapEntry();
		}
		catch (...)
		{
			for (uint32_t i = 0; i < cChild; ++i) delete entries[i];
			delete[] entries;
			throw;
		}

		entries[cChild]->m_index = cChild;
		entries[cChild]->m_original = m_ptrMBR[cChild];
		entries[cChild]->m_combined = m_pTree->m_regionPool.acquire();
		m_ptrMBR[cChild]->getCombinedRegion(*(entries[cChild]->m_combined), r);
		entries[cChild]->m_oa = entries[cChild]->m_original->getArea();
		entries[cChild]->m_ca = entries[cChild]->m_combined->getArea();
		entries[cChild]->m_enlargement = entries[cChild]->m_ca - entries[cChild]->m_oa;

		if (entries[cChild]->m_enlargement < me)
		{
			me = entries[cChild]->m_enlargement;
			best = entries[cChild];
		}
		else if (entries[cChild]->m_enlargement == me && entries[cChild]->m_oa < best->m_oa)
		{
			best = entries[cChild];
		}
	}

	if (me < -std::numeric_limits<double>::epsilon() || me > std::numeric_limits<double>::epsilon())
	{
		uint32_t cIterations;

		if (m_children > m_pTree->m_nearMinimumOverlapFactor)
		{
			// sort entries in increasing order of enlargement.
			::qsort(entries, m_children,
							sizeof(OverlapEntry*),
							OverlapEntry::compareEntries);
			assert(entries[0]->m_enlargement <= entries[m_children - 1]->m_enlargement);

			cIterations = m_pTree->m_nearMinimumOverlapFactor;
		}
		else
		{
			cIterations = m_children;
		}

		// calculate overlap of most important original entries (near minimum overlap cost).
		for (uint32_t cIndex = 0; cIndex < cIterations; ++cIndex)
		{
			double dif = 0.0;
			OverlapEntry* e = entries[cIndex];

			for (uint32_t cChild = 0; cChild < m_children; ++cChild)
			{
				if (e->m_index != cChild)
				{
					double f = e->m_combined->getIntersectingArea(*(m_ptrMBR[cChild]));
					if (f != 0.0) dif += f - e->m_original->getIntersectingArea(*(m_ptrMBR[cChild]));
				}
			} // for (cChild)

			if (dif < leastOverlap)
			{
				leastOverlap = dif;
				best = entries[cIndex];
			}
			else if (dif == leastOverlap)
			{
				if (e->m_enlargement == best->m_enlargement)
				{
					// keep the one with least area.
					if (e->m_original->getArea() < best->m_original->getArea()) best = entries[cIndex];
				}
				else
				{
					// keep the one with least enlargement.
					if (e->m_enlargement < best->m_enlargement) best = entries[cIndex];
				}
			}
		} // for (cIndex)
	}

	uint32_t ret = best->m_index;

	for (uint32_t cChild = 0; cChild < m_children; ++cChild)
	{
		delete entries[cChild];
	}
	delete[] entries;

	return ret;
}

void Index::adjustTree(Node* n, std::stack<id_type>& pathBuffer, bool force)
{
	++(m_pTree->m_stats.m_u64Adjustments);

	// find entry pointing to old node;
	uint32_t child;
	for (child = 0; child < m_children; ++child)
	{
		if (m_pIdentifier[child] == n->m_identifier) break;
	}

	// MBR needs recalculation if either:
	//   1. the NEW child MBR is not contained.
	//   2. the OLD child MBR is touching.
	bool bContained = m_nodeMBR.containsRegion(n->m_nodeMBR);
	bool bTouches = m_nodeMBR.touchesRegion(*(m_ptrMBR[child]));
	bool bRecompute = (! bContained || (bTouches && m_pTree->m_bTightMBRs));

	*(m_ptrMBR[child]) = n->m_nodeMBR;

	if (bRecompute || force)
	{
		for (uint32_t cDim = 0; cDim < m_nodeMBR.m_dimension; ++cDim)
		{
			m_nodeMBR.m_pLow[cDim] = std::numeric_limits<double>::max();
			m_nodeMBR.m_pHigh[cDim] = -std::numeric_limits<double>::max();

			for (uint32_t cChild = 0; cChild < m_children; ++cChild)
			{
				m_nodeMBR.m_pLow[cDim] = std::min(m_nodeMBR.m_pLow[cDim], m_ptrMBR[cChild]->m_pLow[cDim]);
				m_nodeMBR.m_pHigh[cDim] = std::max(m_nodeMBR.m_pHigh[cDim], m_ptrMBR[cChild]->m_pHigh[cDim]);
			}
		}
	}

	m_pTree->writeNode(this);

	if ((bRecompute || force) && (! pathBuffer.empty()))
	{
		id_type cParent = pathBuffer.top(); pathBuffer.pop();
		NodePtr ptrN = m_pTree->readNode(cParent);
		Index* p = static_cast<Index*>(ptrN.get());
		p->adjustTree(this, pathBuffer, force);
	}
}

void Index::adjustTree(Node* n1, Node* n2, std::stack<id_type>& pathBuffer, uint8_t* overflowTable)
{
	++(m_pTree->m_stats.m_u64Adjustments);

	// find entry pointing to old node;
	uint32_t child;
	for (child = 0; child < m_children; ++child)
	{
		if (m_pIdentifier[child] == n1->m_identifier) break;
	}

	// MBR needs recalculation if either:
	//   1. either child MBR is not contained.
	//   2. the OLD child MBR is touching.
	bool bContained1 = m_nodeMBR.containsRegion(n1->m_nodeMBR);
	bool bContained2 = m_nodeMBR.containsRegion(n2->m_nodeMBR);
	bool bContained = bContained1 && bContained2;
	bool bTouches = m_nodeMBR.touchesRegion(*(m_ptrMBR[child]));
	bool bRecompute = (! bContained || (bTouches && m_pTree->m_bTightMBRs));

	*(m_ptrMBR[child]) = n1->m_nodeMBR;

	if (bRecompute)
	{
		for (uint32_t cDim = 0; cDim < m_nodeMBR.m_dimension; ++cDim)
		{
			m_nodeMBR.m_pLow[cDim] = std::numeric_limits<double>::max();
			m_nodeMBR.m_pHigh[cDim] = -std::numeric_limits<double>::max();

			for (uint32_t cChild = 0; cChild < m_children; ++cChild)
			{
				m_nodeMBR.m_pLow[cDim] = std::min(m_nodeMBR.m_pLow[cDim], m_ptrMBR[cChild]->m_pLow[cDim]);
				m_nodeMBR.m_pHigh[cDim] = std::max(m_nodeMBR.m_pHigh[cDim], m_ptrMBR[cChild]->m_pHigh[cDim]);
			}
		}
	}

	// No write necessary here. insertData will write the node if needed.
	//m_pTree->writeNode(this);

	bool bAdjusted = insertData(0, nullptr, n2->m_nodeMBR, n2->m_identifier, pathBuffer, overflowTable);

	// if n2 is contained in the node and there was no split or reinsert,
	// we need to adjust only if recalculation took place.
	// In all other cases insertData above took care of adjustment.
	if ((! bAdjusted) && bRecompute && (! pathBuffer.empty()))
	{
		id_type cParent = pathBuffer.top(); pathBuffer.pop();
		NodePtr ptrN = m_pTree->readNode(cParent);
		Index* p = static_cast<Index*>(ptrN.get());
		p->adjustTree(this, pathBuffer);
	}
}
