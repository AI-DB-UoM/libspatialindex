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

#include <cstring>

#include <spatialindex/SpatialIndex.h>

#include "KDTree.h"
#include "Node.h"
#include "Index.h"
#include "Leaf.h"

using namespace SpatialIndex;
using namespace SpatialIndex::KDTree;

Leaf::~Leaf()
= default;

Leaf::Leaf(SpatialIndex::KDTree::KDTree* pTree, id_type id): Node(pTree, id, 0, pTree->m_leafCapacity)
{
}

NodePtr Leaf::chooseSubtree(const Region&, uint32_t, std::stack<id_type>&)
{
	// should make sure to relinquish other PoolPointer lists that might be pointing to the
	// same leaf.
	return NodePtr(this, &(m_pTree->m_leafPool));
}

NodePtr Leaf::findLeaf(const Region& mbr, id_type id, std::stack<id_type>&)
{
	for (uint32_t cChild = 0; cChild < m_children; ++cChild)
	{
		// should make sure to relinquish other PoolPointer lists that might be pointing to the
		// same leaf.
		if (m_pIdentifier[cChild] == id && mbr == *(m_ptrMBR[cChild])) return NodePtr(this, &(m_pTree->m_leafPool));
	}

	return NodePtr();
}

void Leaf::split(uint32_t dataLength, uint8_t* pData, Region& mbr, id_type id, NodePtr& pLeft, NodePtr& pRight)
{
	++(m_pTree->m_stats.m_u64Splits);

	RstarSplitEntry** dataLow = nullptr;
	try
	{
		dataLow = new RstarSplitEntry*[m_capacity + 1];
	}
	catch (...)
	{
		delete[] dataLow;
		throw;
	}

	m_pDataLength[m_capacity] = dataLength;
	m_pData[m_capacity] = pData;
	m_ptrMBR[m_capacity] = m_pTree->m_regionPool.acquire();
	*(m_ptrMBR[m_capacity]) = mbr;
	m_pIdentifier[m_capacity] = id;

	uint32_t u32Child = 0, cDim, cIndex;
	uint32_t sort_dim_index = m_level % m_pTree->m_dimension;


	for (u32Child = 0; u32Child <= m_capacity; ++u32Child)
	{
		try
		{
			dataLow[u32Child] = new RstarSplitEntry(m_ptrMBR[u32Child].get(), u32Child, sort_dim_index);
		}
		catch (...)
		{
			for (uint32_t i = 0; i < u32Child; ++i) delete dataLow[i];
			delete[] dataLow;
			throw;
		}
	}

	::qsort(dataLow, m_capacity + 1, sizeof(RstarSplitEntry*), RstarSplitEntry::compareLow);

	uint32_t total_entries = m_capacity + 1;
    uint32_t left_node_en = static_cast<uint32_t>(total_entries / 2);
    uint32_t right_node_en = static_cast<uint32_t>(total_entries - left_node_en);

	pLeft = NodePtr(new Leaf(m_pTree, -1), &(m_pTree->m_leafPool));
	for (uint32_t i = 0; i < left_node_en; ++i) {
		try 
		{ 
			uint32_t child_id = dataLow[i]->m_index;
			pLeft->insertEntry(m_pDataLength[child_id], m_pData[child_id], *(m_ptrMBR[child_id]), m_pIdentifier[child_id]);
		}
		catch (Tools::EndOfStreamException&) { break; }
    }
	// m_pTree->writeNode(pLeft.get());

	pRight = NodePtr(new Leaf(m_pTree, -1), &(m_pTree->m_leafPool));
	for (uint32_t i = left_node_en; i < total_entries; ++i) {
		try 
		{ 
			uint32_t child_id = dataLow[i]->m_index;
			pRight->insertEntry(m_pDataLength[child_id], m_pData[child_id], *(m_ptrMBR[child_id]), m_pIdentifier[child_id]);
		}
		catch (Tools::EndOfStreamException&) { break; }
		
    }

}

void Leaf::deleteData(const Region& mbr, id_type id, std::stack<id_type>& pathBuffer)
{
	uint32_t child;

	for (child = 0; child < m_children; ++child)
	{
		if (m_pIdentifier[child] == id && mbr == *(m_ptrMBR[child])) break;
	}

	deleteEntry(child);
	m_pTree->writeNode(this);

	std::stack<NodePtr> toReinsert;
	NodePtr ptrThis(this, &(m_pTree->m_leafPool));
	condenseTree(toReinsert, pathBuffer, ptrThis);
	ptrThis.relinquish();

	// re-insert eliminated nodes.
	while (! toReinsert.empty())
	{
		NodePtr n = toReinsert.top(); toReinsert.pop();
		m_pTree->deleteNode(n.get());

		for (uint32_t cChild = 0; cChild < n->m_children; ++cChild)
		{
			// keep this in the for loop. The tree height might change after insertions.
			uint8_t* overflowTable = new uint8_t[m_pTree->m_stats.m_u32TreeHeight];
			memset(overflowTable, 0, m_pTree->m_stats.m_u32TreeHeight);
			m_pTree->insertData_impl(n->m_pDataLength[cChild], n->m_pData[cChild], *(n->m_ptrMBR[cChild]), n->m_pIdentifier[cChild], n->m_level, overflowTable);
			n->m_pData[cChild] = nullptr;
			delete[] overflowTable;
		}
		if (n.get() == this) n.relinquish();
	}
}
