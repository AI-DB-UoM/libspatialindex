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
#include <cmath>
#include <limits>

#include <spatialindex/SpatialIndex.h>

#include "LearnedIndex.h"
#include "Node.h"
#include "Index.h"
#include "Leaf.h"

// #include <torch/script.h>

using namespace SpatialIndex;
using namespace SpatialIndex::LearnedIndex;

//
// Tools::IObject interface
//
Tools::IObject* Node::clone()
{
	throw Tools::NotSupportedException("IObject::clone should never be called.");
}

//
// Tools::ISerializable interface
//
uint32_t Node::getByteArraySize()
{
	uint32_t model_size = 0;
	if (m_level > 0) {
		model_size = sizeof(uint32_t);
	}
	if (m_modelDataLength > 0) {
		model_size += m_modelDataLength * sizeof(double);
	}
	return
		(sizeof(uint32_t) +
		sizeof(uint32_t) +
		sizeof(uint32_t) +
		(m_children * (m_pTree->m_dimension * sizeof(double) * 2 + sizeof(id_type) + sizeof(uint32_t))) +
		m_totalDataLength +
		model_size +  // for m_modelDataLength
		(2 * m_pTree->m_dimension * sizeof(double)));

}

void Node::loadFromByteArray(const uint8_t* ptr)
{
	m_nodeMBR = m_pTree->m_infiniteRegion;

	// skip the node type information, it is not needed.
	ptr += sizeof(uint32_t);

	memcpy(&m_level, ptr, sizeof(uint32_t));
	ptr += sizeof(uint32_t);

	memcpy(&m_children, ptr, sizeof(uint32_t));
	ptr += sizeof(uint32_t);

	for (uint32_t u32Child = 0; u32Child < m_children; ++u32Child)
	{
		m_ptrMBR[u32Child] = m_pTree->m_regionPool.acquire();
		*(m_ptrMBR[u32Child]) = m_pTree->m_infiniteRegion;

		memcpy(m_ptrMBR[u32Child]->m_pLow, ptr, m_pTree->m_dimension * sizeof(double));
		ptr += m_pTree->m_dimension * sizeof(double);
		memcpy(m_ptrMBR[u32Child]->m_pHigh, ptr, m_pTree->m_dimension * sizeof(double));
		ptr += m_pTree->m_dimension * sizeof(double);
		memcpy(&(m_pIdentifier[u32Child]), ptr, sizeof(id_type));
		ptr += sizeof(id_type);

		memcpy(&(m_pDataLength[u32Child]), ptr, sizeof(uint32_t));
		ptr += sizeof(uint32_t);

		if (m_pDataLength[u32Child] > 0)
		{
			m_totalDataLength += m_pDataLength[u32Child];
			m_pData[u32Child] = new uint8_t[m_pDataLength[u32Child]];
			memcpy(m_pData[u32Child], ptr, m_pDataLength[u32Child]);
			ptr += m_pDataLength[u32Child];
		}
		else
		{
			m_pData[u32Child] = nullptr;
		}

		//m_nodeMBR.combineRegion(*(m_ptrMBR[u32Child]));
	}

	memcpy(m_nodeMBR.m_pLow, ptr, m_pTree->m_dimension * sizeof(double));
	ptr += m_pTree->m_dimension * sizeof(double);
	memcpy(m_nodeMBR.m_pHigh, ptr, m_pTree->m_dimension * sizeof(double));
	ptr += m_pTree->m_dimension * sizeof(double);
	
	if (m_level > 0) {
		memcpy(&m_modelDataLength, ptr, sizeof(uint32_t));
		if (m_modelDataLength > 0) {
			ptr += sizeof(uint32_t);
			memcpy(m_modelData, ptr, m_modelDataLength * sizeof(double));
		}
	}

	
}

void Node::storeToByteArray(uint8_t** data, uint32_t& len)
{
	len = getByteArraySize();

	// std::cerr << "------------ storeToByteArray len------------" << len << std::endl;
	// std::cerr << "------------ m_level------------" << m_level << std::endl;

	*data = new uint8_t[len];
	uint8_t* ptr = *data;

	uint32_t nodeType;

	if (m_level == 0) nodeType = PersistentLeaf;
	else nodeType = PersistentIndex;

	memcpy(ptr, &nodeType, sizeof(uint32_t));
	ptr += sizeof(uint32_t);

	memcpy(ptr, &m_level, sizeof(uint32_t));
	ptr += sizeof(uint32_t);

	memcpy(ptr, &m_children, sizeof(uint32_t));
	ptr += sizeof(uint32_t);

	for (uint32_t u32Child = 0; u32Child < m_children; ++u32Child)
	{
		memcpy(ptr, m_ptrMBR[u32Child]->m_pLow, m_pTree->m_dimension * sizeof(double));
		ptr += m_pTree->m_dimension * sizeof(double);
		memcpy(ptr, m_ptrMBR[u32Child]->m_pHigh, m_pTree->m_dimension * sizeof(double));
		ptr += m_pTree->m_dimension * sizeof(double);
		memcpy(ptr, &(m_pIdentifier[u32Child]), sizeof(id_type));
		ptr += sizeof(id_type);

		memcpy(ptr, &(m_pDataLength[u32Child]), sizeof(uint32_t));
		ptr += sizeof(uint32_t);

		if (m_pDataLength[u32Child] > 0)
		{
			memcpy(ptr, m_pData[u32Child], m_pDataLength[u32Child]);
			ptr += m_pDataLength[u32Child];
		}
	}

	// store the node MBR for efficiency. This increases the node size a little bit.
	memcpy(ptr, m_nodeMBR.m_pLow, m_pTree->m_dimension * sizeof(double));
	ptr += m_pTree->m_dimension * sizeof(double);
	memcpy(ptr, m_nodeMBR.m_pHigh, m_pTree->m_dimension * sizeof(double));
	ptr += m_pTree->m_dimension * sizeof(double);

	// std::cerr << "------------ m_modelDataLength------------" << m_modelDataLength << std::endl;
	
	if (m_level > 0) {
		
		memcpy(ptr, &m_modelDataLength, sizeof(uint32_t));
		ptr += sizeof(uint32_t);
		memcpy(ptr, m_modelData, m_modelDataLength * sizeof(double));
		ptr += m_modelDataLength * sizeof(double);
	}
	// std::cerr << " len:" << len << " (ptr - *data):" << (ptr - *data) << std::endl;

	assert(len == (ptr - *data));

}

//
// SpatialIndex::IEntry interface
//
SpatialIndex::id_type Node::getIdentifier() const
{
	return m_identifier;
}

void Node::getShape(IShape** out) const
{
	*out = new Region(m_nodeMBR);
}

//
// SpatialIndex::INode interface
//
uint32_t Node::getChildrenCount() const
{
	return m_children;
}

SpatialIndex::id_type Node::getChildIdentifier(uint32_t index) const
{
	if (index >= m_children) throw Tools::IndexOutOfBoundsException(index);

	return m_pIdentifier[index];
}

void Node::getChildShape(uint32_t index, IShape** out) const
{
	if (index >= m_children) throw Tools::IndexOutOfBoundsException(index);

	*out = new Region(*(m_ptrMBR[index]));
}

void Node::getChildData(uint32_t index, uint32_t& length, uint8_t** data) const
{
	if (index >= m_children) throw Tools::IndexOutOfBoundsException(index);
	if (m_pData[index] == nullptr)
	{
		length = 0;
		data = nullptr;
	}
	else
	{
		length = m_pDataLength[index];
		*data = m_pData[index];
	}
}

uint32_t Node::getLevel() const
{
	return m_level;
}

bool Node::isLeaf() const
{
	return (m_level == 0);
}

bool Node::isIndex() const
{
	return (m_level != 0);
}

//
// Internal
//

Node::Node()
= default;

Node::Node(SpatialIndex::LearnedIndex::LearnedIndex* pTree, id_type id, uint32_t level, uint32_t capacity) :
	m_pTree(pTree),
	m_level(level),
	m_identifier(id),
	m_children(0),
	m_capacity(capacity),
	m_pData(nullptr),
	m_ptrMBR(nullptr),
	m_pIdentifier(nullptr),
	m_pDataLength(nullptr),
	m_totalDataLength(0),
	m_modelData(nullptr),
	m_modelDataLength(0)
{
	m_nodeMBR.makeInfinite(m_pTree->m_dimension);

	try
	{
		m_pDataLength = new uint32_t[m_capacity + 1];
		m_pData = new uint8_t*[m_capacity + 1];
		m_ptrMBR = new RegionPtr[m_capacity + 1];
		m_pIdentifier = new id_type[m_capacity + 1];
		m_modelData = new double[2]; // Only need two
	}
	catch (...)
	{
		delete[] m_pDataLength;
		delete[] m_pData;
		delete[] m_ptrMBR;
		delete[] m_pIdentifier;
		delete[] m_modelData;
		throw;
	}
}

Node::~Node()
{
	if (m_pData != nullptr)
	{
		for (uint32_t u32Child = 0; u32Child < m_children; ++u32Child)
		{
			if (m_pData[u32Child] != nullptr) delete[] m_pData[u32Child];
		}

		delete[] m_pData;
	}

	delete[] m_pDataLength;
	delete[] m_ptrMBR;
	delete[] m_pIdentifier;
	delete[] m_modelData;
}

Node& Node::operator=(const Node&)
{
	throw Tools::IllegalStateException("operator =: This should never be called.");
}

void Node::insertEntry(uint32_t dataLength, uint8_t* pData, Region& mbr, id_type id)
{
	assert(m_children < m_capacity);
	m_pDataLength[m_children] = dataLength;
	m_pData[m_children] = pData;
	m_ptrMBR[m_children] = m_pTree->m_regionPool.acquire();
	*(m_ptrMBR[m_children]) = mbr;
	m_pIdentifier[m_children] = id;
	m_totalDataLength += dataLength;
	++m_children;
	m_nodeMBR.combineRegion(mbr);
}

// void Node::insertEntry(uint32_t dataLength, uint8_t* pData, Region& mbr, id_type id)
// {
// 	std::cerr << ">> insertEntry called\n";
// 	std::cerr << "  m_children = " << m_children << ", m_capacity = " << m_capacity << "\n";
// 	std::cerr << "  dataLength = " << dataLength << ", id = " << id << "\n";
// 	std::cerr << "  pData ptr  = " << static_cast<void*>(pData) << "\n";
// 	std::cerr << "  m_pTree    = " << m_pTree << "\n";

// 	if (m_children >= m_capacity) {
// 		std::cerr << "[Error] Node over capacity!\n";
// 		std::abort();
// 	}

// 	if (!m_pDataLength || !m_pData || !m_ptrMBR || !m_pIdentifier) {
// 		std::cerr << "[Error] One or more internal arrays are null!\n";
// 		std::abort();
// 	}

// 	m_pDataLength[m_children] = dataLength;
// 	m_pData[m_children] = pData;

// 	// m_ptrMBR[m_children] = m_pTree->m_regionPool.acquire();
// 	// if (!m_pTree->m_regionPool.acquire()) {
// 	// 	std::cerr << "[Error] RegionPool returned nullptr!\n";
// 	// 	std::abort();
// 	// }
// 	// Region* acquired = m_pTree->m_regionPool.acquire();

// 	std::cerr << "  Region acquired at " << acquired << ", assigning value...\n";
// 	*acquired = mbr;  // Potential crash here

// 	m_pIdentifier[m_children] = id;

// 	m_totalDataLength += dataLength;
// 	++m_children;

// 	std::cerr << "  Total data length now = " << m_totalDataLength << "\n";

// 	std::cerr << "  Updating m_nodeMBR...\n";
// 	m_nodeMBR.combineRegion(mbr);  // Optional crash point if mbr has issues

// 	std::cerr << "<< insertEntry done\n";
// }


void Node::insertModel(const std::vector<double>& newModelParams)
{
    delete[] m_modelData;
    m_modelData = nullptr;

    m_modelDataLength = newModelParams.size();

    if (m_modelDataLength > 0)
    {
        m_modelData = new double[m_modelDataLength];
        std::copy(newModelParams.begin(), newModelParams.end(), m_modelData);
    }
}

void Node::deleteEntry(uint32_t index)
{
	assert(index >= 0 && index < m_children);

	// cache it, since I might need it for "touches" later.
	RegionPtr ptrR = m_ptrMBR[index];

	m_totalDataLength -= m_pDataLength[index];
	if (m_pData[index] != nullptr) delete[] m_pData[index];

	if (m_children > 1 && index != m_children - 1)
	{
		m_pDataLength[index] = m_pDataLength[m_children - 1];
		m_pData[index] = m_pData[m_children - 1];
		m_ptrMBR[index] = m_ptrMBR[m_children - 1];
		m_pIdentifier[index] = m_pIdentifier[m_children - 1];
	}

	--m_children;

	// WARNING: index has now changed. Do not use it below here.

	if (m_children == 0)
	{
		m_nodeMBR = m_pTree->m_infiniteRegion;
	}
	else if (m_pTree->m_bTightMBRs && m_nodeMBR.touchesRegion(*ptrR))
	{
		for (uint32_t cDim = 0; cDim < m_nodeMBR.m_dimension; ++cDim)
		{
			m_nodeMBR.m_pLow[cDim] = std::numeric_limits<double>::max();
			m_nodeMBR.m_pHigh[cDim] = -std::numeric_limits<double>::max();

			for (uint32_t u32Child = 0; u32Child < m_children; ++u32Child)
			{
				m_nodeMBR.m_pLow[cDim] = std::min(m_nodeMBR.m_pLow[cDim], m_ptrMBR[u32Child]->m_pLow[cDim]);
				m_nodeMBR.m_pHigh[cDim] = std::max(m_nodeMBR.m_pHigh[cDim], m_ptrMBR[u32Child]->m_pHigh[cDim]);
			}
		}
	}
}

bool Node::insertData(uint32_t dataLength, uint8_t* pData, Region& mbr, id_type id, std::stack<id_type>& pathBuffer, uint8_t* overflowTable)
{
	if (m_children >= m_capacity)
    {
		++(m_pTree->m_stats.m_u64Splits);

        // Copy all current + new entry into dataMid
        std::vector<std::unique_ptr<RstarSplitEntry>> dataMid;
        dataMid.reserve(m_capacity + 1);

        for (uint32_t i = 0; i < m_children; ++i) {
            dataMid.emplace_back(std::make_unique<RstarSplitEntry>(
                m_ptrMBR[i].get(), i, m_pDataLength[i], m_pData[i], 1));
        }

        RegionPtr newRegion = m_pTree->m_regionPool.acquire();
        *newRegion = mbr;
        dataMid.emplace_back(std::make_unique<RstarSplitEntry>(
            newRegion.get(), m_capacity, dataLength, pData, 1));

        // Temporarily store the new entry
        std::vector<double> keys;
        keys.reserve(dataMid.size());

        for (const auto& entry : dataMid) {
            auto* r = entry->m_pRegion;
            keys.push_back((r->m_pLow[1] + r->m_pHigh[1]) / 2.0);
        }

        // Fit linear regression
        double mean_x = 0.0, mean_y = 0.0, covariance = 0.0, variance = 0.0;
        uint64_t n = 0;
        double scale = 1.0 * m_capacity / keys.size();

        for (size_t i = 0; i < keys.size(); ++i) {
            double x = keys[i];
            double y = static_cast<double>(i * scale);

            ++n;
            double dx = x - mean_x;
            mean_x += dx / static_cast<double>(n);
            mean_y += (y - mean_y) / static_cast<double>(n);
            covariance += dx * (y - mean_y);
            variance += dx * (x - mean_x);
        }

        covariance /= static_cast<double>(n - 1);
        variance /= static_cast<double>(n - 1);

        double slope = (variance == 0.0) ? 0.0 : covariance / variance;
        double intercept = mean_y - slope * mean_x;

        std::vector<std::vector<RstarSplitEntry*>> children_se(m_capacity);

        for (const auto& entry : dataMid) {
            double x = (entry->m_pRegion->m_pLow[1] + entry->m_pRegion->m_pHigh[1]) / 2.0;
            double y = slope * x + intercept;
            uint32_t predict_index = std::clamp(static_cast<uint32_t>(y), 0U, m_capacity - 1);
            children_se[predict_index].push_back(entry.get());
        }

        m_children = 0;
        for (uint32_t j = 0; j < m_capacity; ++j) {
            if (children_se[j].empty()) continue;

            NodePtr n(new Leaf(m_pTree, -1), &m_pTree->m_leafPool);

            for (auto* entry : children_se[j]) {
                Region r = *entry->m_pRegion;
                uint8_t* copied = new uint8_t[entry->m_len];
                std::memcpy(copied, entry->m_pData, entry->m_len);
                n->insertEntry(entry->m_len, copied, r, entry->m_id);
            }

            m_pTree->writeNode(n.get());
            insertEntry(0, nullptr, n->m_nodeMBR, n->m_identifier);
        }

        insertModel({slope, intercept});
        m_pTree->writeNode(this);
    }
    else
    {
        insertEntry(dataLength, pData, mbr, id);
        m_pTree->writeNode(this);
    }
    return true;
	// if (m_children >= m_capacity)
	// {
	// 	std::cerr << " insertData m_children >= m_capacity" << std::endl;

	// 	// use all points to train a model
	// 	RstarSplitEntry** dataMid = nullptr;
	// 	try
	// 	{
	// 		dataMid = new RstarSplitEntry*[m_capacity + 1];
	// 	}
	// 	catch (...)
	// 	{
	// 		throw;
	// 	}

	// 	m_pDataLength[m_capacity] = dataLength;
	// 	m_pData[m_capacity] = pData;
	// 	m_ptrMBR[m_capacity] = m_pTree->m_regionPool.acquire();
	// 	*(m_ptrMBR[m_capacity]) = mbr;
	// 	m_pIdentifier[m_capacity] = id;

	// 	uint32_t u32Child = 0, cDim, cIndex;

	// 	for (u32Child = 0; u32Child <= m_capacity; ++u32Child)
	// 	{
	// 		try
	// 		{
	// 			// Region* pr, uint32_t index, uint32_t len, uint8_t* pData, uint32_t dimension):
	// 			dataMid[u32Child] = new RstarSplitEntry(m_ptrMBR[u32Child].get(), u32Child, m_pDataLength[u32Child], m_pData[u32Child], 1);
	// 		}
	// 		catch (...)
	// 		{
	// 			for (uint32_t i = 0; i < u32Child; ++i) delete dataMid[i];
	// 			delete[] dataMid;
	// 			throw;
	// 		}
	// 	}

	// 	::qsort(dataMid, m_capacity + 1, sizeof(RstarSplitEntry*), RstarSplitEntry::compareMid);
	// 	// std::cerr << " after sort" << std::endl;
		
	// 	std::vector<double> keys;
	
	// 	for(uint64_t j = 0; j <= m_capacity; j++) {
	// 		double y = (m_ptrMBR[j]->m_pLow[1] + m_ptrMBR[j]->m_pHigh[1]) / 2;
	// 		keys.push_back(y);
	// 	}
	// 	// std::cerr << " keys.size(): " << keys.size() << std::endl;
	// 	double mean_x = 0.0, mean_y = 0.0, covariance = 0.0, variance = 0.0;
	// 	uint64_t n = 0;  // Count of elements
	// 	int data_size = keys.size();
	// 	double scale = 1.0 * m_capacity / data_size;

	// 	for (int i = 0; i < data_size; ++i) {
	// 		double x = keys[i];
	// 		double y = static_cast<double>(i * scale);

	// 		++n;

	// 		// Update means
	// 		double dx = x - mean_x;
	// 		mean_x += dx / static_cast<double>(n);
	// 		mean_y += (y - mean_y) / static_cast<double>(n);

	// 		// Update covariance and variance
	// 		covariance += dx * (y - mean_y);
	// 		variance += dx * (x - mean_x);
	// 	}


	// 	variance /= static_cast<double>(n - 1);
	// 	covariance /= static_cast<double>(n - 1);

	// 	double slope = 0;
	// 	double intercept = 0;
	// 	if (variance == 0.0) {
	// 		slope = mean_y;  // No variation in x
	// 	} else {
	// 		slope = covariance / variance;  // Slope
	// 		intercept = mean_y - slope * mean_x;  // Intercept
	// 	}

	// 	// std::cerr << " slope: " << slope << " intercept: " << intercept << std::endl;

	// 	std::vector<double> regressionParams = {slope, intercept};

	// 	std::vector<std::vector<RstarSplitEntry *>> children_se(m_capacity);
	// 	RstarSplitEntry* pR;

	// 	uint64_t i = 0;
	// 	while (i < m_capacity) {
	// 		try { 
	// 			pR = dataMid[i]; 
	// 		}
	// 		catch (Tools::EndOfStreamException&) { break; }
	// 		double x = (dataMid[i]->m_pRegion->m_pLow[1] + dataMid[i]->m_pRegion->m_pHigh[1]) / 2;
	// 		double y = slope * x + intercept;
	
	// 		uint32_t predict_index = static_cast<uint32_t>(y);
	// 		// std::cerr << "partition[i]->m_r.m_pLow[1]: " << partition[i]->m_r.m_pLow[1] << " partition[i]->m_r.m_pHigh[1]:" << partition[i]->m_r.m_pHigh[1] << std::endl;
	// 		// std::cerr << "prediction a:" << a << " b:" << b << " x:" << x << " y:" << y << " predict_index:" << predict_index << std::endl;
	
	// 		predict_index = std::max(static_cast<uint32_t>(0), std::min(predict_index, static_cast<uint32_t>(m_capacity - 1)));
	// 		// uint64_t predict_index = i / entries_per_child;
	// 		children_se[predict_index].push_back(pR);
	// 		i++;
	// 	}

	// 	// std::cerr << " assign to list: " << std::endl;

	// 	m_children = 0;
	// 	for (uint32_t j = 0; j < m_capacity; ++j) {
	// 		// Node* n = createNode(m_pTree, children_se[j], 0);
	// 		NodePtr n;
	// 		n = NodePtr(new Leaf(m_pTree, -1), &(m_pTree->m_leafPool));

	// 		// std::cerr << "create new leaf n->m_children:" << n->m_children << std::endl;
	// 		// std::cerr << "create new leaf n->m_capacity:" << n->m_capacity << std::endl;
	// 		// std::cerr << "create new leaf children_se[j].size():" << children_se[j].size() << std::endl;

	// 		// Node(SpatialIndex::LearnedIndex::LearnedIndex* pTree, id_type id, uint32_t level, uint32_t capacity) :
	// 		for (size_t cChild = 0; cChild < children_se[j].size(); ++cChild)
	// 		{
	// 			auto* entry = children_se[j][cChild];
	// 			try 
	// 			{ 
	// 				n->insertEntry(entry->m_len, entry->m_pData, *(entry->m_pRegion), entry->m_id);
	// 			// n->insertEntry(children_se[j][cChild]->m_len, children_se[j][cChild]->m_pData, *children_se[j][cChild]->m_pRegion, children_se[j][cChild]->m_id);
	// 			}
	// 			catch (Tools::EndOfStreamException&) { break; }
	// 		}
	// 		// std::cerr << "insertEntry n->m_identifier:" << n->m_identifier << std::endl;

	// 		// children_se[j].clear();
	// 		// std::cerr << "writeNode " << std::endl;

	// 		m_pTree->writeNode(n.get());
	// 		// std::cerr << "writeNode finish" << std::endl;

	// 		// m_pTree->m_rootID = n->m_identifier;
	// 		// m_pTree->m_stats.m_u32TreeHeight = std::max(m_pTree->m_stats.m_u32TreeHeight, static_cast<uint32_t>(level));
	// 		insertEntry(0, nullptr, n->m_nodeMBR, n->m_identifier);
	// 		// std::cerr << "insertEntry leaf" << std::endl;

	// 		// delete n;
	// 	}

	// 	insertModel(regressionParams);
	// 	// std::cerr << "m_identifiert: " << m_identifier << std::endl;

	// 	m_pTree->writeNode(this);
	// }
	// else
	// {
	// 	// std::cerr << "insertEntry(dataLength, pData, mbr, id) id= " << id << " m_identifier:" << m_identifier << std::endl;

	// 	insertEntry(dataLength, pData, mbr, m_identifier);
	// 	// std::cerr << "directly insert m_identifiert: " << m_identifier << std::endl;
 
	// 	m_pTree->writeNode(this);
	// }
	// return true;
}


void Node::reinsertData(uint32_t dataLength, uint8_t* pData, Region& mbr, id_type id, std::vector<uint32_t>& reinsert, std::vector<uint32_t>& keep)
{
	ReinsertEntry** v = new ReinsertEntry*[m_capacity + 1];

	m_pDataLength[m_children] = dataLength;
	m_pData[m_children] = pData;
	m_ptrMBR[m_children] = m_pTree->m_regionPool.acquire();
	*(m_ptrMBR[m_children]) = mbr;
	m_pIdentifier[m_children] = id;

	PointPtr nc = m_pTree->m_pointPool.acquire();
	m_nodeMBR.getCenter(*nc);
	PointPtr c = m_pTree->m_pointPool.acquire();

	for (uint32_t u32Child = 0; u32Child < m_capacity + 1; ++u32Child)
	{
		try
		{
			v[u32Child] = new ReinsertEntry(u32Child, 0.0);
		}
		catch (...)
		{
			for (uint32_t i = 0; i < u32Child; ++i) delete v[i];
			delete[] v;
			throw;
		}

		m_ptrMBR[u32Child]->getCenter(*c);

		// calculate relative distance of every entry from the node MBR (ignore square root.)
		for (uint32_t cDim = 0; cDim < m_nodeMBR.m_dimension; ++cDim)
		{
			double d = nc->m_pCoords[cDim] - c->m_pCoords[cDim];
			v[u32Child]->m_dist += d * d;
		}
	}

	// sort by increasing order of distances.
	::qsort(v, m_capacity + 1, sizeof(ReinsertEntry*), ReinsertEntry::compareReinsertEntry);

	uint32_t cReinsert = static_cast<uint32_t>(std::floor((m_capacity + 1) * m_pTree->m_reinsertFactor));

	uint32_t cCount;

	for (cCount = 0; cCount < m_capacity + 1; ++cCount)
	{
		if (cCount < m_capacity + 1 - cReinsert)
		{
			// Keep all but cReinsert nodes
			keep.push_back(v[cCount]->m_index);
		}
		else
		{
			// Remove cReinsert nodes which will be
			// reinserted into the tree. Since our array
			// is already sorted in ascending order this
			// matches the order suggested in the paper.
			reinsert.push_back(v[cCount]->m_index);
		}
		delete v[cCount];
	}

	delete[] v;
}

void Node::LearnedIndexSplit(uint32_t dataLength, uint8_t* pData, Region& mbr, id_type id, std::vector<uint32_t>& group1, std::vector<uint32_t>& group2)
{
	
	uint32_t u32Child;
	uint32_t minimumLoad = static_cast<uint32_t>(std::floor(m_capacity * m_pTree->m_fillFactor));

	// use this mask array for marking visited entries.
	uint8_t* mask = new uint8_t[m_capacity + 1];
	memset(mask, 0, m_capacity + 1);

	// insert new data in the node for easier manipulation. Data arrays are always
	// by one larger than node capacity.
	m_pDataLength[m_capacity] = dataLength;
	m_pData[m_capacity] = pData;
	m_ptrMBR[m_capacity] = m_pTree->m_regionPool.acquire();
	*(m_ptrMBR[m_capacity]) = mbr;
	m_pIdentifier[m_capacity] = id;
	// m_totalDataLength does not need to be increased here.

	// initialize each group with the seed entries.
	uint32_t seed1, seed2;
	pickSeeds(seed1, seed2);

	group1.push_back(seed1);
	group2.push_back(seed2);

	mask[seed1] = 1;
	mask[seed2] = 1;

	// find MBR of each group.
	RegionPtr mbr1 = m_pTree->m_regionPool.acquire();
	*mbr1 = *(m_ptrMBR[seed1]);
	RegionPtr mbr2 = m_pTree->m_regionPool.acquire();
	*mbr2 = *(m_ptrMBR[seed2]);

	// count how many entries are left unchecked (exclude the seeds here.)
	uint32_t cRemaining = m_capacity + 1 - 2;

	while (cRemaining > 0)
	{
		if (minimumLoad - group1.size() == cRemaining)
		{
			// all remaining entries must be assigned to group1 to comply with minimun load requirement.
			for (u32Child = 0; u32Child < m_capacity + 1; ++u32Child)
			{
				if (mask[u32Child] == 0)
				{
					group1.push_back(u32Child);
					mask[u32Child] = 1;
					--cRemaining;
				}
			}
		}
		else if (minimumLoad - group2.size() == cRemaining)
		{
			// all remaining entries must be assigned to group2 to comply with minimun load requirement.
			for (u32Child = 0; u32Child < m_capacity + 1; ++u32Child)
			{
				if (mask[u32Child] == 0)
				{
					group2.push_back(u32Child);
					mask[u32Child] = 1;
					--cRemaining;
				}
			}
		}
		else
		{
			// For all remaining entries compute the difference of the cost of grouping an
			// entry in either group. When done, choose the entry that yielded the maximum
			// difference. In case of linear split, select any entry (e.g. the first one.)
			uint32_t sel;
			double md1 = 0.0, md2 = 0.0;
			double m = -std::numeric_limits<double>::max();
			double d1, d2, d;
			double a1 = mbr1->getArea();
			double a2 = mbr2->getArea();

			RegionPtr a = m_pTree->m_regionPool.acquire();
			RegionPtr b = m_pTree->m_regionPool.acquire();

			for (u32Child = 0; u32Child < m_capacity + 1; ++u32Child)
			{
				if (mask[u32Child] == 0)
				{
					mbr1->getCombinedRegion(*a, *(m_ptrMBR[u32Child]));
					d1 = a->getArea() - a1;
					mbr2->getCombinedRegion(*b, *(m_ptrMBR[u32Child]));
					d2 = b->getArea() - a2;
					d = std::abs(d1 - d2);

					if (d > m)
					{
						m = d;
						md1 = d1; md2 = d2;
						sel = u32Child;
						// if (m_pTree->m_treeVariant== RV_LINEAR || m_pTree->m_treeVariant == RV_RSTAR) break;
					}
				}
			}

			// determine the group where we should add the new entry.
			int32_t group = -1;

			if (md1 < md2)
			{
				group1.push_back(sel);
				group = 1;
			}
			else if (md2 < md1)
			{
				group2.push_back(sel);
				group = 2;
			}
			else if (a1 < a2)
			{
				group1.push_back(sel);
				group = 1;
			}
			else if (a2 < a1)
			{
				group2.push_back(sel);
				group = 2;
			}
			else if (group1.size() < group2.size())
			{
				group1.push_back(sel);
				group = 1;
			}
			else if (group2.size() < group1.size())
			{
				group2.push_back(sel);
				group = 2;
			}
			else
			{
				group1.push_back(sel);
				group = 1;
			}
			mask[sel] = 1;
			--cRemaining;
			if (group == 1)
			{
				mbr1->combineRegion(*(m_ptrMBR[sel]));
			}
			else
			{
				mbr2->combineRegion(*(m_ptrMBR[sel]));
			}
		}
	}

	delete[] mask;
}

void Node::rlLearnedIndexSplit(uint32_t dataLength, uint8_t* pData, Region& mbr, id_type id, std::vector<uint32_t>& group1, std::vector<uint32_t>& group2)
{
	uint32_t minimumLoad = static_cast<uint32_t>(std::floor(m_capacity * m_pTree->m_fillFactor));

	RstarSplitEntry** dataLeft = nullptr;
	RstarSplitEntry** dataBottom = nullptr;

	try
	{
		dataLeft = new RstarSplitEntry*[m_capacity + 1];
		dataBottom = new RstarSplitEntry*[m_capacity + 1];
	}
	catch (...)
	{
		delete[] dataLeft;
		throw;
	}

	m_pDataLength[m_capacity] = dataLength;
	m_pData[m_capacity] = pData;
	m_ptrMBR[m_capacity] = m_pTree->m_regionPool.acquire();
	*(m_ptrMBR[m_capacity]) = mbr;
	m_pIdentifier[m_capacity] = id;
	// m_totalDataLength does not need to be increased here.
	uint32_t u32Child = 0, cDim, cIndex;
	for (u32Child = 0; u32Child <= m_capacity; ++u32Child)
	{
		try
		{
			dataLeft[u32Child] = new RstarSplitEntry(m_ptrMBR[u32Child].get(), u32Child, 0);
			dataBottom[u32Child] = new RstarSplitEntry(m_ptrMBR[u32Child].get(), u32Child, 1);
		}
		catch (...)
		{
			for (uint32_t i = 0; i < u32Child; ++i) delete dataLeft[i];
			delete[] dataLeft;
			delete[] dataBottom;
			throw;
		}

	}

	// 2. PrepareSplitLocations
	std::vector<SplitLocation> split_locations((m_capacity - 2 * minimumLoad + 2) * 2);
	::qsort(dataLeft, m_capacity + 1, sizeof(RstarSplitEntry*), RstarSplitEntry::compareLow);
	Region prefix, suffix;
	prefix = *(dataLeft[0]->m_pRegion);
	suffix = *(dataLeft[m_capacity - minimumLoad + 1]->m_pRegion);
	int loc = 0;
	for (u32Child = 1; u32Child < minimumLoad - 1; ++u32Child)
	{
		prefix.combineRegion(*(dataLeft[u32Child]->m_pRegion));
	}

	for (u32Child = m_capacity - minimumLoad + 2; u32Child <= m_capacity; ++u32Child)
	{
		suffix.combineRegion(*(dataLeft[u32Child]->m_pRegion));
	}

	for(u32Child = minimumLoad - 1; u32Child < m_capacity - minimumLoad + 1; ++u32Child)
	{
		prefix.combineRegion(*(dataLeft[u32Child]->m_pRegion));
		Region remaining(suffix);
		for(int j = u32Child + 1; j < m_capacity - minimumLoad + 1; j++){
			remaining.combineRegion(*(dataLeft[j]->m_pRegion));
		}
		split_locations[loc].perimeter1 = std::max(prefix.getMargin(), remaining.getMargin());
		split_locations[loc].perimeter2 = std::min(prefix.getMargin(), remaining.getMargin());
		split_locations[loc].area1 = std::max(prefix.getArea(), remaining.getArea());
		split_locations[loc].area2 = std::min(prefix.getArea(), remaining.getArea());
		split_locations[loc].overlap = prefix.getIntersectingArea(remaining);
		split_locations[loc].location = u32Child;
		split_locations[loc].dimension = 0;
		loc += 1;
	}

	::qsort(dataBottom, m_capacity + 1, sizeof(RstarSplitEntry*), RstarSplitEntry::compareLow);
	prefix = *(dataBottom[0]->m_pRegion);
	suffix = *(dataBottom[m_capacity - minimumLoad + 1]->m_pRegion);


	
	for (u32Child = 1; u32Child < minimumLoad - 1; ++u32Child)
	{
		prefix.combineRegion(*(dataBottom[u32Child]->m_pRegion));
	}

	for (u32Child = m_capacity - minimumLoad + 2; u32Child <= m_capacity; ++u32Child)
	{
		suffix.combineRegion(*(dataBottom[u32Child]->m_pRegion));
	}

	for(u32Child = minimumLoad - 1; u32Child < m_capacity - minimumLoad + 1; ++u32Child)
	{
		prefix.combineRegion(*(dataBottom[u32Child]->m_pRegion));
		Region remaining(suffix);
		for(int j = u32Child + 1; j < m_capacity - minimumLoad + 1; j++){
			remaining.combineRegion(*(dataBottom[j]->m_pRegion));
		}
		split_locations[loc].perimeter1 = std::max(prefix.getMargin(), remaining.getMargin());
		split_locations[loc].perimeter2 = std::min(prefix.getMargin(), remaining.getMargin());
		split_locations[loc].area1 = std::max(prefix.getArea(), remaining.getArea());
		split_locations[loc].area2 = std::min(prefix.getArea(), remaining.getArea());
		split_locations[loc].overlap = prefix.getIntersectingArea(remaining);
		split_locations[loc].location = u32Child;
		split_locations[loc].dimension = 1;
		loc += 1;
	}
	int non_overlap_split_num = 0;
	for (int i = 0; i < split_locations.size(); i++)
	{
		if (split_locations[i].overlap == 0)
		{
			non_overlap_split_num += 1;
		}
	}
	if (non_overlap_split_num <= 1)
	{
		LearnedIndexSplit(dataLength, pData, mbr, id, group1, group2);
	}
	else
	{

		// std::cerr << "------------LearnedIndexSplit------------" << std::endl;
		std::vector<uint32_t> LearnedIndex_group1, LearnedIndex_group2;
		LearnedIndexSplit(dataLength, pData, mbr, id, LearnedIndex_group1, LearnedIndex_group2);
		Region LearnedIndex_mbr1 = *(m_ptrMBR[LearnedIndex_group1[0]]);
		Region LearnedIndex_mbr2 = *(m_ptrMBR[LearnedIndex_group2[0]]);
		for (int i = 1; i < LearnedIndex_group1.size(); i++)
		{
			LearnedIndex_mbr1.combineRegion(*(m_ptrMBR[LearnedIndex_group1[i]]));
		}
		for (int i = 1; i < LearnedIndex_group2.size(); i++)
		{
			LearnedIndex_mbr2.combineRegion(*(m_ptrMBR[LearnedIndex_group2[i]]));
		}
		// std::cerr << "MBR1: " << "m_xmin" << mbr1.m_pLow[0] << " m_xmax" << mbr1.m_pHigh[0] << " m_ymin " << mbr1.m_pLow[1] << " m_ymax" << mbr1.m_pHigh[1] << std::endl;
		// std::cerr << "MBR2: " << "m_xmin" << mbr2.m_pLow[0] << " m_xmax" << mbr2.m_pHigh[0] << " m_ymin " << mbr2.m_pLow[1] << " m_ymax" << mbr2.m_pHigh[1] << std::endl;

		// std::cerr << "MBR Area: " << LearnedIndex_mbr1.getArea() + LearnedIndex_mbr2.getArea() << std::endl;
		// std::cerr << "MBR Margin: " << LearnedIndex_mbr1.getMargin() + LearnedIndex_mbr2.getMargin() << std::endl;
		
		// else
		// {
		std::vector<uint32_t> candidate_split_action(m_pTree->m_rlr_action_space_size);

		std::vector<std::pair<double, int> > zero_ovlp_splits;
		for (int i = 0; i < split_locations.size(); i++)
		{
			if (split_locations[i].overlap == 0)
			{
				// double perimeter = std::max(split_locations[i].perimeter1, split_locations[i].perimeter2);
				double area = std::max(split_locations[i].area1, split_locations[i].area2);
				zero_ovlp_splits.emplace_back(area, i);
			}
		}
		sort(zero_ovlp_splits.begin(), zero_ovlp_splits.end());
		double maxArea = std::numeric_limits<double>::min();
		double minArea = std::numeric_limits<double>::max();
		double maxPerimeter = std::numeric_limits<double>::min();
		double minPerimeter = std::numeric_limits<double>::max();

		// Init states
		std::vector<double> states(m_pTree->m_rlr_action_space_size * m_pTree->m_rlr_scale);

		for (int i = 0; i < m_pTree->m_rlr_action_space_size; i++)
		{
			int idx = zero_ovlp_splits[i].second;
			candidate_split_action[i] = idx;
			states[i * 4 + 0] = split_locations[idx].area1;
			states[i * 4 + 1] = split_locations[idx].area2;
			states[i * 4 + 2] = split_locations[idx].perimeter1;
			states[i * 4 + 3] = split_locations[idx].perimeter2;
			maxArea = std::max(maxArea, states[i * 4 + 0]);
			minArea = std::min(minArea, states[i * 4 + 1]);
			maxPerimeter = std::max(maxPerimeter, states[i * 4 + 2]);
			minPerimeter = std::min(minPerimeter, states[i * 4 + 3]);
		}
		for (int i = 0; i < 2; i++)
		{
			states[i * 4 + 0] = (states[i * 4 + 0] - minArea) / (maxArea - minArea + 0.1);
			states[i * 4 + 1] = (states[i * 4 + 1] - minArea) / (maxArea - minArea + 0.1);
			states[i * 4 + 2] = (states[i * 4 + 2] - minPerimeter) / (maxPerimeter - minPerimeter + 0.1);
			states[i * 4 + 3] = (states[i * 4 + 3] - minPerimeter) / (maxPerimeter - minPerimeter + 0.1);
		}

		// Invoke the libtorch model
		uint32_t best = m_pTree->splitModelForward(states); 
		best = std::min(m_pTree->m_rlr_action_space_size - 1, best);
		// Do split
		int dim = split_locations[candidate_split_action[best]].dimension;
		int split_loc = split_locations[candidate_split_action[best]].location;
		// std::cerr << "------------RLLearnedIndexSplit------------" << std::endl;
		// std::cerr << "split best: " << best << " dim: " << dim << " split_loc: " << split_loc << std::endl;

		if (dim == 0) 
		{
			for (u32Child = 0; u32Child <= split_loc; ++u32Child)
			{
				group1.push_back(dataLeft[u32Child]->m_index);
			}
			for (u32Child = split_loc + 1; u32Child <= m_capacity; ++u32Child)
			{
				group2.push_back(dataLeft[u32Child]->m_index);
			}
		}
		else
		{
			for (u32Child = 0; u32Child <= split_loc; ++u32Child)
			{
				group1.push_back(dataBottom[u32Child]->m_index);
			}
			for (u32Child = split_loc + 1; u32Child <= m_capacity; ++u32Child)
			{
				group2.push_back(dataBottom[u32Child]->m_index);
			}
		}

		Region mbr1 = *(m_ptrMBR[group1[0]]);
		Region mbr2 = *(m_ptrMBR[group2[0]]);
		for (int i = 1; i < group1.size(); i++)
		{
			mbr1.combineRegion(*(m_ptrMBR[group1[i]]));
		}
		for (int i = 1; i < group2.size(); i++)
		{
			mbr2.combineRegion(*(m_ptrMBR[group2[i]]));
		}

		// std::cerr << "MBR Area: " << mbr1.getArea() + mbr2.getArea() << std::endl;
		// std::cerr << "MBR Margin: " << mbr1.getMargin() + mbr2.getMargin() << std::endl;
		if (LearnedIndex_mbr1.getMargin() + LearnedIndex_mbr2.getMargin() < mbr1.getMargin() + mbr2.getMargin()) {
			group1 = LearnedIndex_group1;
			group2 = LearnedIndex_group2;
		}

		// std::cerr << "MBR1: " << "m_xmin" << mbr1.m_pLow[0] << " m_xmax" << mbr1.m_pHigh[0] << " m_ymin " << mbr1.m_pLow[1] << " m_ymax" << mbr1.m_pHigh[1] << std::endl;
		// std::cerr << "MBR2: " << "m_xmin" << mbr2.m_pLow[0] << " m_xmax" << mbr2.m_pHigh[0] << " m_ymin " << mbr2.m_pLow[1] << " m_ymax" << mbr2.m_pHigh[1] << std::endl;

	}

}

void Node::rstarSplit(uint32_t dataLength, uint8_t* pData, Region& mbr, id_type id, std::vector<uint32_t>& group1, std::vector<uint32_t>& group2)
{
	RstarSplitEntry** dataLow = nullptr;
	RstarSplitEntry** dataHigh = nullptr;

	try
	{
		dataLow = new RstarSplitEntry*[m_capacity + 1];
		dataHigh = new RstarSplitEntry*[m_capacity + 1];
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
	// m_totalDataLength does not need to be increased here.

	uint32_t nodeSPF = static_cast<uint32_t>(
		std::floor((m_capacity + 1) * m_pTree->m_splitDistributionFactor));
	uint32_t splitDistribution = (m_capacity + 1) - (2 * nodeSPF) + 2;

	uint32_t u32Child = 0, cDim, cIndex;

	for (u32Child = 0; u32Child <= m_capacity; ++u32Child)
	{
		try
		{
			dataLow[u32Child] = new RstarSplitEntry(m_ptrMBR[u32Child].get(), u32Child, 0);
		}
		catch (...)
		{
			for (uint32_t i = 0; i < u32Child; ++i) delete dataLow[i];
			delete[] dataLow;
			delete[] dataHigh;
			throw;
		}

		dataHigh[u32Child] = dataLow[u32Child];
	}

	double minimumMargin = std::numeric_limits<double>::max();
	uint32_t splitAxis = std::numeric_limits<uint32_t>::max();
	uint32_t sortOrder = std::numeric_limits<uint32_t>::max();

	// chooseSplitAxis.
	for (cDim = 0; cDim < m_pTree->m_dimension; ++cDim)
	{
		::qsort(dataLow, m_capacity + 1, sizeof(RstarSplitEntry*), RstarSplitEntry::compareLow);
		::qsort(dataHigh, m_capacity + 1, sizeof(RstarSplitEntry*), RstarSplitEntry::compareHigh);

		// calculate sum of margins and overlap for all distributions.
		double marginl = 0.0;
		double marginh = 0.0;

		Region bbl1, bbl2, bbh1, bbh2;

		for (u32Child = 1; u32Child <= splitDistribution; ++u32Child)
		{
			uint32_t l = nodeSPF - 1 + u32Child;

			bbl1 = *(dataLow[0]->m_pRegion);
			bbh1 = *(dataHigh[0]->m_pRegion);

			for (cIndex = 1; cIndex < l; ++cIndex)
			{
				bbl1.combineRegion(*(dataLow[cIndex]->m_pRegion));
				bbh1.combineRegion(*(dataHigh[cIndex]->m_pRegion));
			}

			bbl2 = *(dataLow[l]->m_pRegion);
			bbh2 = *(dataHigh[l]->m_pRegion);

			for (cIndex = l + 1; cIndex <= m_capacity; ++cIndex)
			{
				bbl2.combineRegion(*(dataLow[cIndex]->m_pRegion));
				bbh2.combineRegion(*(dataHigh[cIndex]->m_pRegion));
			}

			marginl += bbl1.getMargin() + bbl2.getMargin();
			marginh += bbh1.getMargin() + bbh2.getMargin();
		} // for (u32Child)

		double margin = std::min(marginl, marginh);

		// keep minimum margin as split axis.
		if (margin < minimumMargin)
		{
			minimumMargin = margin;
			splitAxis = cDim;
			sortOrder = (marginl < marginh) ? 0 : 1;
		}

		// increase the dimension according to which the data entries should be sorted.
		for (u32Child = 0; u32Child <= m_capacity; ++u32Child)
		{
			dataLow[u32Child]->m_sortDim = cDim + 1;
		}
	} // for (cDim)

	for (u32Child = 0; u32Child <= m_capacity; ++u32Child)
	{
		dataLow[u32Child]->m_sortDim = splitAxis;
	}

	::qsort(dataLow, m_capacity + 1, sizeof(RstarSplitEntry*), (sortOrder == 0) ? RstarSplitEntry::compareLow : RstarSplitEntry::compareHigh);

	double ma = std::numeric_limits<double>::max();
	double mo = std::numeric_limits<double>::max();
	uint32_t splitPoint = std::numeric_limits<uint32_t>::max();

	Region bb1, bb2;

	for (u32Child = 1; u32Child <= splitDistribution; ++u32Child)
	{
		uint32_t l = nodeSPF - 1 + u32Child;

		bb1 = *(dataLow[0]->m_pRegion);

		for (cIndex = 1; cIndex < l; ++cIndex)
		{
			bb1.combineRegion(*(dataLow[cIndex]->m_pRegion));
		}

		bb2 = *(dataLow[l]->m_pRegion);

		for (cIndex = l + 1; cIndex <= m_capacity; ++cIndex)
		{
			bb2.combineRegion(*(dataLow[cIndex]->m_pRegion));
		}

		double o = bb1.getIntersectingArea(bb2);

		if (o < mo)
		{
			splitPoint = u32Child;
			mo = o;
			ma = bb1.getArea() + bb2.getArea();
		}
		else if (o == mo)
		{
			double a = bb1.getArea() + bb2.getArea();

			if (a < ma)
			{
				splitPoint = u32Child;
				ma = a;
			}
		}
	} // for (u32Child)

	uint32_t l1 = nodeSPF - 1 + splitPoint;

	for (cIndex = 0; cIndex < l1; ++cIndex)
	{
		group1.push_back(dataLow[cIndex]->m_index);
		delete dataLow[cIndex];
	}

	for (cIndex = l1; cIndex <= m_capacity; ++cIndex)
	{
		group2.push_back(dataLow[cIndex]->m_index);
		delete dataLow[cIndex];
	}

	delete[] dataLow;
	delete[] dataHigh;
}

void Node::pickSeeds(uint32_t& index1, uint32_t& index2)
{
	double separation = -std::numeric_limits<double>::max();
	double inefficiency = -std::numeric_limits<double>::max();
	uint32_t cDim, u32Child, cIndex;

	// switch (m_pTree->m_treeVariant)
	// {
	// 	case RV_LINEAR:
	// 	case RV_RLLearnedIndex:
	// 	case RV_RSTAR:
	// 		for (cDim = 0; cDim < m_pTree->m_dimension; ++cDim)
	// 		{
	// 			double leastLower = m_ptrMBR[0]->m_pLow[cDim];
	// 			double greatestUpper = m_ptrMBR[0]->m_pHigh[cDim];
	// 			uint32_t greatestLower = 0;
	// 			uint32_t leastUpper = 0;
	// 			double width;

	// 			for (u32Child = 1; u32Child <= m_capacity; ++u32Child)
	// 			{
	// 				if (m_ptrMBR[u32Child]->m_pLow[cDim] > m_ptrMBR[greatestLower]->m_pLow[cDim]) greatestLower = u32Child;
	// 				if (m_ptrMBR[u32Child]->m_pHigh[cDim] < m_ptrMBR[leastUpper]->m_pHigh[cDim]) leastUpper = u32Child;

	// 				leastLower = std::min(m_ptrMBR[u32Child]->m_pLow[cDim], leastLower);
	// 				greatestUpper = std::max(m_ptrMBR[u32Child]->m_pHigh[cDim], greatestUpper);
	// 			}

	// 			width = greatestUpper - leastLower;
	// 			if (width <= 0) width = 1;

	// 			double f = (m_ptrMBR[greatestLower]->m_pLow[cDim] - m_ptrMBR[leastUpper]->m_pHigh[cDim]) / width;

	// 			if (f > separation)
	// 			{
	// 				index1 = leastUpper;
	// 				index2 = greatestLower;
	// 				separation = f;
	// 			}
	// 		}  // for (cDim)

	// 		if (index1 == index2)
	// 		{
	// 			if (index2 == 0) ++index2;
	// 			else --index2;
	// 		}

	// 		break;
	// 	case RV_QUADRATIC:
	// 		// for each pair of Regions (account for overflow Region too!)
	// 		for (u32Child = 0; u32Child < m_capacity; ++u32Child)
	// 		{
	// 			double a = m_ptrMBR[u32Child]->getArea();

	// 			for (cIndex = u32Child + 1; cIndex <= m_capacity; ++cIndex)
	// 			{
	// 				// get the combined MBR of those two entries.
	// 				Region r;
	// 				m_ptrMBR[u32Child]->getCombinedRegion(r, *(m_ptrMBR[cIndex]));

	// 				// find the inefficiency of grouping these entries together.
	// 				double d = r.getArea() - a - m_ptrMBR[cIndex]->getArea();

	// 				if (d > inefficiency)
	// 				{
	// 					inefficiency = d;
	// 					index1 = u32Child;
	// 					index2 = cIndex;
	// 				}
	// 			}  // for (cIndex)
	// 		} // for (u32Child)

	// 		break;
	// 	default:
	// 		throw Tools::NotSupportedException("Node::pickSeeds: Tree variant not supported.");
	// }
}

void Node::condenseTree(std::stack<NodePtr>& toReinsert, std::stack<id_type>& pathBuffer, NodePtr& ptrThis)
{
	uint32_t minimumLoad = static_cast<uint32_t>(std::floor(m_capacity * m_pTree->m_fillFactor));

	if (pathBuffer.empty())
	{
		// eliminate root if it has only one child.
		if (m_level != 0 && m_children == 1)
		{
			NodePtr ptrN = m_pTree->readNode(m_pIdentifier[0]);
			m_pTree->deleteNode(ptrN.get());
			ptrN->m_identifier = m_pTree->m_rootID;
			m_pTree->writeNode(ptrN.get());

			m_pTree->m_stats.m_nodesInLevel.pop_back();
			m_pTree->m_stats.m_u32TreeHeight -= 1;
			// HACK: pending deleteNode for deleted child will decrease nodesInLevel, later on.
			m_pTree->m_stats.m_nodesInLevel[m_pTree->m_stats.m_u32TreeHeight - 1] = 2;
		}
		else
		{
			// due to data removal.
			if (m_pTree->m_bTightMBRs)
			{
				for (uint32_t cDim = 0; cDim < m_nodeMBR.m_dimension; ++cDim)
				{
					m_nodeMBR.m_pLow[cDim] = std::numeric_limits<double>::max();
					m_nodeMBR.m_pHigh[cDim] = -std::numeric_limits<double>::max();

					for (uint32_t u32Child = 0; u32Child < m_children; ++u32Child)
					{
						m_nodeMBR.m_pLow[cDim] = std::min(m_nodeMBR.m_pLow[cDim], m_ptrMBR[u32Child]->m_pLow[cDim]);
						m_nodeMBR.m_pHigh[cDim] = std::max(m_nodeMBR.m_pHigh[cDim], m_ptrMBR[u32Child]->m_pHigh[cDim]);
					}
				}
			}

            // write parent node back to storage.
			m_pTree->writeNode(this);
		}
	}
	else
	{
		id_type cParent = pathBuffer.top(); pathBuffer.pop();
		NodePtr ptrParent = m_pTree->readNode(cParent);
		Index* p = static_cast<Index*>(ptrParent.get());

		// find the entry in the parent, that points to this node.
		uint32_t child;

		for (child = 0; child != p->m_children; ++child)
		{
			if (p->m_pIdentifier[child] == m_identifier) break;
		}

		if (m_children < minimumLoad)
		{
			// used space less than the minimum
			// 1. eliminate node entry from the parent. deleteEntry will fix the parent's MBR.
			p->deleteEntry(child);
			// 2. add this node to the stack in order to reinsert its entries.
			toReinsert.push(ptrThis);
		}
		else
		{
			// adjust the entry in 'p' to contain the new bounding region of this node.
			*(p->m_ptrMBR[child]) = m_nodeMBR;

			// global recalculation necessary since the MBR can only shrink in size,
			// due to data removal.
			if (m_pTree->m_bTightMBRs)
			{
				for (uint32_t cDim = 0; cDim < p->m_nodeMBR.m_dimension; ++cDim)
				{
					p->m_nodeMBR.m_pLow[cDim] = std::numeric_limits<double>::max();
					p->m_nodeMBR.m_pHigh[cDim] = -std::numeric_limits<double>::max();

					for (uint32_t u32Child = 0; u32Child < p->m_children; ++u32Child)
					{
						p->m_nodeMBR.m_pLow[cDim] = std::min(p->m_nodeMBR.m_pLow[cDim], p->m_ptrMBR[u32Child]->m_pLow[cDim]);
						p->m_nodeMBR.m_pHigh[cDim] = std::max(p->m_nodeMBR.m_pHigh[cDim], p->m_ptrMBR[u32Child]->m_pHigh[cDim]);
					}
				}
			}
		}

		// write parent node back to storage.
		m_pTree->writeNode(p);

		p->condenseTree(toReinsert, pathBuffer, ptrParent);
	}
}
