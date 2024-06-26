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

#pragma once

namespace SpatialIndex
{
	namespace KDTree
	{
		enum KDTreeVariant
		{
			KD_NORMAL = 0x0,
			KD_GREEDY = 0x1,
			QD_NORMAL = 0x2
		};

		enum LoadMethod
		{
			LOAD_KD = 0x0,
			LOAD_KD_GREEDY = 0x1,
			LOAD_QD = 0x2
		};

		enum PersistenObjectIdentifier
		{
			PersistentIndex = 0x1,
			PersistentLeaf = 0x2
		};

		enum RangeQueryType
		{
			ContainmentQuery = 0x1,
			IntersectionQuery = 0x2
		};

		class SIDX_DLL Data : public IData, public Tools::ISerializable
		{
		public:
			Data(uint32_t len, uint8_t* pData, Region& r, id_type id);
			~Data() override;

			Data* clone() override;
			id_type getIdentifier() const override;
			void getShape(IShape** out) const override;
			void getData(uint32_t& len, uint8_t** data) const override;
			uint32_t getByteArraySize() override;
			void loadFromByteArray(const uint8_t* data) override;
			void storeToByteArray(uint8_t** data, uint32_t& len) override;

			id_type m_id;
			Region m_region;
			uint8_t* m_pData;
			uint32_t m_dataLength;
		}; // Data

		SIDX_DLL ISpatialIndex* returnKDTree(IStorageManager& ind, Tools::PropertySet& in);
		SIDX_DLL ISpatialIndex* createNewKDTree(
			IStorageManager& sm,
			double fillFactor,
			uint32_t indexCapacity,
			uint32_t leafCapacity,
			uint32_t dimension,
			KDTreeVariant rv,
			id_type& indexIdentifier
		);
		// SIDX_DLL ISpatialIndex* createNewQDTree(
		// 	IStorageManager& sm,
		// 	double fillFactor,
		// 	uint32_t indexCapacity,
		// 	uint32_t leafCapacity,
		// 	uint32_t dimension,
		// 	KDTreeVariant rv,
		// 	id_type& indexIdentifier,
		// 	const std::string& modelPath
		// );
		SIDX_DLL ISpatialIndex* createNewKDTree(
			LoadMethod m,
			IDataStream& stream,
			IStorageManager& sm,
			double fillFactor,
			uint32_t indexCapacity,
			uint32_t leafCapacity,
			uint32_t dimension,
			KDTreeVariant rv,
			id_type& indexIdentifier
		);
		SIDX_DLL ISpatialIndex* createNewKDTree(
			LoadMethod m,
			IDataStream& stream,
			IStorageManager& sm,
			double fillFactor,
			uint32_t indexCapacity,
			uint32_t leafCapacity,
			uint32_t dimension,
			KDTreeVariant rv,
			id_type& indexIdentifier,
			IDataStream& queryStream
		);
		SIDX_DLL ISpatialIndex* createNewQDTree(
			LoadMethod m,
			IDataStream& stream,
			IStorageManager& sm,
			double fillFactor,
			uint32_t indexCapacity,
			uint32_t leafCapacity,
			uint32_t dimension,
			KDTreeVariant rv,
			id_type& indexIdentifier,
			IDataStream& queryStream,
			const std::string& modelPath,
			int action_space_sizes
		);
		SIDX_DLL ISpatialIndex* createNewKDTree(
			LoadMethod m,
			IDataStream& stream,
			IStorageManager& sm,
			Tools::PropertySet& ps,
			id_type& indexIdentifier
		);
		SIDX_DLL ISpatialIndex* loadKDTree(IStorageManager& in, id_type indexIdentifier);
	}
}
