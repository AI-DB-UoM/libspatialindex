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

#include <memory>

namespace SpatialIndex
{
	namespace KDTree
	{
		class ExternalSorter
		{
		public:
			class Record
			{
			public:
				Record();
				Record(const Region& r, id_type id, uint32_t len, uint8_t* pData, uint32_t s);
				~Record();

				bool operator<(const Record& r) const;

				void storeToFile(Tools::TemporaryFile& f);
				void loadFromFile(Tools::TemporaryFile& f);

				struct SortAscending
				{
					bool operator()(Record* const r1, Record* const r2)
					{
						if (*r1 < *r2) return true;
						else return false;
					}
				};

			public:
				Region m_r;
				id_type m_id;
				uint32_t m_len;
				uint8_t* m_pData{nullptr};
				uint32_t m_s; // sort index
			};

		public:
			ExternalSorter(uint32_t u32PageSize, uint32_t u32BufferPages);
			virtual ~ExternalSorter();

			void insert(Record* r);
			void sort();
			void finishLoad();
			Record* getNextRecord();
			uint64_t getTotalEntries() const;

		private:
			class PQEntry
			{
			public:
				PQEntry(Record* r, uint32_t u32Index) : m_r(r), m_u32Index(u32Index) {}

				struct SortAscending
				{
					bool operator()(const PQEntry& e1, const PQEntry& e2)
					{
						if (*(e1.m_r) < *(e2.m_r)) return true;
						else return false;
					}
				};

				Record* m_r;
				uint32_t m_u32Index;
			};

		private:
			bool m_bInsertionPhase;
			uint32_t m_u32PageSize;
			uint32_t m_u32BufferPages;
            std::shared_ptr<Tools::TemporaryFile> m_sortedFile;
			std::list<std::shared_ptr<Tools::TemporaryFile> > m_runs;
			std::vector<Record*> m_buffer;
			uint64_t m_u64TotalEntries;
			uint32_t m_stI;
		};

		class BulkLoader
		{
		public:
			void topDownPartitioning(
				KDTree* pTree,
				IDataStream& stream,
				uint32_t bindex,
				uint32_t bleaf,
				uint32_t pageSize, // The number of node entries per page.
				uint32_t numberOfPages // The total number of pages to use.
			);

			void topDownGreedyPartitioning(
				KDTree* pTree,
				IDataStream& stream,
				IDataStream& queryStream,
				uint32_t bindex,
				uint32_t bleaf,
				uint32_t pageSize, // The number of node entries per page.
				uint32_t numberOfPages // The total number of pages to use.
			);

			std::vector<std::pair<uint32_t, double>> generageCandidateCutPos(
				std::vector<Region>& regions
			);

			uint64_t calculateSkip(
				std::vector<std::pair<Region, uint64_t>>& candidateSplits, 
				std::vector<Region>& regions,
				uint32_t dimension,
				double value
			);

		protected:
			void partition(
				KDTree* pTree,
				std::shared_ptr<ExternalSorter> es,
				uint32_t dimension,
				uint32_t indexSize,
				uint32_t leafSize,
				uint32_t level,
				// uint32_t depth,
				std::shared_ptr<ExternalSorter> es2,
				uint32_t pageSize,
				uint32_t numberOfPages
			);

			bool greedyPartition(
				SpatialIndex::KDTree::KDTree* pTree,
				std::vector<ExternalSorter::Record *> tupleSet,
				uint32_t bleaf,
				uint32_t bindex,
				uint32_t level,
				Region parentMBR,
				std::vector<ExternalSorter::Record *>& tupleSet2,
				std::vector<Region>& queryRegions,
				std::vector<std::pair<uint32_t, double>>& candidateCutPos
			);

			Node* createNode(
				KDTree* pTree,
				std::vector<ExternalSorter::Record*>& e,
				uint32_t level
			);
		};
	}
}
