#include <spatialindex/SpatialIndex.h>

using namespace SpatialIndex;

#define INSERT 1
#define DELETE 0
#define QUERY 2

class MyDataStream : public IDataStream
{
public:
	MyDataStream(std::string inputFile) : m_pNext(nullptr)
	{
		m_fin.open(inputFile.c_str());

		if (! m_fin)
			throw Tools::IllegalArgumentException("Input file not found.");

		readNextEntry();
	}

	~MyDataStream() override
	{
		if (m_pNext != nullptr) delete m_pNext;
	}

	IData* getNext() override
	{
		if (m_pNext == nullptr) return nullptr;

		RTree::Data* ret = m_pNext;
		m_pNext = nullptr;
		readNextEntry();
		return ret;
	}

	bool hasNext() override
	{
		return (m_pNext != nullptr);
	}

	uint32_t size() override
	{
		throw Tools::NotSupportedException("Operation not supported.");
	}

	void rewind() override
	{
		if (m_pNext != nullptr)
		{
			delete m_pNext;
			m_pNext = nullptr;
		}

		m_fin.seekg(0, std::ios::beg);
		readNextEntry();
	}

	void readNextEntry()
	{
		id_type id;
		uint32_t op;
		double low[3], high[3];

		// m_fin >> op >> id >> low[0] >> low[1] >> high[0] >> high[1];
		m_fin >> op >> id >> low[0] >> low[1] >> low[2] >> high[0] >> high[1] >> high[2];

		if (m_fin.good())
		{
			if (op != INSERT)
				throw Tools::IllegalArgumentException(
					"The data input should contain insertions only."
				);

			Region r(low, high, 3);
			m_pNext = new RTree::Data(sizeof(double), reinterpret_cast<uint8_t*>(low), r, id);
				// Associate a bogus data array with every entry for testing purposes.
				// Once the data array is given to RTRee:Data a local copy will be created.
				// Hence, the input data array can be deleted after this operation if not
				// needed anymore.
		}
	}

	std::ifstream m_fin;
	RTree::Data* m_pNext;
};

int main(int argc, char** argv)
{
	try
	{
		if (argc != 7)
		{
			std::cerr << "Usage: " << argv[0] << " input_file tree_file capacity utilization." << std::endl;
			return -1;
		}

		std::string baseName = argv[2];
		double utilization = atof(argv[4]);

		IStorageManager* diskfile = StorageManager::createNewDiskStorageManager(baseName, 4096);
			// Create a new storage manager with the provided base name and a 4K page size.

		StorageManager::IBuffer* file = StorageManager::createNewRandomEvictionsBuffer(*diskfile, 0, true);
			// applies a main memory random buffer on top of the persistent storage manager
			// (LRU buffer, etc can be created the same way).

		MyDataStream stream(argv[1]);

		// Create and bulk load a new RTree with dimensionality 2, using "file" as
		// the StorageManager and the RSTAR splitting policy.
		id_type indexIdentifier;
		ISpatialIndex* tree;
		std::string packMethods = argv[5];
		if (!packMethods.compare("tgs")) {
			std::cout << "bulkloading using tgs...\n";
			tree = RTree::createAndBulkLoadNewRTreeLearned(
				RTree::BLM_TGS, stream, *file, utilization, atoi(argv[3]), atoi(argv[3]), 3, SpatialIndex::RTree::RV_RSTAR, indexIdentifier, argv[6]);
		}
		else if (!packMethods.compare("str")){
			std::cout << "bulkloading using str...\n";
			tree = RTree::createAndBulkLoadNewRTree(
			RTree::BLM_STR, stream, *file, utilization, atoi(argv[3]), atoi(argv[3]), 3, SpatialIndex::RTree::RV_RSTAR, indexIdentifier);
		}

		std::cerr << *tree;
		std::cerr << "Buffer hits: " << file->getHits() << std::endl;
		std::cerr << "Index ID: " << indexIdentifier << std::endl;

		bool ret = tree->isIndexValid();
		if (ret == false) std::cerr << "ERROR: Structure is invalid!" << std::endl;
		else std::cerr << "The stucture seems O.K." << std::endl;

		delete tree;
		delete file;
		delete diskfile;
			// delete the buffer first, then the storage manager
			// (otherwise the the buffer will fail trying to write the dirty entries).
	}
	catch (Tools::Exception& e)
	{
		std::cerr << "******ERROR******" << std::endl;
		std::string s = e.what();
		std::cerr << s << std::endl;
		return -1;
	}

	return 0;
}
