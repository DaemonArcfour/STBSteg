#define _CRT_SECURE_NO_WARNINGS
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <bitset>
#include <vector>
#include <map>
#include <stb_image.h>
#include <stb_image_write.h>
#include <math.h>
#include <shlwapi.h>
#include <algorithm>
#include <chrono>
#include <omp.h>
#pragma comment(lib, "Shlwapi.lib");

#define STEGPLUS_SIGNATURE_RESERVE 9
#define M_PI           3.14159265358979323846

const int checkblock = 0;

enum eSteganographyModes
{
	STEG_LSB,
	STEG_DCT
};

enum eStreamChannels
{
	CHANNEL_RED = 0x01,
	CHANNEL_GREEN = 0x02,
	CHANNEL_BLUE = 0x04,
	CHANNEL_ALPHA = 0x08
};

class CDiscreteCosineTransform
{
private:
	inline double IDCTCalcSingle(double i, double j, double** pDCTBlock, double iBlockSize)
	{
		double IDCTSingle = 0;
		double fVar = (1.0f / sqrt(2.0f * iBlockSize));
		for (double x = 0; x < iBlockSize; x++)
		{
			for (double y = 0; y < iBlockSize; y++)
			{
				double Ci = 1.0;
				double Cj = 1.0;
				if (x == 0)
					Ci = 1.0f / sqrt(2.0f);
				if (y == 0)
					Cj = 1.0f / sqrt(2.0f);

				IDCTSingle += Ci * Cj * pDCTBlock[(int)x][(int)y]
					* cos(((2.0f * i + 1.0f) * x * M_PI) / (2.0f * iBlockSize))
					* cos(((2.0f * j + 1.0f) * y * M_PI) / (2.0f * iBlockSize));
			}
		}
		IDCTSingle = fVar * IDCTSingle;
		return IDCTSingle;
	}

	inline double DCTCalcSingle(double i, double j, double*** pBlock, double iBlockSize)
	{
		double Ci = 1;
		double Cj = 1;
		if (i == 0)
			Ci = 1.0f / sqrt(2.0f);
		if (j == 0)
			Cj = 1.0f / sqrt(2.0f);

		double DCTSingle = 0;
		double fVar = (1.0f / sqrt(2.0f * iBlockSize)) * Ci * Cj;
		for (double x = 0; x < iBlockSize; x++)
		{
			for (double y = 0; y < iBlockSize; y++)
			{
				DCTSingle += *pBlock[(int)x][(int)y] * cos(((2.0f * x + 1.0f) * i * M_PI) / (2.0f * iBlockSize))
					* cos(((2.0f * y + 1.0f) * j * M_PI) / (2.0f * iBlockSize));
			}
		}
		DCTSingle *= fVar;
		return DCTSingle;
	}

	inline double DCTCalcSingle(double i, double j, double** pBlock, double iBlockSize)
	{
		double Ci = 1;
		double Cj = 1;
		if (i == 0)
			Ci = 1.0f / sqrt(2.0f);
		if (j == 0)
			Cj = 1.0f / sqrt(2.0f);

		double DCTSingle = 0;
		double fVar = (1.0f / sqrt(2.0f * iBlockSize)) * Ci * Cj;
		for (double x = 0; x < iBlockSize; x++)
		{
			for (double y = 0; y < iBlockSize; y++)
			{
				DCTSingle += pBlock[(int)x][(int)y] * cos(((2.0f * x + 1.0f) * i * M_PI) / (2.0f * iBlockSize))
					* cos(((2.0f * y + 1.0f) * j * M_PI) / (2.0f * iBlockSize));
			}
		}
		DCTSingle *= fVar;
		return DCTSingle;
	}

public:
	void DCT(double*** pBlock, double** pDCTBlock, int iBlockSize)
	{
		#pragma omp parallel for
		for (int i = 0; i < iBlockSize; i++)
		{
				for (int j = 0; j < iBlockSize; j++)
				{
					pDCTBlock[(int)i][(int)j] = DCTCalcSingle(i, j, pBlock, iBlockSize);
				}
		}
	}

	void DCT(double** pBlock, double** pDCTBlock, int iBlockSize)
	{
		#pragma omp parallel for
		for (int i = 0; i < iBlockSize; i++)
		{
				for (int j = 0; j < iBlockSize; j++)
				{
					pDCTBlock[(int)i][(int)j] = DCTCalcSingle(i, j, pBlock, iBlockSize);
				}
		}
	}

	void IDCT(double** pBlock, double** pDCTBlock, int iBlockSize)
	{
		#pragma omp parallel for
		for (int i = 0; i < iBlockSize; i++)
		{
				for (int j = 0; j < iBlockSize; j++)
				{
					pBlock[(int)i][(int)j] = IDCTCalcSingle(i, j, pDCTBlock, iBlockSize);
				}
		}
	}
};

template<typename T>
class CMatrix2x2
{
private:
	T** m_pMatrix;
	int32_t m_iMatrixWidth;
	int32_t m_iMatrixHeight;

	T** TGenerateMatrix(int32_t iWidth, int32_t iHeight)
	{
		T** Matrix = new T * [iWidth];
		for (int i = 0; i < iWidth; ++i)
		{
			Matrix[i] = new T[iHeight];
		}

		return Matrix;
	}

	T*** TGeneratePointerMatrix(int32_t iWidth, int32_t iHeight)
	{
		T*** pMatrix = new T * *[iWidth];
		for (int i = 0; i < iWidth; ++i)
		{
			pMatrix[i] = new T * [iHeight];
		}

		for (int i = 0; i < iWidth; i++)
		{
			for (int j = 0; j < iHeight; j++)
			{
				pMatrix[i][j] = nullptr;
			}
		}

		return pMatrix;
	}

	T*** TGeneratePointerMesh(T** pMatrix, int32_t x, int32_t y, int32_t iWidth, int32_t iHeight)
	{
		T*** ppMatrix = TGeneratePointerMatrix(iWidth, iHeight);
		for (int i = 0; i < iWidth; i++)
		{
			for (int j = 0; j < iHeight; j++)
			{
				ppMatrix[i][j] = &pMatrix[x + i][y + j];

				if (*ppMatrix[i][j] > INT_MAX || *ppMatrix[i][j] < -INT_MAX)
					printf("%d %d\n", x + i, y + j);
			}
		}

		return ppMatrix;
	}
public:
	CMatrix2x2()
	{
		m_iMatrixWidth = 0;
		m_iMatrixHeight = 0;
		m_pMatrix = nullptr;
	}

	CMatrix2x2(int32_t iWidth, int32_t iHeight)
	{
		CMatrix2x2();
		GenerateMatrix(iWidth, iHeight);
	}

	void GenerateMatrix(int32_t iWidth, int32_t iHeight)
	{
		ClearMatrix();
		m_pMatrix = TGenerateMatrix(iWidth, iHeight);
		m_iMatrixWidth = iWidth;
		m_iMatrixHeight = iHeight;
	}

	T*** GetMatrixSegment(int32_t x, int32_t y, int32_t iWidth, int32_t iHeight)
	{
		if (m_pMatrix == nullptr)
			return nullptr;

		if (iWidth - x > m_iMatrixWidth || iHeight - y > m_iMatrixHeight)
			return nullptr;

		return TGeneratePointerMesh(m_pMatrix, x, y, iWidth, iHeight);
	}

	T** GetMatrix()
	{
		return m_pMatrix;
	}

	int32_t GetMatrixWidth()
	{
		return m_iMatrixWidth;
	}

	int32_t GetMatrixHeight()
	{
		return m_iMatrixHeight;
	}

	void ClearMatrix()
	{
		if (m_pMatrix != nullptr)
		{
			// TODO: 1. Memory leak removal.
			//for (int i = 0; i < m_iMatrixWidth; i++)
			//	delete[] m_pMatrix[i];
			m_pMatrix = nullptr;
			m_iMatrixWidth = 0;
			m_iMatrixHeight = 0;
		}
	}
};

class CChannelBlock
{
private:
	int32_t m_iDimensions;
	CDiscreteCosineTransform m_cDCTInterface;
	double*** m_pBlock;
	CMatrix2x2<double> m_pBlockDCT;
	CMatrix2x2<double> m_pBlockIDCT;

public:
	CChannelBlock()
	{
		m_iDimensions = 0;
		m_pBlock = nullptr;
	}

	CChannelBlock(double*** pBlock, int32_t iDimensions)
	{
		m_pBlock = pBlock;
		m_iDimensions = iDimensions;
		CalculateDCT();
	}

	inline void CalculateDCT()
	{
		m_pBlockDCT.GenerateMatrix(m_iDimensions, m_iDimensions);
		m_cDCTInterface.DCT(m_pBlock, m_pBlockDCT.GetMatrix(), 8);
	}

	CMatrix2x2<double> GetBlockDCT()
	{
		return m_pBlockDCT;
	}

	void UpdateDCT()
	{
		m_pBlockIDCT.GenerateMatrix(m_iDimensions, m_iDimensions);
		m_cDCTInterface.IDCT(m_pBlockIDCT.GetMatrix(), m_pBlockDCT.GetMatrix(), 8);
		
		for (int i = 0; i < m_iDimensions; i++)
		{
			for (int j = 0; j < m_iDimensions; j++)
			{
				*m_pBlock[i][j] = stbi__clamp(m_pBlockIDCT.GetMatrix()[i][j]);
			}
		}
		
	}
};

class CImage
{
private:
	uint8_t* m_uImageData;
	std::map<eStreamChannels, CMatrix2x2<double>> m_vImageChannelMatrix;
	std::map<eStreamChannels, std::vector<CChannelBlock>> m_vChannelBlocks;
	int32_t m_iImageWidth;
	int32_t m_iImageHeight;
	int32_t m_iComponents;
	char m_sFileExtension[255];
	bool m_bLoaded;

public:
	CImage() { m_bLoaded = false; };
	CImage(const char* path)
	{
		LoadPicture(path);
	}

	void LoadPicture(const char* path)
	{
		m_uImageData = stbi_load(path, &m_iImageWidth, &m_iImageHeight, &m_iComponents, 0);

		m_vImageChannelMatrix[CHANNEL_RED] = CMatrix2x2<double>(m_iImageWidth, m_iImageHeight);
		m_vImageChannelMatrix[CHANNEL_GREEN] = CMatrix2x2<double>(m_iImageWidth, m_iImageHeight);
		m_vImageChannelMatrix[CHANNEL_BLUE] = CMatrix2x2<double>(m_iImageWidth, m_iImageHeight);

		if (m_iComponents > 3)
			m_vImageChannelMatrix[CHANNEL_ALPHA] = CMatrix2x2<double>(m_iImageWidth, m_iImageHeight);

		SeparateChannels();
	
		strcpy(m_sFileExtension, PathFindExtensionA(path));
		if (m_uImageData)
			m_bLoaded = true;
	}

	void SeparateChannels()
	{
		double** RedChannel = m_vImageChannelMatrix[CHANNEL_RED].GetMatrix();
		double** GreenChannel = m_vImageChannelMatrix[CHANNEL_GREEN].GetMatrix();
		double** BlueChannel = m_vImageChannelMatrix[CHANNEL_BLUE].GetMatrix();
		double** AlphaChannel = m_vImageChannelMatrix[CHANNEL_ALPHA].GetMatrix();

		for (int h1 = 0; h1 < m_iImageHeight; h1++)
		{
			for (int w1 = 0; w1 < m_iImageWidth; w1++)
			{
				stbi_uc* pixel = m_uImageData + (m_iComponents * (h1 * m_iImageWidth + w1));
				RedChannel[h1][w1] = pixel[0];
				GreenChannel[h1][w1] = pixel[1];
				BlueChannel[h1][w1] = pixel[2];
				if (m_iComponents > 3)
					AlphaChannel[h1][w1] = pixel[3];

			}
		}

		DivideChannelsIntoBlocks();
	}

	void DivideChannelsIntoBlocks()
	{
		/*
			TODO:
			1. Account for width being smaller than height (use the largest one for the top loop)
			2. Account for width or height being undividable by 8 (reduce until they do)
		*/

		if (m_iImageWidth % 8 == 0 && m_iImageHeight % 8 == 0)
		{
			for (int h1 = 0; h1 < m_iImageHeight; h1 += 8)
			{
				for (int w1 = 0; w1 < m_iImageWidth; w1 += 8)
				{
					m_vChannelBlocks[CHANNEL_RED].push_back(CChannelBlock(m_vImageChannelMatrix[CHANNEL_RED].GetMatrixSegment(h1, w1, 8, 8), 8));
					m_vChannelBlocks[CHANNEL_GREEN].push_back(CChannelBlock(m_vImageChannelMatrix[CHANNEL_GREEN].GetMatrixSegment(h1, w1, 8, 8), 8));
					m_vChannelBlocks[CHANNEL_BLUE].push_back(CChannelBlock(m_vImageChannelMatrix[CHANNEL_BLUE].GetMatrixSegment(h1, w1, 8, 8), 8));
					if (m_iComponents > 3)
						m_vChannelBlocks[CHANNEL_ALPHA].push_back(CChannelBlock(m_vImageChannelMatrix[CHANNEL_ALPHA].GetMatrixSegment(h1, w1, 8, 8), 8));
				}
			}
		}
	}

	CMatrix2x2<double> GetImageChannel(eStreamChannels Channel)
	{
		return m_vImageChannelMatrix[Channel];
	}

	std::vector<CChannelBlock> GetImageChannelBlocks(eStreamChannels Channel)
	{
		return m_vChannelBlocks[Channel];
	}

	void UpdateImageDTC()
	{
		for (int CurrentChannel = CHANNEL_RED; CurrentChannel <= CHANNEL_ALPHA; CurrentChannel *= 0x02)
		{
			for (int i = 0; i < m_vChannelBlocks[(eStreamChannels)CurrentChannel].size(); i++)
			{
				m_vChannelBlocks[(eStreamChannels)CurrentChannel].at(i).UpdateDCT();
}
		}
	}

	void UpdateChannels()
	{
		UpdateImageDTC();
		double** RedChannel = m_vImageChannelMatrix[CHANNEL_RED].GetMatrix();
		double** GreenChannel = m_vImageChannelMatrix[CHANNEL_GREEN].GetMatrix();
		double** BlueChannel = m_vImageChannelMatrix[CHANNEL_BLUE].GetMatrix();
		double** AlphaChannel = m_vImageChannelMatrix[CHANNEL_ALPHA].GetMatrix();

		for (int h1 = 0; h1 < m_iImageHeight; h1++)
		{
			for (int w1 = 0; w1 < m_iImageWidth; w1++)
			{
				stbi_uc* pixel = m_uImageData + (m_iComponents * (h1 * m_iImageWidth + w1));
				pixel[0] = stbi__clamp(RedChannel[h1][w1]);
				pixel[1] = stbi__clamp(GreenChannel[h1][w1]);
				pixel[2] = stbi__clamp(BlueChannel[h1][w1]);
				if (m_iComponents > 3)
					pixel[3] = stbi__clamp(AlphaChannel[h1][w1]);
			}
		}
	}

	int32_t GetImageWidth()
	{
		if (!m_bLoaded)
			return 0;

		return m_iImageWidth;
	}

	int32_t GetImageHeight()
	{
		if (!m_bLoaded)
			return 0;

		return m_iImageHeight;
	}

	int32_t GetImageComponents()
	{
		if (!m_bLoaded)
			return 0;

		return m_iComponents;
	}

	int32_t GetImagePixelCount()
	{
		return m_iImageWidth * m_iImageHeight;
	}

	uint8_t* GetImageData()
	{
		if (!m_bLoaded)
			return nullptr;

		return m_uImageData;
	}

	const char* GetImageFormat()
	{
		return m_sFileExtension;
	}

	bool IsLoaded()
	{
		return m_bLoaded;
	}
};


class CBitContainer
{
private:
	union BitUnion
	{
		std::bitset<8> bits;
		uint8_t bytes;
		BitUnion() {};
	};

	std::vector<unsigned char> m_vBytes;
	BitUnion m_BitUnion;
	int32_t m_iBitShift;
	int16_t m_uBitQuantity = 1;
public:
	CBitContainer()
	{
		Clear();
	}

	void Clear()
	{
		m_vBytes.clear();
		ClearBitField();
	}

	void ClearBitField()
	{
		m_iBitShift = 0;
		for (int i = 0; i < m_BitUnion.bits.size() - 0; i++)
			m_BitUnion.bits.set(i, 0);
	}

	void SetBitQuantity(int16_t uBitQuantity)
	{
		m_uBitQuantity = uBitQuantity;
	}

	void ReadBits(uint8_t input)
	{
		input &= (0b11111111 >> (8 - m_uBitQuantity));
		m_BitUnion.bits |= input;
		m_iBitShift++;
		if (m_iBitShift == (8 / m_uBitQuantity))
		{
			PushBitsToBytes();
		}
		else
			m_BitUnion.bits <<= m_uBitQuantity;

	}

	void PushBitsToBytes()
	{
		m_vBytes.push_back(m_BitUnion.bytes);
		for (int i = 0; i < m_BitUnion.bits.size() - 0; i++)
			m_BitUnion.bits.set(i, 0);
		m_iBitShift = 0;
	}

	std::vector<unsigned char> GetBytes()
	{
		return m_vBytes;
	}
};

class CBitEncoder
{
private:
	uint8_t* m_pBuffer;
	uint16_t m_uBitSequenceSize;
	uint32_t m_iLength;
	uint32_t m_iCurrentByte;
	uint32_t m_iCurrentShift;
	bool	 m_bOutOfBits = false;
public:
	CBitEncoder() {};

	CBitEncoder(uint8_t* pBuf, uint32_t len, uint16_t uBitSequenceSize)
	{
		Init(pBuf, len, uBitSequenceSize);
	}

	void Init(uint8_t* pBuf, uint32_t len, uint16_t uBitSequenceSize)
	{
		m_iLength = len;
		m_pBuffer = new uint8_t[len];
		memcpy(m_pBuffer, pBuf, len);
		m_uBitSequenceSize = uBitSequenceSize;
		m_iCurrentByte = 0;
		m_iCurrentShift = 8 - m_uBitSequenceSize;
	}

	void Clear()
	{
		m_iLength = 0;
		if (m_pBuffer)
			delete[] m_pBuffer;
		m_iCurrentByte = 0;
		m_iCurrentShift = 8 - m_uBitSequenceSize;
	}

	uint32_t GetBitsSize()
	{
		return (m_iLength * 8) / m_uBitSequenceSize;
	}

	uint8_t GetNextBitSeq()
	{
		uint8_t BitSequence = m_pBuffer[m_iCurrentByte] >> m_iCurrentShift;
		BitSequence &= (0b11111111 >> (8 - m_uBitSequenceSize));
		m_iCurrentShift -= m_uBitSequenceSize;
		if (m_iCurrentShift == (-1 * m_uBitSequenceSize))
		{
			m_iCurrentShift = 8 - m_uBitSequenceSize;
			if (m_iCurrentByte < m_iLength)
				m_iCurrentByte++;

			else
				m_bOutOfBits = true;
		}

		return BitSequence;
	}
};

class CDataChunk
{
private:
	uint8_t* m_pData;
	uint32_t m_iDataSize;
	bool m_bReferenceMode = false;
public:
	CDataChunk()
	{
		m_pData = nullptr;
		m_iDataSize = 0;
	}

	CDataChunk(uint8_t* pData, uint32_t iDataSize)
	{
		ReferToData(pData, iDataSize);
	}

	void SaveData(uint8_t* pData, uint32_t iDataSize)
	{
		m_iDataSize = iDataSize;
		m_pData = new uint8_t[m_iDataSize];
		memcpy(m_pData, pData, iDataSize);
		m_bReferenceMode = false;
	}

	void SaveData()
	{
		if (!m_bReferenceMode)
			return;

		uint8_t* pStoreData = new uint8_t[m_iDataSize];
		memcpy(pStoreData, m_pData, m_iDataSize);
		m_pData = new uint8_t[m_iDataSize];
		memcpy(m_pData, pStoreData, m_iDataSize);
		delete[] pStoreData;
		m_bReferenceMode = false;
	}

	void Clear()
	{
		if (!m_bReferenceMode)
			delete[] m_pData;

		m_pData = nullptr;
		m_iDataSize = 0;
	}

	void ReferToData(uint8_t* pData, uint32_t iDataSize)
	{
		m_pData = pData;
		m_iDataSize = iDataSize;
		m_bReferenceMode = true;
	}

	uint32_t GetDataSize()
	{
		return m_iDataSize;
	}

	uint8_t* GetData()
	{
		return m_pData;
	}

	CDataChunk operator+ (CDataChunk cFirstChunk) const
	{
		uint8_t* pData = new uint8_t[m_iDataSize + cFirstChunk.GetDataSize()];
		int32_t iDataSize = m_iDataSize + cFirstChunk.GetDataSize();
		memcpy(pData, m_pData, m_iDataSize);
		memcpy(pData + m_iDataSize, cFirstChunk.GetData(), cFirstChunk.GetDataSize());
		CDataChunk cUnitedChunks(pData, iDataSize);
		cUnitedChunks.SaveData();
		return cUnitedChunks;
	}

	CDataChunk& operator= (CDataChunk cFirstChunk)
	{
		ReferToData(cFirstChunk.GetData(), cFirstChunk.GetDataSize());
		return *this;
	}

	CDataChunk& operator+= (CDataChunk cFirstChunk)
	{
		return *this = *this + cFirstChunk;
	}
};

class CFragmentedData
{
private:
	uint8_t* m_pData;
	std::vector<CDataChunk> m_pChunks;
	uint32_t m_iDataSize;
	uint32_t m_iChunkSize;
	uint32_t m_iReservedSpace;
	uint16_t m_uChunks;
public:
	CFragmentedData() {};

	CFragmentedData(uint8_t* pData, uint32_t iDataSize, uint16_t uChunks, uint32_t iReservedSpace)
	{
		Init(pData, iDataSize, uChunks, iReservedSpace);
	}

	void Init(uint8_t* pData, uint32_t iDataSize, uint16_t uChunks, uint32_t iReservedSpace)
	{
		m_pData = pData;
		m_iDataSize = iDataSize;
		m_uChunks = uChunks;
		m_iReservedSpace = iReservedSpace;
		FragmentData();
	}

	void Clear()
	{
		m_pData = nullptr;
		m_iDataSize = 0;
		m_uChunks = 0;
		m_iChunkSize = 0;
		m_pChunks.clear();
	}

	void FragmentData()
	{
		m_iChunkSize = (m_iDataSize - m_iReservedSpace) / m_uChunks;
		uint32_t iLeftData = ((m_iDataSize - m_iReservedSpace) % m_uChunks) + m_iReservedSpace;
		uint32_t iPointerShift = 0;
		m_pChunks.push_back(CDataChunk(m_pData, m_iChunkSize + iLeftData));
		iPointerShift += m_pChunks.at(0).GetDataSize();

		for (int i = 1; i < m_uChunks; i++)
		{
			m_pChunks.push_back(CDataChunk(m_pData + iPointerShift, m_iChunkSize));
			iPointerShift += m_pChunks.at(i).GetDataSize();
		}
	}

	std::vector<CDataChunk> GetDataChunks()
	{
		return m_pChunks;
	}
};

class CLeastSignificantBitsInterface
{
private:
	CBitContainer m_cBitStream;
	CBitEncoder m_cBitEncoder;
	CFragmentedData m_fDataFragments;
	uint32_t m_iSizeAvailablePerChannel;
	uint32_t m_iAvailableSpace;
	uint8_t* m_pDataBuffer = nullptr;
	uint8_t m_iTotalChannelsAvailable;
	const uint16_t m_uBitsPerPixel = 2;
	const uint16_t m_uMaxChannels = 4;
public:
	CLeastSignificantBitsInterface() { m_cBitStream.SetBitQuantity(2); };

	uint32_t CalculateChannelSize(int ImageWidth, int ImageHeight)
	{
		m_iSizeAvailablePerChannel = 0;
		m_iSizeAvailablePerChannel = ((ImageWidth * ImageHeight) * m_uBitsPerPixel) / 8;
		return m_iSizeAvailablePerChannel;
	}

	uint32_t CalculateAvailableSpace(uint8_t AvailableChannels)
	{
		m_iAvailableSpace = 0;
		m_iTotalChannelsAvailable = 0;
		if (AvailableChannels & CHANNEL_RED)
			m_iAvailableSpace += m_iSizeAvailablePerChannel;

		if (AvailableChannels & CHANNEL_GREEN)
			m_iAvailableSpace += m_iSizeAvailablePerChannel;

		if (AvailableChannels & CHANNEL_BLUE)
			m_iAvailableSpace += m_iSizeAvailablePerChannel;

		if (AvailableChannels & CHANNEL_ALPHA)
			m_iAvailableSpace += m_iSizeAvailablePerChannel;

		m_iTotalChannelsAvailable = m_iAvailableSpace / m_iSizeAvailablePerChannel;

		return m_iAvailableSpace;
	}

	bool StreamDoubletsToImage(CDataChunk cDataChunk, CImage& cImage, uint8_t AvailableChannels)
	{
		CalculateChannelSize(cImage.GetImageWidth(), cImage.GetImageHeight());
		CalculateAvailableSpace(AvailableChannels);

		if (m_iAvailableSpace < cDataChunk.GetDataSize())
			return false; // TODO: 1. Add an exception if there is not enough space to store our data.

		m_fDataFragments.Init(cDataChunk.GetData(), cDataChunk.GetDataSize(), m_iTotalChannelsAvailable, STEGPLUS_SIGNATURE_RESERVE);

		int32_t iCurrentChunk = 0;
		int32_t iCurrentChannel = 0;
		for (int CurrentChannel = CHANNEL_RED; CurrentChannel <= CHANNEL_ALPHA; CurrentChannel *= 0x02)
		{
			if (AvailableChannels & CurrentChannel)
			{
				CDataChunk CurrentChunk = m_fDataFragments.GetDataChunks().at(iCurrentChunk);
				StreamDoubletsToChannel(CurrentChunk.GetData(), CurrentChunk.GetDataSize(), cImage, iCurrentChannel);
				iCurrentChunk++;
			}
			iCurrentChannel++;
		}

		m_fDataFragments.Clear();
		return true;
	}

	void StreamDoubletsToChannel(uint8_t* pDataBuffer, uint32_t DataBufferSize, CImage cImage, int Channel)
	{
		int DoubletsWritten = 0;
		m_cBitEncoder.Init(pDataBuffer, DataBufferSize, 2);
		/*
			TODO:
			1. Reduce amount of interclass calls.
		*/

		for (int h1 = 0; h1 < cImage.GetImageHeight(); h1++)
		{
			for (int w1 = 0; w1 < cImage.GetImageWidth(); w1++)
			{
				DoubletsWritten++;
				stbi_uc* pixel = cImage.GetImageData() + (cImage.GetImageComponents() * (h1 * cImage.GetImageWidth() + w1));
				pixel[Channel] &= 0b11111100;
				pixel[Channel] |= m_cBitEncoder.GetNextBitSeq();
				if (DoubletsWritten > m_cBitEncoder.GetBitsSize()) // TODO: 1. Is this really necessary?
					return;
			}
		}
		m_cBitEncoder.Clear();
	}

	void GetDoubletsFromImage(CImage cImage, uint32_t DesiredDataSize, uint8_t AvailableChannels)
	{
		m_cBitStream.Clear();
		CalculateChannelSize(cImage.GetImageWidth(), cImage.GetImageHeight());
		CalculateAvailableSpace(AvailableChannels);

		int iChunkSize = (DesiredDataSize - STEGPLUS_SIGNATURE_RESERVE) / m_iTotalChannelsAvailable; // 2
		int iLeftOver = ((DesiredDataSize - STEGPLUS_SIGNATURE_RESERVE) % m_iTotalChannelsAvailable) + STEGPLUS_SIGNATURE_RESERVE; // 9
		std::vector<int> vChunkSizes;
		vChunkSizes.push_back(iChunkSize + iLeftOver);
		for (int i = 1; i < m_iTotalChannelsAvailable; i++)
			vChunkSizes.push_back(iChunkSize);

		int32_t iCurrentChunk = 0;
		int32_t iCurrentChannel = 0;

		for (int CurrentChannel = CHANNEL_RED; CurrentChannel <= CHANNEL_ALPHA; CurrentChannel *= 0x02) // redo loop
		{
			if (AvailableChannels & CurrentChannel)
			{
				GetDoubletsFromChannel(cImage, vChunkSizes.at(iCurrentChunk), iCurrentChannel);
				m_cBitStream.ClearBitField();
				iCurrentChunk++;
			}

			iCurrentChannel++;
		}
	}

	void GetDoubletsFromChannel(CImage cImage, uint32_t DesiredDataSize, int Channel)
	{
		int DoubletsRead = 0;
		/*
			TODO:
			1. Reduce the number of interclass calls.
		*/

		for (int h1 = 0; h1 < cImage.GetImageHeight(); h1++)
		{
			for (int w1 = 0; w1 < cImage.GetImageWidth(); w1++)
			{
				stbi_uc* pixel = cImage.GetImageData() + (cImage.GetImageComponents() * (h1 * cImage.GetImageWidth() + w1));
				m_cBitStream.ReadBits(pixel[Channel]);
				DoubletsRead++;
				if (DoubletsRead > (DesiredDataSize * 8) / 2)
					return;
			}
		}
	}

	CDataChunk GetDataFromDoublets()
	{
		//	TODO: 1. This somehow causes CLI to crash, fix that.
		if (m_pDataBuffer)
			delete[] m_pDataBuffer;
		m_pDataBuffer = nullptr;

		m_pDataBuffer = new uint8_t[m_cBitStream.GetBytes().size()];
		ZeroMemory(m_pDataBuffer, m_cBitStream.GetBytes().size());

		for (int i = 0; i < m_cBitStream.GetBytes().size(); i++)
			m_pDataBuffer[i] = m_cBitStream.GetBytes().at(i);

		CDataChunk DataBuffer(m_pDataBuffer, m_cBitStream.GetBytes().size());
		DataBuffer.SaveData();
		return DataBuffer;
	}

	std::vector<uint8_t> GetBytesFromDoublets()
	{
		return m_cBitStream.GetBytes();
	}
};

class CXORCrypt
{
private:
	char m_cXORKey[128];
public:
	CXORCrypt() { ZeroMemory(m_cXORKey, sizeof(m_cXORKey)); }

	void SetKey(const char* Key, int32_t iSize)
	{
		if (iSize > sizeof(m_cXORKey) || iSize < 0)
			return;

		memcpy(m_cXORKey, Key, iSize);
	}

	void xorData(uint8_t* pData, int32_t iDataSize)
	{
		for (int i = 0; i < iDataSize; i++)
			pData[i] = pData[i] ^ m_cXORKey[i % (sizeof(m_cXORKey) / sizeof(char))];
	}
};

class CDCTInterface
{
private:
	CBitContainer m_cBitStream;
	CBitEncoder m_cBitEncoder;
	CFragmentedData m_fDataFragments;
	uint32_t m_iSizeAvailablePerChannel;
	uint32_t m_iAvailableSpace;
	uint8_t* m_pDataBuffer = nullptr;
	uint8_t m_iTotalChannelsAvailable;
	const uint16_t m_uBitsPerBlock = 1;
	const uint16_t m_uBlockSize = 8;
	const uint16_t m_uPseudoRandomFactor = 200;
	const uint16_t m_uMaxChannels = 4;
public:
	CDCTInterface() { m_cBitStream.SetBitQuantity(1); }
	uint32_t CalculateChannelSize(int ImageWidth, int ImageHeight)
	{
		/*
			TODO: 1. Account for ImageWidth || ImageHeight being not dividable by 8.
		*/

		m_iSizeAvailablePerChannel = 0;
		m_iSizeAvailablePerChannel = (((ImageWidth * ImageHeight) * m_uBitsPerBlock) / 8) / m_uBlockSize;
		return m_iSizeAvailablePerChannel;
	}

	uint32_t CalculateAvailableSpace(uint8_t AvailableChannels)
	{
		m_iAvailableSpace = 0;
		m_iTotalChannelsAvailable = 0;
		if (AvailableChannels & CHANNEL_RED)
			m_iAvailableSpace += m_iSizeAvailablePerChannel;

		if (AvailableChannels & CHANNEL_GREEN)
			m_iAvailableSpace += m_iSizeAvailablePerChannel;

		if (AvailableChannels & CHANNEL_BLUE)
			m_iAvailableSpace += m_iSizeAvailablePerChannel;

		if (AvailableChannels & CHANNEL_ALPHA)
			m_iAvailableSpace += m_iSizeAvailablePerChannel;

		m_iTotalChannelsAvailable = m_iAvailableSpace / m_iSizeAvailablePerChannel;

		return m_iAvailableSpace;
	}

	void StreamBitsToChannel(uint8_t* pDataBuffer, uint32_t DataBufferSize, CImage cImage, eStreamChannels Channel)
	{
		m_cBitEncoder.Init(pDataBuffer, DataBufferSize, 1);
		std::vector<CChannelBlock> ChannelBlocks = cImage.GetImageChannelBlocks(Channel);
		int iBitsWritten = 0;
		for (int i = 0; i < ChannelBlocks.size(); i++)
		{
			double** mSingleBlock = ChannelBlocks.at(i).GetBlockDCT().GetMatrix();
			bool bSingleBit = m_cBitEncoder.GetNextBitSeq();
			float iFrequencyDifference = abs(abs(mSingleBlock[4][4]) - abs(mSingleBlock[5][5]));

			if (bSingleBit)
			{
			//	if (iFrequencyDifference >= m_uPseudoRandomFactor)
			//	{
				//if(mSingleBlock[4][4] >= mSingleBlock[5][5])
					mSingleBlock[5][5] = abs(mSingleBlock[4][4]) + m_uPseudoRandomFactor + 5;
			//	}
			}

			else
			{
			//	if (iFrequencyDifference < m_uPseudoRandomFactor)
			//	{
				//if (mSingleBlock[5][5] < mSingleBlock[4][4])
					mSingleBlock[4][4] = abs(mSingleBlock[5][5]) + m_uPseudoRandomFactor + 5;
			//	}
			}

			if (iBitsWritten > m_cBitEncoder.GetBitsSize()) // TODO: 1. Is this really necessary?
			{
				return;
			}

		}

		m_cBitEncoder.Clear();

	}

	bool StreamBitsToImage(CDataChunk cDataChunk, CImage& cImage, uint8_t AvailableChannels)
	{
		CalculateChannelSize(cImage.GetImageWidth(), cImage.GetImageHeight());
		CalculateAvailableSpace(AvailableChannels);

		if (m_iAvailableSpace < cDataChunk.GetDataSize())
			return false; // TODO: 1. Add an exception if there is not enough space to store our data.

		m_fDataFragments.Init(cDataChunk.GetData(), cDataChunk.GetDataSize(), m_iTotalChannelsAvailable, STEGPLUS_SIGNATURE_RESERVE);
		
		int32_t iCurrentChunk = 0;
		for (int CurrentChannel = CHANNEL_RED; CurrentChannel <= CHANNEL_ALPHA; CurrentChannel *= 0x02)
		{
			if (AvailableChannels & CurrentChannel)
			{
				CDataChunk CurrentChunk = m_fDataFragments.GetDataChunks().at(iCurrentChunk);
				StreamBitsToChannel(CurrentChunk.GetData(), CurrentChunk.GetDataSize(), cImage, (eStreamChannels)CurrentChannel);
				iCurrentChunk++;
			}
		}
		m_fDataFragments.Clear();
		return true;

	}

	void GetBitsFromImage(CImage cImage, uint32_t DesiredDataSize, uint8_t AvailableChannels)
	{
		m_cBitStream.Clear();
		CalculateChannelSize(cImage.GetImageWidth(), cImage.GetImageHeight());
		CalculateAvailableSpace(AvailableChannels);

		int iChunkSize = (DesiredDataSize - STEGPLUS_SIGNATURE_RESERVE) / m_iTotalChannelsAvailable; // 2
		int iLeftOver = ((DesiredDataSize - STEGPLUS_SIGNATURE_RESERVE) % m_iTotalChannelsAvailable) + STEGPLUS_SIGNATURE_RESERVE; // 9
		std::vector<int> vChunkSizes;
		vChunkSizes.push_back(iChunkSize + iLeftOver);
		for (int i = 1; i < m_iTotalChannelsAvailable; i++)
			vChunkSizes.push_back(iChunkSize);

		int32_t iCurrentChunk = 0;

		double** nSingleBlock = cImage.GetImageChannelBlocks(CHANNEL_RED).at(checkblock).GetBlockDCT().GetMatrix();
		for (int CurrentChannel = CHANNEL_RED; CurrentChannel <= CHANNEL_ALPHA; CurrentChannel *= 0x02) // redo loop
		{
			if (AvailableChannels & CurrentChannel)
			{
				GetBitsFromChannel(cImage, vChunkSizes.at(iCurrentChunk), (eStreamChannels)CurrentChannel);
				m_cBitStream.ClearBitField();
				iCurrentChunk++;
			}
		}
	}

	void GetBitsFromChannel(CImage cImage, uint32_t DesiredDataSize, eStreamChannels Channel)
	{
		std::vector<CChannelBlock> ChannelBlocks = cImage.GetImageChannelBlocks(Channel);
		int32_t BitsRead = 0;

		for (int i = 0; i < ChannelBlocks.size(); i++)
		{
			double** mSingleBlock = ChannelBlocks.at(i).GetBlockDCT().GetMatrix();
			float iFrequencyDifference = abs(mSingleBlock[4][4]) - abs(mSingleBlock[5][5]);

			if (mSingleBlock[4][4] <= mSingleBlock[5][5])
			{
				m_cBitStream.ReadBits(1);
			}

			if (mSingleBlock[4][4] > mSingleBlock[5][5])
			{
				m_cBitStream.ReadBits(0);
			}

			BitsRead++;
			if (BitsRead > (DesiredDataSize * 8))
				return;

		}
	}

	CDataChunk GetDataFromBits()
	{
		//	TODO: 1. This somehow causes CLI to crash, fix that.
		if (m_pDataBuffer)
			delete[] m_pDataBuffer;
		m_pDataBuffer = nullptr;

		m_pDataBuffer = new uint8_t[m_cBitStream.GetBytes().size()];
		ZeroMemory(m_pDataBuffer, m_cBitStream.GetBytes().size());

		for (int i = 0; i < m_cBitStream.GetBytes().size(); i++)
			m_pDataBuffer[i] = m_cBitStream.GetBytes().at(i);

		CDataChunk DataBuffer(m_pDataBuffer, m_cBitStream.GetBytes().size());
		DataBuffer.SaveData();
		return DataBuffer;
	}

};

class CSteganographyPlusData
{
private:
	CDataChunk m_cStegPlusData;
	CXORCrypt m_cEncryptionInterface;
	char m_cStegPlusSignature[5] = "STG+";
public:
	CSteganographyPlusData()
	{
		Clear();
	}

	void EncapsulateData(uint8_t* pData, int32_t iDataSize, bool bEncrypt = false, const char* cKey = "")
	{
		if (!bEncrypt)
		{
			m_cStegPlusData += CDataChunk((uint8_t*)&iDataSize, sizeof(int32_t)) + CDataChunk(pData, iDataSize);
		}

		else
		{
			CDataChunk DataSize, DataSegment;
			DataSize.ReferToData((uint8_t*)&iDataSize, sizeof(int32_t));
			DataSize.SaveData();

			DataSegment.ReferToData(pData, iDataSize);
			DataSegment.SaveData();

			m_cEncryptionInterface.SetKey(cKey, strlen(cKey));
			m_cEncryptionInterface.xorData(m_cStegPlusData.GetData(), m_cStegPlusData.GetDataSize());
			m_cEncryptionInterface.xorData(DataSize.GetData(), DataSize.GetDataSize());
			m_cEncryptionInterface.xorData(DataSegment.GetData(), DataSegment.GetDataSize());
			m_cStegPlusData += DataSize + DataSegment;
		}
	}

	bool CheckStegPlusSignature(uint8_t* pData)
	{
		if (pData[0] == 'S' && pData[1] == 'T' && pData[2] == 'G' && pData[3] == '+' && pData[4] == '\0')
			return true;

		return false;
	}

	bool DecapsulateData(CImage* cImage, void* cStegPlusInterface, uint8_t fChannels, eSteganographyModes iMode, bool bEncrypted, const char* cKey = "")
	{
		CLeastSignificantBitsInterface* LSBInterface = nullptr;
		CDCTInterface* DCTInterface = nullptr;
		uint8_t* pSig = nullptr;

		if (iMode == STEG_LSB)
		{
			LSBInterface = (CLeastSignificantBitsInterface*)cStegPlusInterface;
			LSBInterface->GetDoubletsFromImage(*cImage, STEGPLUS_SIGNATURE_RESERVE, fChannels);
			pSig = LSBInterface->GetDataFromDoublets().GetData();
		}

		if (iMode == STEG_DCT)
		{
			DCTInterface = (CDCTInterface*)cStegPlusInterface;
			DCTInterface->GetBitsFromImage(*cImage, STEGPLUS_SIGNATURE_RESERVE, fChannels);
			pSig = DCTInterface->GetDataFromBits().GetData();
		}

		if (bEncrypted)
		{

			m_cEncryptionInterface.SetKey(cKey, strlen(cKey));
			m_cEncryptionInterface.xorData(pSig, 5);
			m_cEncryptionInterface.xorData(pSig + 5, 4);
		}

		if (CheckStegPlusSignature(pSig))
		{

			int32_t iRealSize = *(int32_t*)(pSig + 5);
			uint8_t* pData = nullptr;
			if (iMode == STEG_LSB)
			{
				LSBInterface->GetDoubletsFromImage(*cImage, iRealSize + STEGPLUS_SIGNATURE_RESERVE, fChannels);
				pData = LSBInterface->GetDataFromDoublets().GetData() + STEGPLUS_SIGNATURE_RESERVE;
			}

			if (iMode == STEG_DCT)
			{
				DCTInterface->GetBitsFromImage(*cImage, iRealSize + STEGPLUS_SIGNATURE_RESERVE, fChannels);
				pData = DCTInterface->GetDataFromBits().GetData() + STEGPLUS_SIGNATURE_RESERVE;
			}

			if (bEncrypted)
				m_cEncryptionInterface.xorData(pData, iRealSize);

			m_cStegPlusData = CDataChunk(pData, iRealSize);
			return true;
		}
		else
		{
			return false; // TODO: Add an exception if the sig is not found
		}
	}

	CDataChunk GetDataChunk()
	{
		return m_cStegPlusData;
	}

	void Clear()
	{
		m_cStegPlusData.Clear();
		m_cStegPlusData.ReferToData((uint8_t*)m_cStegPlusSignature, 5);
	}
};

class CSteganographyPlus
{
private:
	CLeastSignificantBitsInterface m_cLSBInterface;
	CDCTInterface m_cDCTInterface;
	CSteganographyPlusData m_cStegPlusData;
	CImage m_cImage;
	bool m_bUseEncryption;
public:
	CSteganographyPlus()
	{
		m_bUseEncryption = false;
	}

	CImage GetImage()
	{
		return m_cImage;
	}

	void LoadPicture(const char* path)
	{
		m_cImage.LoadPicture(path);
	}

	void EncodeData(uint8_t* pData, int32_t iDataSize, uint8_t fChannels, eSteganographyModes iMode = STEG_LSB, bool bEncrypt = false, const char* cKey = "")
	{
		if (!m_cImage.IsLoaded())
			return;

		if (fChannels & CHANNEL_ALPHA && m_cImage.GetImageComponents() < 4)
			return;

		int32_t iUsedChannels = 0;
		for (int i = 0x1; i <= CHANNEL_ALPHA; i *= 0x2)
			if (fChannels & i)
				iUsedChannels++;

		m_cStegPlusData.EncapsulateData(pData, iDataSize, bEncrypt, cKey);
		if (iMode == STEG_LSB)
		{
			if (iDataSize + STEGPLUS_SIGNATURE_RESERVE > (((m_cImage.GetImageWidth() * m_cImage.GetImageHeight()) * 8) / 2) * iUsedChannels)
			{
				printf("Not enough space!\n");
				return;
			}

			m_cLSBInterface.StreamDoubletsToImage(m_cStegPlusData.GetDataChunk(), m_cImage, fChannels);
		}

		if (iMode == STEG_DCT)
		{
			if (iDataSize + STEGPLUS_SIGNATURE_RESERVE > ((m_cImage.GetImageWidth() * m_cImage.GetImageHeight()) * 8) * iUsedChannels)
			{
				printf("Not enough space!\n");
				return;
			}
			m_cDCTInterface.StreamBitsToImage(m_cStegPlusData.GetDataChunk(), m_cImage, fChannels);
		}
		m_cStegPlusData.Clear();
	}

	void DecodeData(uint8_t fChannels, eSteganographyModes iMode = STEG_LSB, bool bEncrypt = false, const char* cKey = "")
	{
		if (!m_cImage.IsLoaded())
			return;

		if (fChannels & CHANNEL_ALPHA && m_cImage.GetImageComponents() < 4)
			return;

		if (iMode == STEG_LSB)
			m_cStegPlusData.DecapsulateData(&m_cImage, &m_cLSBInterface, fChannels, iMode, bEncrypt, cKey);

		if (iMode == STEG_DCT)
			m_cStegPlusData.DecapsulateData(&m_cImage, &m_cDCTInterface, fChannels, iMode, bEncrypt, cKey);

	}

	CDataChunk GetDataChunk()
	{
		return m_cStegPlusData.GetDataChunk();
	}

	void SavePicture(const char* path)
	{
		if (strcmp(m_cImage.GetImageFormat(), ".png") == 0)
		{
			stbi_write_png(path, m_cImage.GetImageWidth(), m_cImage.GetImageHeight(), m_cImage.GetImageComponents(), m_cImage.GetImageData(), 0);
		}

		else if (strcmp(m_cImage.GetImageFormat(), ".jpg") == 0)
		{
			stbi_write_jpg(path, m_cImage.GetImageWidth(), m_cImage.GetImageHeight(), m_cImage.GetImageComponents(), m_cImage.GetImageData(), 0);
		}
		
		else if (strcmp(m_cImage.GetImageFormat(), ".bmp") == 0)
		{
			stbi_write_bmp(path, m_cImage.GetImageWidth(), m_cImage.GetImageHeight(), m_cImage.GetImageComponents(), m_cImage.GetImageData());
		}

		else if (strcmp(m_cImage.GetImageFormat(), ".tga") == 0)
		{
			stbi_write_tga(path, m_cImage.GetImageWidth(), m_cImage.GetImageHeight(), m_cImage.GetImageComponents(), m_cImage.GetImageData());
		}

		else
		{
			printf("Unknown image format.\n");
		}
	}
};


void main()
{
	using std::chrono::high_resolution_clock;
	using std::chrono::duration_cast;
	using std::chrono::duration;
	using std::chrono::milliseconds;
	using std::chrono::seconds;
	auto t1 = high_resolution_clock::now();
	CSteganographyPlus StegPlus;
	CSteganographyPlus StegPlus2;
	CDiscreteCosineTransform cdct;
	uint8_t Channels =  CHANNEL_RED | CHANNEL_GREEN | CHANNEL_BLUE;
	uint8_t StegMethod = STEG_DCT;

	StegPlus.LoadPicture("anton.png");
	StegPlus.EncodeData((uint8_t*)"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30", strlen("1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30") + 1, (eStreamChannels)Channels, (eSteganographyModes)StegMethod, true, "1337");
	StegPlus.GetImage().UpdateChannels();
	StegPlus.SavePicture("anton.png");
	
	StegPlus2.LoadPicture("anton.png");
	StegPlus2.DecodeData((eStreamChannels)Channels, (eSteganographyModes)StegMethod, true, "1337");

	auto t2 = high_resolution_clock::now();
	auto ms_int = duration_cast<seconds>(t2 - t1);
	duration<double, std::milli> ms_double = t2 - t1;
	std::cout << ms_int.count() << "sec\n";
	printf("%s\n", ((const char*)StegPlus2.GetDataChunk().GetData()));
	system("pause");
}