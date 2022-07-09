// This sample shows how to write a simple unit test for a function,
// using Google C++ testing framework.
// It is based on https://github.com/google/googletest/blob/main/docs/index.md

#include <limits.h>
#include "dy4.h"
#include "fourier.h"
#include "genfunc.h"
#include "gtest/gtest.h"

namespace {

class IDFT_Fixture: public ::testing::Test {

	public:

	IDFT_Fixture( ) {
	}

	void SetUp( ) {
		generateRandomSamples(data, no_samples, max_value, precision);
	}

	void TearDown( ) {
	}

	~IDFT_Fixture( )  {
	}

	unsigned int no_samples = NFFT;
	unsigned short int max_value = 10;
	unsigned char precision = 2;
	std::vector<float> data;
	std::vector<std::complex<float>> dataFreq;
	std::vector<std::complex<float>> dataInverse;
};

TEST_F(IDFT_Fixture, DISABLED_DFT_IDFT_Equal) {

	DFT(data, dataFreq);
	IDFT(dataFreq, dataInverse);

	ASSERT_EQ(data.size(), dataInverse.size()) << "Original/inverse vectors are of unequal length";

	for (unsigned int i = 0; i < data.size(); ++i) {
		EXPECT_EQ(data[i], std::real(dataInverse[i])) << "Original/inverse vectors differ at index " << i;
	}
}

TEST_F(IDFT_Fixture, DFT_IDFT_NEAR){

	DFT(data, dataFreq);
	IDFT(dataFreq, dataInverse);

	ASSERT_EQ(data.size(), dataInverse.size()) << "Original/inverse vectors are of unequal length";

	for (unsigned int i = 0; i < data.size(); ++i) {
		EXPECT_NEAR(data[i], std::real(dataInverse[i]), 10-2) << "Original/inverse vectors differ at index " << i;
	}
}


} // end of namespace
