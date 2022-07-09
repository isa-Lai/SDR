// This sample shows how to write a simple unit test for a function,
// using Google C++ testing framework.
// It is based on https://github.com/google/googletest/blob/main/docs/index.md

#include <limits.h>
#include "dy4.h"
#include "fourier.h"
#include "genfunc.h"
#include "gtest/gtest.h"

namespace {

class FFT_Fixture: public ::testing::Test {

	public:

	FFT_Fixture( ) {
	}

	void SetUp( ) {
		generateRandomSamples(data, no_samples, max_value, precision);

		twiddles.resize(int(NFFT/2));
		compute_twiddles(twiddles);

		dataComplex.resize(data.size());
		for (int i=0; i<data.size(); i++) {
			std::complex<float> val(data[i], 0.0);
			dataComplex[i] = val;
		}

		dataFourierOriginal.resize(data.size());
		dataFourierCompare.resize(data.size());
	}

	void TearDown( ) {
	}

	~FFT_Fixture( )  {
	}

	unsigned int no_samples = NFFT;
	unsigned short int max_value = 10;
	unsigned char precision = 2;
	float epsilon = 10-2;
	std::vector<float> data;
	std::vector<std::complex<float>> dataComplex;
	std::vector<std::complex<float>> dataFourierOriginal;
	std::vector<std::complex<float>> dataFourierCompare;
	std::vector<std::complex<float>> twiddles;
};

TEST_F(FFT_Fixture, DFT_FFT_NEAR){

	DFT(data, dataFourierOriginal);
	FFT_recursive(dataComplex, dataFourierCompare);

	ASSERT_EQ(dataFourierOriginal.size(), dataFourierCompare.size()) << "The two Fourier vectors are of unequal length";

	for (unsigned int i = 0; i < dataFourierOriginal.size(); ++i) {
		EXPECT_NEAR(std::abs(dataFourierOriginal[i]), std::abs(dataFourierCompare[i]), epsilon) << "Original/inverse vectors differ at index " << i;
	}

}

TEST_F(FFT_Fixture, FFT_improved_NEAR){

	FFT_recursive(dataComplex, dataFourierOriginal);
	FFT_improved(dataComplex, dataFourierCompare, twiddles, 1);

	ASSERT_EQ(dataFourierOriginal.size(), dataFourierCompare.size()) << "The two Fourier vectors are of unequal length";

	for (unsigned int i = 0; i < dataFourierOriginal.size(); ++i) {
		EXPECT_NEAR(std::abs(dataFourierOriginal[i]), std::abs(dataFourierCompare[i]), epsilon) << "Original/inverse vectors differ at index " << i;
	}

}

TEST_F(FFT_Fixture, FFT_optimized_NEAR){

	FFT_recursive(dataComplex, dataFourierOriginal);
	FFT_optimized(dataComplex, dataFourierCompare, twiddles);

	ASSERT_EQ(dataFourierOriginal.size(), dataFourierCompare.size()) << "The two Fourier vectors are of unequal length";

	for (unsigned int i = 0; i < dataFourierOriginal.size(); ++i) {
		EXPECT_NEAR(std::abs(dataFourierOriginal[i]), std::abs(dataFourierCompare[i]), epsilon) << "Original/inverse vectors differ at index " << i;
	}

}

} // end of namespace
