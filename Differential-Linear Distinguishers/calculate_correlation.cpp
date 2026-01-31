#include <iostream>
#include <bitset>
#include <vector>
#include <random>
#include <cassert>
#include <ctime>
#include <omp.h>  // For parallel programming using OpenMP
#include <fstream>
using namespace std;

// Thread-local Mersenne Twister random number generator (64-bit)
thread_local std::mt19937_64 rng(std::random_device{}());


const int KEY_SIZE = 128;
const int TWEAK_SIZE = 64;
const int BLOCK32_SIZE = 32; 
const int BLOCK40_SIZE = 40;


// Helper function to extract a range of bits from a bitset
// Parameters:
//   bs: The source bitset
//   start: Starting position of the bits to extract
//   len: Number of bits to extract
// Returns: Extracted bits as a uint64_t value
template<size_t N>
uint64_t extract_bits(const bitset<N>& bs, size_t start, size_t len)
{
    if (start + len > N)
        throw out_of_range("Bit range out of bounds");
    uint64_t value = 0;
    for (size_t i = 0; i < len; i++)
        value |= static_cast<uint64_t>(bs[start + i]) << i;
    return value;
}


// Helper function to create a bitset from a uint64_t value
// Parameters:
//   value: The source uint64_t value
//   bits: Number of bits to use from the value (default: N)
// Returns: A bitset<N> containing the specified bits
template<size_t N>
bitset<N> create_bitset(uint64_t value, size_t bits = N)
{
    bitset<N> bs;
    for (size_t i = 0; i < min<size_t>(bits, 64); i++)
        bs[i] = (value >> i) & 1;
    return bs;
}


struct LinearParams {
    int alpha;
    int beta0;
    int beta1;
    int beta2;
};
LinearParams L32_params = { 11, 5, 9, 12 };
LinearParams L32_prime_params = { 11, 1, 26, 30 };
LinearParams L40_params = { 17, 1, 9, 30 };
LinearParams L64_params = { 3, 1, 26, 50 };
LinearParams L128_params = { 17, 7, 11, 14 };



// ChiChi non-linear transformation for even-dimensional bitsets (Definition 2)
// Input: x - bitset to transform (size N must be even)
// Output: Transformed bitset with specialized rules for different bit positions
template<size_t N>
bitset<N> ChiChi(const bitset<N>& x)
{
    const size_t m = N / 2;  // Half the size of the bitset

    bitset<N> y;
    // Transform bits in first segment (0 to m-4)
    for (size_t i = 0; i < m - 3; i++)
        y[i] = x[i] ^ (~x[i + 1] & x[i + 2]);
    // Transform bits in second segment (m+1 to N-3)
    for (size_t i = m + 1; i < N - 2; i++)
        y[i] = x[i] ^ (~x[i + 1] & x[i + 2]);

    // Special transformation rules for boundary bits
    y[m - 3] = x[m] ^ (~x[m - 2] & x[0]);
    y[m - 2] = x[m - 1] ^ (~x[0] & x[1]);
    y[m - 1] = ~x[m - 3] ^ (~x[m] & ~x[m + 1]);
    y[m] = x[m - 2] ^ (~x[m + 1] & x[m + 2]);
    y[N - 2] = x[N - 2] ^ (~x[N - 1] & x[m - 1]);
    y[N - 1] = x[N - 1] ^ (~x[m - 1] & x[m]);
    return y;
}


// Linear layer transformation
// Inputs:
//   x - bitset to transform
//   params - LinearParams defining the transformation rules
// Output: Transformed bitset where y[i] = x[idx0] ^ x[idx1] ^ x[idx2]
//         with idx0/1/2 = (alpha*i + beta0/1/2) % N
template<size_t N>
bitset<N> linear_layer(const bitset<N>& x, const LinearParams& params)
{
    bitset<N> y;
    for (size_t i = 0; i < N; i++)
    {
        size_t idx0 = (params.alpha * i + params.beta0) % N;
        size_t idx1 = (params.alpha * i + params.beta1) % N;
        size_t idx2 = (params.alpha * i + params.beta2) % N;
        y[i] = x[idx0] ^ x[idx1] ^ x[idx2];
    }
    return y;
}



// ChiLow-(32 + τ) decryption function
// Inputs:
//   ROUNDS - Number of encryption rounds
//   X1, X2 - Ciphertext blocks to decrypt
// Output: Decrypted plaintext block
bitset<BLOCK32_SIZE> Chilow32_decrypt(int ROUNDS, bitset<BLOCK32_SIZE> X1, bitset<BLOCK32_SIZE> X2)
{
    uniform_int_distribution<uint32_t> dist32(0, UINT32_MAX);
    for (int i = 0; i < ROUNDS - 1; i++)
    {
        // Apply non-linear ChiChi transformation
        X1 = ChiChi<BLOCK32_SIZE>(X1);
        X2 = ChiChi<BLOCK32_SIZE>(X2);

        // Apply linear layer transformation
        X1 = linear_layer<BLOCK32_SIZE>(X1, L32_params);
        X2 = linear_layer<BLOCK32_SIZE>(X2, L32_params);

        // Generate random tweak and XOR with both blocks
        bitset<32> T_state(dist32(rng));
        for (int j = 0; j < 32; j++)
        {
            X1[j] = X1[j] ^ T_state[j];
            X2[j] = X2[j] ^ T_state[j];
        }
    }

    // Final round transformations without tweak
    X1 = ChiChi<BLOCK32_SIZE>(X1);
    X2 = ChiChi<BLOCK32_SIZE>(X2);

    // Combine X1 and X2 by XOR and apply final linear layer
    bitset<BLOCK32_SIZE> X;
    for (int i = 0; i < 32; i++)
        X[i] = X1[i] ^ X2[i];

    X = linear_layer<BLOCK32_SIZE>(X, L32_params);

    return X;
}



// ChiLow-(40) decryption function
// Inputs:
//   ROUNDS - Number of encryption rounds
//   X1, X2 - Ciphertext blocks to decrypt
// Output: Decrypted plaintext block
bitset<BLOCK40_SIZE> Chilow40_decrypt(int ROUNDS, bitset<BLOCK40_SIZE> X1, bitset<BLOCK40_SIZE> X2)
{
    uniform_int_distribution<uint64_t> dist40(0, 1099511627775);
    for (int i = 0; i < ROUNDS - 1; i++)
    {
        // Apply non-linear ChiChi transformation
        X1 = ChiChi<BLOCK40_SIZE>(X1);
        X2 = ChiChi<BLOCK40_SIZE>(X2);

        // Apply linear layer transformation
        X1 = linear_layer<BLOCK40_SIZE>(X1, L40_params);
        X2 = linear_layer<BLOCK40_SIZE>(X2, L40_params);

        // Generate random tweak and XOR with both blocks
        bitset<40> T_state(dist40(rng));
        for (int j = 0; j < 40; j++)
        {
            X1[j] = X1[j] ^ T_state[j];
            X2[j] = X2[j] ^ T_state[j];
        }
    }

    // Final round transformations without tweak
    X1 = ChiChi<BLOCK40_SIZE>(X1);
    X2 = ChiChi<BLOCK40_SIZE>(X2);

    // Combine X1 and X2 by XOR and apply final linear layer
    bitset<BLOCK40_SIZE> X;
    for (int i = 0; i < 40; i++)
        X[i] = X1[i] ^ X2[i];

    X = linear_layer<BLOCK40_SIZE>(X, L40_params);
    return X;
}


// Function to calculate correlation for 32-bit blocks
// Parameters:
//   nr - Number of rounds
//   count - Number of test cases
//   Delta - Bit positions for initial vector (IV)
//   lambda - Bit positions to check in plaintext
void D32(int nr, int count, vector<int> Delta, vector<int> lambda)
{
    int counts = 0;
    bitset<32> IV(0);
    for (auto m : Delta)
        IV.set(m);

    clock_t start = clock();
#pragma omp parallel for
    for (int test = 0; test < count; test++)
    {
        uniform_int_distribution<uint32_t> dist32(0, UINT32_MAX);

        bitset<32> C(dist32(rng));
        bitset<32> C2;
        for (int j = 0; j < 32; j++)
            C2[j] = C[j] ^ IV[j];

        bitset<32> plaintext = Chilow32_decrypt(nr, C, C2);

        char t = 0;
        for (auto m : lambda)
            t = t ^ plaintext[m];

#pragma omp critical
        {
            if (t == 0)
                counts++;
            else
                counts--;
        }
    }
    clock_t end = clock();

    double elapsed_seconds = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    cout << "Elapsed time: " << elapsed_seconds << " s\n";

    double correlation = counts / (double)count;
    cout << "Correlation: " << correlation << "  =  2^{" << log2(abs(counts)) - log2(count) << "}" << endl;
};




// Function to calculate correlation for 40-bit blocks
// Parameters:
//   nr - Number of rounds
//   count - Number of test cases
//   Delta - Bit positions for initial vector (IV)
//   lambda - Bit positions to check in plaintext
void D40(int nr, int count, vector<int> Delta, vector<int> lambda)
{
    int counts = 0;
    bitset<40> IV(0);
    for (auto m : Delta)
        IV.set(m);

    clock_t start = clock();
#pragma omp parallel for
    for (int test = 0; test < count; test++)
    {
        uniform_int_distribution<uint64_t> dist40(0, 1099511627775);


        bitset<40> C(dist40(rng));
        bitset<40> C2;
        for (int j = 0; j < 40; j++)
            C2[j] = C[j] ^ IV[j];

        bitset<40> plaintext = Chilow40_decrypt(nr, C, C2);

        char t = 0;
        for (auto m : lambda)
            t = t ^ plaintext[m];

#pragma omp critical
        {
            if (t == 0)
                counts++;
            else
                counts--;
        }
    }
    clock_t end = clock();

    double elapsed_seconds = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    cout << "Elapsed time: " << elapsed_seconds << " s\n";

    double correlation = counts / (double)count;
    cout << "Correlation: " << correlation << "  =  2^{" << log2(abs(counts)) - log2(count) << "}" << endl;
};




int main()
{

    int count = 1073741824;  // Number of test cases (2^30)


    // Test case 1: 32-bit block
    // Delta: 0x01403000 (bit positions 7,9,18,19)
    // lambda: 0x00000804 (bit positions 20,29)
    // Expected correlation: 2^-10.39

    cout << "Calculate the correlation of the differential-linear pair in the middle part of D32. " << endl;
    cout << "Delta:  0x01403000     lambda:  0x00000804" << endl;
    vector<int> Delta1 = { 7,9,18,19 };
    vector<int> lambda1 = { 20,29 };
    D32(3, count, Delta1, lambda1);
    cout << endl << endl;


    // Test case 2: 40-bit block
    // Delta: 0x0000010104 (bit positions 23,31,37)
    // lambda: 0x0001000100 (bit positions 15,31)
    // Expected correlation: 2^-11.17
    cout << "Calculate the correlation of the differential-linear pair in the middle part of D40. " << endl;
    cout << "Delta:  0x0000010104   lambda:  0x0001000100" << endl;
    vector<int> Delta2 = { 23, 31, 37 };
    vector<int> lambda2 = { 15, 31 };
    D40(3, count, Delta2, lambda2);
    cout << endl << endl;

    return 0;

}


