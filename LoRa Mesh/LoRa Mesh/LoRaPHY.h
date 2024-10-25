#ifndef LoRaPHY_h
#define LoRaPHY_h

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <bitset>
#include <algorithm>
#include <fftw3.h>

class LoRaPHY {
public:
    // Constructor
    LoRaPHY(double rf, int spreading_factor, double bandwidth, double sampling_freq)
        : rf_freq(rf), sf(spreading_factor), bw(bandwidth), fs(sampling_freq),
          cr(1), payload_len(0), has_header(true), crc(true), ldr(false), preamble_len(8) {
        initWhiteningSequence();
        
        // Calculate sample_num based on BW, Fs, and SF
        sample_num = static_cast<int>((fs / bw) * pow(2, sf));

        // Calculate fft_len based on sample_num (use next power of two if needed)
        fft_len = nextPowerOfTwo(sample_num);

        // Initialize chirps
        downchirp.resize(sample_num);
        upchirp.resize(sample_num);

        // Fill downchirp and upchirp with the appropriate values
        initializeChirps();
    }

    // Method to encode a message
    std::vector<int> encode(const std::string& message) {
        std::vector<int> encoded_data;
        for (char ch : message) {
            int byte = static_cast<int>(ch);
            int gray_encoded_byte = grayEncode(encodeByte(byte));
            encoded_data.push_back(gray_encoded_byte);
        }
        return encoded_data;
    }

    // Method to modulate encoded symbols
    std::vector<Sample> modulate(const std::vector<int>& symbols) {
        // Calculate the total number of samples needed
        size_t total_samples = symbols.size() * sample_num;

        // Resize the samples vector to the correct length
        std::vector<Sample> samples(total_samples);

        // Fill the samples vector with modulated data for each symbol
        for (size_t i = 0; i < symbols.size(); ++i) {
            int symbol = symbols[i];

            // Generate the chirp for the symbol and place it in the samples vector
            std::vector<Sample> chirp = generateChirp(symbol, true); // Assuming true means upchirp

            // Copy the chirp into the samples vector at the correct position
            std::copy(chirp.begin(), chirp.end(), samples.begin() + i * sample_num);
        }

        return samples;
    }

//    // Function to generate an upchirp and downchirp for alignment and demodulation
//    std::vector<Sample> generateUpChirp(int samples_per_symbol) {
//        std::vector<Sample> upchirp = generateChirp(samples_per_symbol, 1);
//        return upchirp;
//    };
//    std::vector<Sample> generateDownChirp(int samples_per_symbol) {
//        std::vector<Sample> downchirp = generateChirp(samples_per_symbol, 0);
//        return downchirp;
//    };

    // Function to convert the decoded symbols back into a character array
    char* convertToCharArray(const std::vector<int>& decoded_data) {
        char* message = new char[decoded_data.size() + 1];
        for (size_t i = 0; i < decoded_data.size(); ++i) {
            message[i] = static_cast<char>(decoded_data[i]);
        }
        message[decoded_data.size()] = '\0';
        return message;
    }
        
    std::vector<int> grayCoding(const std::vector<int>& din) {
        std::vector<int> symbols(din.size());
        for (size_t i = 0; i < 8 && i < din.size(); ++i) {
            symbols[i] = std::floor(din[i] / 4);
        }
        return symbols;
    }

    std::vector<int> diagDeinterleave(const std::vector<int>& symbols, int ppm) {
        std::vector<int> codewords;
        for (int symbol : symbols) {
            std::vector<int> binary = intToBinary(symbol, ppm);
            codewords.insert(codewords.end(), binary.begin(), binary.end());
        }
        return codewords;
    }

    std::vector<int> hammingDecode(std::vector<int>& codewords, int rdd) {
        std::vector<int> nibbles;
        int p1 = bitReduce(codewords, {8, 4, 3, 1});
        int p2 = bitReduce(codewords, {7, 6, 3, 2});
        int p4 = bitReduce(codewords, {5, 6, 3, 2});
        int p8 = bitReduce(codewords, {5, 4, 3, 2});

        int error_pos = p1 + 2 * p2 + 4 * p4 + 8 * p8;
        if (error_pos > 0 && error_pos <= codewords.size()) {
            codewords[error_pos - 1] ^= 1;
        }

        int nibble = ((codewords[2] >> 1) & 1) +
                     ((codewords[3] >> 3) & 1) +
                     ((codewords[4] >> 5) & 1) +
                     ((codewords[5] >> 7) & 1);
        nibbles.push_back(nibble);

        return nibbles;
    }

    std::vector<int> dewhiten(const std::vector<int>& bytes) {
        std::vector<int> bytes_w(bytes.size());
        for (size_t i = 0; i < bytes.size(); ++i) {
            bytes_w[i] = bytes[i] ^ whitening_seq[i];
        }
        return bytes_w;
    }

    int dechirp(int x, bool is_up, std::vector<Sample>& signal) {
        std::vector<Sample> chirp = is_up ? upchirp : downchirp;

        if (signal.size() < x + sample_num || chirp.size() < sample_num) {
            std::cerr << "Insufficient data in signal or chirp vectors!" << std::endl;
            return -1;
        }

        std::vector<Sample> sym(signal.begin() + x, signal.begin() + x + sample_num);
        std::vector<Sample> dechirped(sample_num);

        for (int i = 0; i < sample_num; ++i) {
            dechirped[i].I = sym[i].I * chirp[i].I + sym[i].Q * chirp[i].Q;
            dechirped[i].Q = sym[i].Q * chirp[i].I - sym[i].I * chirp[i].Q;
        }

        fftw_complex* in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fft_len);
        fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fft_len);
        if (in == nullptr || out == nullptr) {
            std::cerr << "Memory allocation for FFTW arrays failed!" << std::endl;
            if (in) fftw_free(in);
            if (out) fftw_free(out);
            return -1;
        }

        fftw_plan plan = fftw_plan_dft_1d(fft_len, in, out, FFTW_FORWARD, FFTW_MEASURE);
        if (plan == nullptr) {
            std::cerr << "FFTW plan creation failed!" << std::endl;
            fftw_free(in);
            fftw_free(out);
            return -1;
        }

        for (int i = 0; i < sample_num; ++i) {
            in[i][0] = dechirped[i].I;
            in[i][1] = dechirped[i].Q;
        }
        for (int i = sample_num; i < fft_len; ++i) {
            in[i][0] = 0.0;
            in[i][1] = 0.0;
        }

        fftw_execute(plan);

        std::vector<double> fft_result(fft_len);
        for (int i = 0; i < fft_len; ++i) {
            fft_result[i] = std::sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]);
        }

        auto max_it = std::max_element(fft_result.begin(), fft_result.end());
        int pk = std::distance(fft_result.begin(), max_it);

        fftw_destroy_plan(plan);
        fftw_free(in);
        fftw_free(out);

        // Normalize the peak index to the LoRa symbol range
        pk = (pk * (1 << sf)) / fft_len;

        return pk;
    }

    std::pair<std::vector<int>, bool> decode(const std::vector<int>& symbols) {
        std::vector<int> data;
        bool checksum_valid = false;

        // Step 1: Gray coding
        std::vector<int> gray_coded_symbols = grayCoding(symbols);

        // Step 2: Diagonal deinterleaving
        std::vector<int> codewords = diagDeinterleave(gray_coded_symbols, 8);

        // Step 3: Hamming decoding
        std::vector<int> nibbles = hammingDecode(codewords, cr);

        // Step 4: Dewhitening the decoded data
        data = dewhiten(nibbles);

        // Step 5: Validate checksum
        checksum_valid = validateChecksum(data);

        return std::make_pair(data, checksum_valid);
    }

    std::vector<int> demodulate(std::vector<Sample>& signal) {
        std::vector<int> symbols;
        // Ensure the size is at least 4096 elements
        if (signal.size() < 4096) {
            // Resize the vector and fill new elements with zeros
            signal.resize(4096, Sample{0.0f, 0.0f});
        }
        // Ensure we have enough signal data for demodulation
        if (signal.size() < sample_num) {
            std::cerr << "Insufficient signal size for demodulation!" << std::endl;
            return symbols;
        }

        for (int i = 0; i + sample_num <= signal.size(); i += sample_num) {
            int pk = dechirp(i, false, signal);
            if (pk >= 0 && pk < (1 << sf)) { // Ensure pk is within the valid range
                std::cout << "Symbol at index " << i << ": " << pk << std::endl; // Debug output
                symbols.push_back(pk);
            } else {
                std::cerr << "Error: Invalid symbol value " << pk << " at index " << i << std::endl;
            }
        }

        return symbols;
    }

private:
    double rf_freq;           // Carrier frequency
    int sf;                   // Spreading factor
    double bw;                // Bandwidth
    double fs;                // Sampling frequency
    int cr;                   // Coding rate
    int payload_len;          // Payload length
    bool has_header;          // Explicit header
    bool crc;                 // CRC Check enabled
    bool ldr;                 // Low Data Rate Optimization
    int preamble_len;         // Preamble length
    std::vector<uint8_t> whitening_seq;
    int fft_len;
    int sample_num;
    std::vector<Sample> downchirp;
    std::vector<Sample> upchirp;

    // Method to initialize the whitening sequence
    void initWhiteningSequence() {
        whitening_seq = {
            0xff, 0xfe, 0xfc, 0xf8, 0xf0, 0xe1, 0xc2, 0x85, 0x0b, 0x17, 0x2f, 0x5e,
            0xbc, 0x78, 0xf1, 0xe3, 0xc6, 0x8d, 0x1a, 0x34, 0x68, 0xd0, 0xa0, 0x40,
            0x80, 0x01, 0x02, 0x04, 0x08, 0x11, 0x23, 0x47, 0x8e, 0x1c, 0x38, 0x71,
            0xe2, 0xc4, 0x89, 0x12, 0x25, 0x4b, 0x97, 0x2e, 0x5c, 0xb8, 0x70, 0xe0,
            0xc0, 0x81, 0x03, 0x06, 0x0c, 0x19, 0x32, 0x64, 0xc9, 0x92, 0x24, 0x49
        };
    }

    // Method to encode a single byte
    int encodeByte(int byte) {
        int fec_encoded = applyFEC(byte);
        int interleaved = interleaveBits(fec_encoded);
        int whitened = interleaved ^ whitening_seq[interleaved % whitening_seq.size()];
        return whitened;
    }

    // Method to apply forward error correction (FEC)
    int applyFEC(int byte) {
        std::bitset<8> data_bits(byte);
        std::bitset<13> encoded_bits;

        for (int i = 0; i < 4; ++i) {
            encoded_bits[i] = data_bits[i];
        }

        for (int i = 4; i < 4 + cr; ++i) {
            encoded_bits[i] = encoded_bits[0] ^ encoded_bits[1] ^ encoded_bits[2] ^ encoded_bits[3];
        }

        return static_cast<int>(encoded_bits.to_ulong());
    }

    // Method to interleave bits based on spreading factor and coding rate
    int interleaveBits(int data) {
        std::bitset<13> bits(data);
        std::bitset<13> interleaved;

        for (int i = 0; i < 13; ++i) {
            int new_pos = (i * sf + cr) % 13;
            interleaved[new_pos] = bits[i];
        }

        return static_cast<int>(interleaved.to_ulong());
    }

    // Gray code encoding
    int grayEncode(int binary) {
        return binary ^ (binary >> 1);
    }

    // Gray decoding function
    int grayDecode(int gray) {
        int binary = gray;
        for (int shift = 1; shift < 8; shift <<= 1) {
            binary ^= (gray >> shift);
        }
        return binary;
    }
    
    void initializeChirps() {
        double symbol_duration = static_cast<double>(sample_num) / fs; // Duration of one symbol
        double chirp_rate = bw / symbol_duration; // Chirp rate (BW/T)

        // Generate the upchirp
        for (int i = 0; i < sample_num; ++i) {
            double t = static_cast<double>(i) / fs; // Time for the current sample
            double phase = 2 * M_PI * (0.5 * chirp_rate * t * t);
            upchirp[i] = {static_cast<float>(std::cos(phase)), static_cast<float>(std::sin(phase))};
        }

        // Generate the downchirp
        for (int i = 0; i < sample_num; ++i) {
            double t = static_cast<double>(i) / fs; // Time for the current sample
            double phase = 2 * M_PI * (0.5 * chirp_rate * t * t);
            downchirp[i] = {static_cast<float>(std::cos(-phase)), static_cast<float>(std::sin(-phase))};
        }
    }

    // Method to generate a chirp for a given symbol
    std::vector<Sample> generateChirp(int samples_per_symbol, bool is_up) {
        std::vector<Sample> chirp(samples_per_symbol);

        // Number of points in a full cycle (2^SF)
        int N = std::pow(2, sf);

        // Symbol duration in seconds and frequency slope
        double T = static_cast<double>(N) / bw;  // Symbol period
        double k = is_up ? bw / T : -bw / T;     // Frequency slope (up or down chirp)
        double f0 = is_up ? -bw / 2 : bw / 2;    // Initial frequency (up or down chirp)

        // Time vector (normalized to sampling rate)
        for (int i = 0; i < samples_per_symbol; ++i) {
            double t = static_cast<double>(i) / fs;  // Time for the current sample

            // Compute the instantaneous phase: phase(t) = 2π(f0 * t + 0.5 * k * t^2)
            double phase = 2 * M_PI * (f0 * t + 0.5 * k * t * t);

            // Calculate the I and Q components based on the phase
            chirp[i].I = std::cos(phase);
            chirp[i].Q = std::sin(phase);
        }

        return chirp;
    }
    
    int bitReduce(const std::vector<int>& codewords, const std::vector<int>& indices) {
        int result = 0;
        for (int index : indices) {
            result ^= codewords[index - 1]; // MATLAB uses 1-based indexing; adjust for 0-based in C++
        }
        return result;
    }

    float magnitude(const Sample& sample) {
        return std::sqrt(sample.I * sample.I + sample.Q * sample.Q);
    }

    std::vector<int> intToBinary(int num, int length) {
        std::vector<int> binary(length, 0);
        for (int i = 0; i < length; ++i) {
            binary[length - 1 - i] = (num >> i) & 1;
        }
        return binary;
    }

    uint16_t calculateCRC(const std::vector<int>& data) {
        uint16_t crc = 0xFFFF;
        const uint16_t polynomial = 0x1021; // x^16 + x^12 + x^5 + 1

        for (auto byte : data) {
            crc ^= (byte << 8);
            for (int i = 0; i < 8; ++i) {
                if (crc & 0x8000) {
                    crc = (crc << 1) ^ polynomial;
                } else {
                    crc <<= 1;
                }
            }
        }
        return crc;
    }

    bool validateChecksum(const std::vector<int>& data) {
        // Extract the payload and the CRC from the data
        std::vector<int> payload(data.begin(), data.end() - 2);
        uint16_t received_crc = (data[data.size() - 2] << 8) | data[data.size() - 1];

        // Calculate the CRC for the payload
        uint16_t calculated_crc = calculateCRC(payload);

        // Return true if the CRCs match, false otherwise
        return received_crc == calculated_crc;
    }
    
    int nextPowerOfTwo(int n) {
        return pow(2, ceil(log2(n)));
    }

};

#endif /* LoRaPHY_h */
