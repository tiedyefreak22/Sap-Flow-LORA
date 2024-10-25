//
//  LoRaPHY.h
//  LoRa Mesh
//
//  Created by Kevin Hardin on 10/23/24.
//

#ifndef LoRaPHY_h
#define LoRaPHY_h

#include <iostream>
#include <vector>
#include <complex>
#include <fstream>
#include <numeric>
#include <cmath>
#include <stdexcept>
#include <bitset>

// LoRaPHY Class Definition
class LoRaPHY {
public:
    // Properties equivalent to MATLAB attributes
    double rf_freq;           // Carrier frequency
    int sf;                   // Spreading factor (e.g., 7, 8, 9, 10, 11, 12)
    double bw;                // Bandwidth (e.g., 125kHz, 250kHz, 500kHz)
    double fs;                // Sampling frequency
    int cr;                   // Code rate: (1 for 4/5, 2 for 4/6, etc.)
    int payload_len;          // Payload length
    bool has_header;          // Explicit header (true) or implicit header (false)
    bool crc;                 // CRC Check enabled (true/false)
    bool ldr;                 // Low Data Rate Optimization enabled (true/false)
    std::vector<uint8_t> whitening_seq;  // Whitening sequence
    std::vector<std::vector<int>> header_checksum_matrix; // Header checksum matrix
    int preamble_len;         // Preamble length
    std::vector<std::complex<float>> sig;  // Signal vector (complex data)

    // Constructor
    LoRaPHY(double rf, int spreading_factor, double bandwidth, double sampling_freq)
        : rf_freq(rf), sf(spreading_factor), bw(bandwidth), fs(sampling_freq),
          cr(1), payload_len(0), has_header(true), crc(true), ldr(false), preamble_len(8) {
        initWhiteningSequence();
        initHeaderChecksumMatrix();
    }

    // Method to initialize whitening sequence
    void initWhiteningSequence() {
        whitening_seq = {
            0xff, 0xfe, 0xfc, 0xf8, 0xf0, 0xe1, 0xc2, 0x85, 0x0b, 0x17, 0x2f, 0x5e,
            0xbc, 0x78, 0xf1, 0xe3, 0xc6, 0x8d, 0x1a, 0x34, 0x68, 0xd0, 0xa0, 0x40,
            0x80, 0x01, 0x02, 0x04, 0x08, 0x11, 0x23, 0x47, 0x8e, 0x1c, 0x38, 0x71,
            0xe2, 0xc4, 0x89, 0x12, 0x25, 0x4b, 0x97, 0x2e, 0x5c, 0xb8, 0x70, 0xe0,
            0xc0, 0x81, 0x03, 0x06, 0x0c, 0x19, 0x32, 0x64, 0xc9, 0x92, 0x24, 0x49
            // Complete the sequence as needed
        };
    }

    // Method to initialize the header checksum matrix
    void initHeaderChecksumMatrix() {
        header_checksum_matrix = {
            {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
            {1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1},
            {0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0},
            {0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1},
            {0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1}
            // Add remaining rows as per MATLAB matrix
        };
    }

    // Method to decode symbols
    std::pair<std::vector<int>, int> decode(const std::vector<int>& symbols) {
        std::vector<int> data;
        int checksum = 0;

        for (const int &symbol : symbols) {
            int decoded_value = decodeSymbol(symbol);
            data.push_back(decoded_value);
        }

        checksum = computeChecksum(data);
        return {data, checksum};
    }

    // Sync method with CFO (Carrier Frequency Offset) elimination
    int sync(int x) {
        bool found = false;
        while (x < sig.size() - sf) {
            auto up_peak = dechirp(x);
            auto down_peak = dechirp(x, false);

            if (std::abs(down_peak.real()) > std::abs(up_peak.real())) {
                found = true;
                break;
            }
            x += sf;
        }

        if (!found) {
            throw std::runtime_error("Sync failed: downchirp not detected.");
        }

        auto pkd = dechirp(x, false);
        int to = (pkd.imag() > bw / 2) ? std::round((pkd.imag() - 1 - bw) / calculateZeroPaddingRatio())
                                        : std::round((pkd.imag() - 1) / calculateZeroPaddingRatio());

        x += to;

        auto pku = dechirp(x - 4 * sf);
        int preamble_bin = std::round(pku.imag());

        float cfo = (preamble_bin > bw / 2) ? (preamble_bin - bw - 1) * bw / sf
                                            : (preamble_bin - 1) * bw / sf;

        adjustForCFO(cfo);
        return x;
    }

    // Method to demodulate the signal
   std::vector<int> demodulate(const std::vector<std::complex<float>>& received_signal) {
        std::vector<int> symbols;
        int index = 0;

        // Iterate over the received signal
        while (index < received_signal.size() - sf) {
            // Calculate the down chirp for the current segment
            std::complex<float> downchirp_result(0.0, 0.0);
            for (int i = 0; i < sf; ++i) {
                // Generate the down chirp signal
                float phase = -2 * M_PI * i * i / sf; // Down chirp (negative frequency ramp)
                std::complex<float> downchirp(std::cos(phase), std::sin(phase));

                // Multiply the received signal by the conjugate of the down chirp
                downchirp_result += received_signal[index + i] * std::conj(downchirp);
            }

            // Calculate the symbol value by converting the dechirped signal
            int symbol = std::round(std::arg(downchirp_result) * sf / (2 * M_PI));

            // Normalize symbol value within expected range
            symbol = (symbol + sf) % sf;

            symbols.push_back(symbol);
            index += sf;
        }
        return symbols;
    }

    // Method to encode data
    std::vector<int> encode(const std::vector<int>& data) {
        std::vector<int> encoded_data;

        for (const int &byte : data) {
            int encoded_byte = encodeByte(byte);
            encoded_data.push_back(encoded_byte);
        }

        return encoded_data;
    }

    // Method to modulate encoded data
    Sample* modulate(const std::vector<int>& symbols) {
        std::vector<std::complex<float>> modulated_signal;

        for (int symbol : symbols) {
            for (int i = 0; i < sf; ++i) {
                float phase = 2 * M_PI * symbol * i / sf;
                std::complex<float> chirp(std::cos(phase), std::sin(phase));
                modulated_signal.push_back(chirp);
            }
        }
        Sample* new_modulated_signal;
        for (int i = 0; i < sizeof(modulated_signal); i++) {
            new_modulated_signal[i].I = modulated_signal[i].real();
            new_modulated_signal[i].Q = modulated_signal[i].imag();
        }

        return new_modulated_signal;
    }

    // File read and write methods remain the same...

private:
    // Helper method to decode a single symbol
    int decodeSymbol(int symbol) {
        int decoded_value = (symbol ^ whitening_seq[symbol % whitening_seq.size()]);
        return decoded_value;
    }

    // Example Hamming (7,4) encoder
    int hamming74Encode(int nibble) {
        // Hamming (7,4) code takes 4 bits and adds 3 parity bits
        int p1 = (nibble >> 0 & 1) ^ (nibble >> 1 & 1) ^ (nibble >> 3 & 1);
        int p2 = (nibble >> 0 & 1) ^ (nibble >> 2 & 1) ^ (nibble >> 3 & 1);
        int p3 = (nibble >> 1 & 1) ^ (nibble >> 2 & 1) ^ (nibble >> 3 & 1);

        int encoded = (nibble & 0x0F) | (p1 << 4) | (p2 << 5) | (p3 << 6);
        return encoded;
    }

    // Adjusted encodeByte method for LoRa encoding
    int encodeByte(int byte) {
        // Step 1: Apply forward error correction (FEC)
        int fec_encoded = applyFEC(byte);

        // Step 2: Interleave bits based on spreading factor and coding rate
        int interleaved = interleaveBits(fec_encoded);

        // Step 3: Apply whitening using the whitening sequence
        int whitened = interleaved ^ whitening_seq[interleaved % whitening_seq.size()];

        // Return the fully encoded byte
        return whitened;
    }

    // Method to apply forward error correction (FEC) based on LoRa coding rate
    int applyFEC(int byte) {
        std::bitset<8> data_bits(byte); // Represent byte as 8 bits
        std::bitset<13> encoded_bits;

        // LoRa uses CR from 4/5 to 4/8. We'll apply additional parity bits based on CR.
        for (int i = 0; i < 4; ++i) {
            encoded_bits[i] = data_bits[i];
        }

        // Calculate parity bits based on coding rate (CR)
        for (int i = 4; i < 4 + cr; ++i) {
            // Simple parity calculation using XOR (expand this for more robust CRC if needed)
            encoded_bits[i] = encoded_bits[0] ^ encoded_bits[1] ^ encoded_bits[2] ^ encoded_bits[3];
        }

        // Combine data and parity bits into the final encoded symbol
        int fec_encoded = static_cast<int>(encoded_bits.to_ulong());
        return fec_encoded;
    }

    // Method to interleave bits based on spreading factor and coding rate
    int interleaveBits(int data) {
        std::bitset<13> bits(data);
        std::bitset<13> interleaved;

        // Interleaving as per LoRa standard: interleave according to SF and CR
        // LoRa defines a specific pattern for how bits should be interleaved.
        // Example pattern for illustration:
        for (int i = 0; i < 13; ++i) {
            int new_pos = (i * sf + cr) % 13; // Simplified formula for demonstration
            interleaved[new_pos] = bits[i];
        }

        return static_cast<int>(interleaved.to_ulong());
    }

    // Method for dechirping the signal
    std::complex<float> dechirp(int x, bool up = true) {
        std::complex<float> result(0.0, 0.0);
        if (x < 0 || x >= sig.size()) {
            throw std::out_of_range("Index out of bounds in dechirp method.");
        }

        for (int i = 0; i < sf; ++i) {
            float phase = calculatePhase(x, i, up);
            std::complex<float> chirp(std::cos(phase), std::sin(phase));
            result += sig[x + i] * std::conj(chirp);
        }
        return result;
    }

    // Adjust for Carrier Frequency Offset (CFO)
    void adjustForCFO(float cfo) {
        for (size_t i = 0; i < sig.size(); ++i) {
            float phase_shift = 2 * M_PI * cfo * i / fs;
            std::complex<float> correction(std::cos(phase_shift), std::sin(phase_shift));
            sig[i] *= correction;
        }
    }

    // Method to compute CRC
    bool computeCRC(const std::vector<int>& data) {
        uint16_t crc = 0xFFFF;

        for (auto byte : data) {
            crc ^= (byte << 8);
            for (int i = 0; i < 8; ++i) {
                if (crc & 0x8000) {
                    crc = (crc << 1) ^ 0x1021;
                } else {
                    crc <<= 1;
                }
            }
        }
        return (crc == 0);
    }

    // Compute checksum for data
    int computeChecksum(const std::vector<int>& data) {
        int checksum = std::accumulate(data.begin(), data.end(), 0) % 256;
        return checksum;
    }

    // Calculate phase for dechirping
    float calculatePhase(int x, int i, bool up) {
        float base_freq = (up ? 1 : -1) * bw / fs;
        return 2 * M_PI * base_freq * (x + i);
    }

    // Calculate the zero-padding ratio
    float calculateZeroPaddingRatio() {
        return (ldr) ? 4.0f / sf : 1.0f;
    }
};


#endif /* LoRaPHY_h */
