#ifndef LoRaPHY_h
#define LoRaPHY_h

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <bitset>

class LoRaPHY {
public:
    // Constructor
    LoRaPHY(double rf, int spreading_factor, double bandwidth, double sampling_freq)
        : rf_freq(rf), sf(spreading_factor), bw(bandwidth), fs(sampling_freq),
          cr(1), payload_len(0), has_header(true), crc(true), ldr(false), preamble_len(8) {
        initWhiteningSequence();
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
        std::vector<Sample> modulated_signal;
        int samples_per_symbol = std::pow(2, sf);

        for (int symbol : symbols) {
            // Generate chirp for this symbol
            auto chirp = generateChirp(samples_per_symbol, 1);
            modulated_signal.insert(modulated_signal.end(), chirp.begin(), chirp.end());
        }

        return modulated_signal;
    }

    // Function to generate an upchirp and downchirp for alignment and demodulation
    std::vector<Sample> generateUpChirp(int samples_per_symbol) {
        std::vector<Sample> upchirp = generateChirp(samples_per_symbol, 1);
        return upchirp;
    };
    std::vector<Sample> generateDownChirp(int samples_per_symbol) {
        std::vector<Sample> downchirp = generateChirp(samples_per_symbol, 0);
        return downchirp;
    };

    // Function to align the signal based on the preamble
    int alignSignal(const std::vector<Sample>& received_signal, int samples_per_symbol) {
        int signal_length = received_signal.size();
        int best_index = -1;
        double max_correlation = 0.0;

        auto upchirp = generateUpChirp(samples_per_symbol);

        // Search for the highest correlation with the upchirp to find the preamble
        for (int i = 0; i <= signal_length - samples_per_symbol; ++i) {
            std::complex<double> correlation_sum(0.0, 0.0);

            for (int j = 0; j < samples_per_symbol; ++j) {
                std::complex<double> received(received_signal[i + j].I, received_signal[i + j].Q);
                std::complex<double> up(upchirp[j].I, upchirp[j].Q);

                correlation_sum += received * std::conj(up);
            }

            double magnitude = std::abs(correlation_sum);
            if (magnitude > max_correlation) {
                max_correlation = magnitude;
                best_index = i;
            }
        }

        if (best_index == -1) {
            std::cerr << "Error: Preamble not found in the signal." << std::endl;
            return -1;
        }

        std::cout << "Preamble found at index: " << best_index << std::endl;
        return best_index;
    }

    // Function to demodulate the signal using downchirps and compensate for CFO
    std::vector<int> demodulate(const std::vector<Sample>& received_signal, int samples_per_symbol, int num_symbols, double fs, double bw, double cfo) {
        std::vector<int> symbols;
        int preamble_symbols = 8; // Number of preamble symbols

        // Align the signal to the start of the preamble
        int start_index = alignSignal(received_signal, samples_per_symbol);
        if (start_index == -1) {
            std::cerr << "Error: Could not align the signal for demodulation." << std::endl;
            return symbols; // Return empty if alignment fails
        }

        // Iterate over the number of expected symbols after alignment
        for (int sym_idx = 0; sym_idx < preamble_symbols + num_symbols; ++sym_idx) {
            int symbol_start = start_index + sym_idx * samples_per_symbol;

            // Ensure we don't go out of bounds when extracting the symbol segment
            if (symbol_start + samples_per_symbol > received_signal.size()) {
                std::cerr << "Error: Not enough samples in the signal to demodulate all symbols." << std::endl;
                break;
            }

            // Extract the segment of the received signal corresponding to one symbol
            std::vector<Sample> symbol_segment(
                received_signal.begin() + symbol_start,
                received_signal.begin() + symbol_start + samples_per_symbol
            );

            auto downchirp = generateDownChirp(samples_per_symbol);

            std::complex<double> accumulated_phase(0.0, 0.0);
            for (int i = 0; i < samples_per_symbol; ++i) {
                std::complex<double> received(symbol_segment[i].I, symbol_segment[i].Q);
                std::complex<double> down(downchirp[i].I, downchirp[i].Q);

                // Apply CFO correction: phase shift based on the estimated CFO
                double phase_shift = -2 * M_PI * cfo * i / fs; // Negative sign to compensate for CFO
                std::complex<double> correction(std::cos(phase_shift), std::sin(phase_shift));
                received *= correction;

                accumulated_phase += received * std::conj(down);
            }

            int symbol_value = std::round(std::arg(accumulated_phase) * std::pow(2, sf) / (2 * M_PI));
            symbol_value = (symbol_value + static_cast<int>(std::pow(2, sf))) % static_cast<int>(std::pow(2, sf));

            symbols.push_back(symbol_value);
        }

        return symbols;
    }


    // Estimate CFO using the preamble
    double estimateCFO(const std::vector<Sample>& received_signal, int samples_per_symbol, int preamble_symbols, double fs, double bw) {
        std::complex<double> accumulated_phase(0.0, 0.0);

        // Generate the downchirp used for correlation
        std::vector<Sample> downchirp = generateDownChirp(samples_per_symbol);

        // Iterate over each preamble symbol and calculate the correlation
        for (int sym_idx = 0; sym_idx < preamble_symbols; ++sym_idx) {
            std::complex<double> phase_sum(0.0, 0.0);

            for (int i = 0; i < samples_per_symbol; ++i) {
                std::complex<double> received(received_signal[sym_idx * samples_per_symbol + i].I,
                                              received_signal[sym_idx * samples_per_symbol + i].Q);
                std::complex<double> down(downchirp[i].I, downchirp[i].Q);

                // Correlate the received signal with the downchirp
                phase_sum += received * std::conj(down);
            }

            // Debug: Output the phase angle for each preamble symbol
            double phase_angle = std::arg(phase_sum);
            std::cout << "Symbol " << sym_idx << " phase angle: " << phase_angle << std::endl;

            // Accumulate the overall phase
            accumulated_phase += phase_sum;
        }

        // Calculate the average phase shift per sample
        double avg_phase_angle = std::arg(accumulated_phase) / preamble_symbols;
        double cfo_estimate = avg_phase_angle / (2 * M_PI * samples_per_symbol / fs);

        // Threshold to check if CFO is effectively zero
        double cfo_threshold = 1e-6; // Adjust the threshold based on precision
        if (std::abs(cfo_estimate) < cfo_threshold) {
            cfo_estimate = 0.0; // Set CFO to zero if below the threshold
        }

        std::cout << "[CFO Estimate] Estimated CFO: " << cfo_estimate << " Hz" << std::endl;
        return cfo_estimate;
    }

    // Hamming (7, 4) decoding based on CR = 1
    uint8_t applyHamming74Decode(int encoded_value) {
        std::bitset<7> encoded_bits(encoded_value);

        int parity_bit1 = encoded_bits[4];
        int parity_bit2 = encoded_bits[5];
        int parity_bit3 = encoded_bits[6];

        int syndrome1 = encoded_bits[0] ^ encoded_bits[1] ^ encoded_bits[3] ^ parity_bit1;
        int syndrome2 = encoded_bits[0] ^ encoded_bits[2] ^ encoded_bits[3] ^ parity_bit2;
        int syndrome3 = encoded_bits[1] ^ encoded_bits[2] ^ encoded_bits[3] ^ parity_bit3;

        int syndrome = (syndrome1 << 2) | (syndrome2 << 1) | syndrome3;

        if (syndrome != 0 && syndrome <= 7) {
            encoded_bits.flip(syndrome - 1);
        }

        uint8_t corrected_data = (encoded_bits[0] << 0) | (encoded_bits[1] << 1) |
                                 (encoded_bits[2] << 2) | (encoded_bits[3] << 3);

        return corrected_data;
    }

    // Function to decode symbols using Gray decoding and FEC
    std::vector<uint8_t> decode(const std::vector<int>& symbols) {
        std::vector<uint8_t> decoded_data;

        for (const int& symbol : symbols) {
            int gray_decoded_symbol = grayDecode(symbol);
            uint8_t decoded_byte = applyHamming74Decode(gray_decoded_symbol);
            decoded_data.push_back(decoded_byte);
        }

        return decoded_data;
    }

    // Function to convert the decoded symbols back into a character array
    char* convertToCharArray(const std::vector<uint8_t>& decoded_data) {
        char* message = new char[decoded_data.size() + 1];
        for (size_t i = 0; i < decoded_data.size(); ++i) {
            message[i] = static_cast<char>(decoded_data[i]);
        }
        message[decoded_data.size()] = '\0';
        return message;
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

};

#endif /* LoRaPHY_h */
