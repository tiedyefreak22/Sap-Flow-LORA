#ifndef LORAPHY_H
#define LORAPHY_H

#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include <fftw3.h>
#include <iostream>
#include <bitset>

class LoRaPHY {
public:
    LoRaPHY(double rf_freq, int sf, double bw, double fs, uint8_t sync_word = 0x34)
    : rf_freq(rf_freq), sf(sf), bw(bw), fs(fs), sync_word(sync_word), cfo(0.0), symbol_timing_offset(0), crc_enabled(true), has_header(true), whitening_sequence(generate_whiten_sequence(255)) {
        initialize_header_checksum_matrix();
    }

    // Modulation
    Sample* modulate(const std::vector<int>& symbols, int& length) {
        Sample* preamble_signal = generate_preamble(sync_word);
        int preamble_length = preamble_len * (2 * std::pow(2, sf));

        std::vector<int> header = create_header(symbols.size());
        header = encode(whiten(header));

        Sample* header_signal = modulate_symbols(header, length);
        int data_length;
        Sample* data_signal = modulate_symbols(encode(whiten(symbols)), data_length);

        Sample* full_signal = new Sample[preamble_length + length + data_length];
        std::copy(preamble_signal, preamble_signal + preamble_length, full_signal);
        std::copy(header_signal, header_signal + length, full_signal + preamble_length);
        std::copy(data_signal, data_signal + data_length, full_signal + preamble_length + length);

        length += preamble_length + data_length;
        delete[] preamble_signal;
        delete[] header_signal;
        delete[] data_signal;

        return full_signal;
    }

    // Demodulation
    std::pair<std::vector<int>, bool> demodulate(const Sample* signal, int length) {
        auto [header_symbols, data_start] = demodulate_symbols(signal, 8, true);
        std::cout << "Demodulated Header Symbols: ";
        for (int val : header_symbols) std::cout << val << " ";
        std::cout << std::endl;

        auto dewhitened_header = dewhiten(header_symbols);
        std::cout << "Dewhitened Header: ";
        for (int val : dewhitened_header) std::cout << val << " ";
        std::cout << std::endl;

        auto [decoded_header, valid_header] = decode(dewhitened_header);

        std::cout << "Decoded Header: ";
        for (int val : decoded_header) std::cout << val << " ";
        std::cout << std::endl;
        if (!valid_header) {
            std::cerr << "Header is invalid, transmission may be corrupted." << std::endl;
            return {{}, false};
        }

        int payload_length = decoded_header[1];
        auto [data_symbols, _] = demodulate_symbols(signal, payload_length, false, data_start);
        std::cout << "Demodulated Data Symbols: ";
        for (int val : data_symbols) std::cout << val << " ";
        std::cout << std::endl;

        auto dewhitened_data = dewhiten(data_symbols);
        auto [decoded_data, valid_crc] = decode(dewhitened_data);

        std::cout << "Decoded Data: ";
        for (int val : decoded_data) std::cout << val << " ";
        std::cout << std::endl;
        if (!valid_crc) {
            std::cerr << "CRC check failed, data may be corrupted." << std::endl;
            return {{}, false};
        }

        return {decoded_data, true};
    }

    // Encoding
    std::vector<int> encode(const std::vector<int>& data) {
        auto encoded_data = hamming_encode(data);
        encoded_data = gray_code(encoded_data);
        encoded_data = diagonal_interleave(encoded_data);
        return encoded_data;
    }
    
    // Decoding
    std::pair<std::vector<int>, bool> decode(const std::vector<int>& symbols) {
        auto deinterleaved_symbols = diagonal_deinterleave(symbols);
        auto gray_decoded_symbols = gray_decode(deinterleaved_symbols);
        auto hamming_decoded_data = hamming_decode(gray_decoded_symbols);
        bool checksum_valid = crc_enabled ? compute_crc(hamming_decoded_data) : true;
        return {hamming_decoded_data, checksum_valid};
    }

    // Whitening function
    std::vector<int> whiten(const std::vector<int>& data) {
        std::vector<int> whitened_data(data.size());
        for (size_t i = 0; i < data.size(); ++i) {
            whitened_data[i] = data[i] ^ whitening_sequence[i % whitening_sequence.size()];
        }
        return whitened_data;
    }

    // Dewhitening function
    std::vector<int> dewhiten(const std::vector<int>& data) {
        return whiten(data);
    }

private:
    double rf_freq;
    int sf;
    double bw;
    double fs;
    uint8_t sync_word;
    int preamble_len = 6;
    double cfo;
    int symbol_timing_offset;
    bool crc_enabled;
    bool has_header;
    std::vector<int> whitening_sequence;
    int header_checksum_matrix[5][12];

    // Initialize header checksum matrix
    void initialize_header_checksum_matrix() {
        int matrix[5][12] = {
            {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
            {1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1},
            {0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0},
            {0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1},
            {0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1}
        };
        std::copy(&matrix[0][0], &matrix[0][0] + 5 * 12, &header_checksum_matrix[0][0]);
    }

    // Create header
    std::vector<int> create_header(int payload_length) {
        std::vector<int> header(5, 0);
        header[0] = 0;
        header[1] = payload_length;
        header[2] = 3;
        header[3] = 0;

        std::vector<int> data_bits;
        for (int i = 0; i < 3; ++i) {
            for (int j = 3; j >= 0; --j) {
                data_bits.push_back((header[i] >> j) & 1);
            }
        }
        std::bitset<12> checksum = calculate_matrix_checksum(data_bits);
        header[4] = static_cast<int>(checksum.to_ulong());

        std::cout << "Created Header: ";
        for (int val : header) std::cout << val << " ";
        std::cout << std::endl;

        return header;
    }

    // Calculate checksum using matrix
    std::bitset<12> calculate_matrix_checksum(const std::vector<int>& data_bits) {
        std::bitset<12> checksum;
        for (int col = 0; col < 12; ++col) {
            int sum = 0;
            for (int row = 0; row < 5; ++row) {
                sum += data_bits[row] * header_checksum_matrix[row][col];
            }
            checksum[col] = sum % 2;
        }
        return checksum;
    }

    // Generate whitening sequence
    std::vector<int> generate_whiten_sequence(int length) {
        std::vector<int> sequence(length);
        uint8_t reg = 0xFF;
        for (int i = 0; i < length; ++i) {
            sequence[i] = reg;
            reg = (reg << 1) ^ ((reg & 0x80) ? 0xE1 : 0);
        }
        return sequence;
    }

    // Generate preamble
    Sample* generate_preamble(uint8_t sync_word) {
        int samples_per_symbol = 2 * std::pow(2, sf);
        Sample* preamble_signal = new Sample[samples_per_symbol * preamble_len];
        Sample* sync_chirp = chirp(true, sf, bw, fs, sync_word, 0, samples_per_symbol);

        for (int i = 0; i < preamble_len - 1; ++i) {
            Sample* upchirp = chirp(true, sf, bw, fs, 0, 0, samples_per_symbol);
            std::copy(upchirp, upchirp + samples_per_symbol, preamble_signal + i * samples_per_symbol);
            delete[] upchirp;
        }
        std::copy(sync_chirp, sync_chirp + samples_per_symbol, preamble_signal + (preamble_len - 1) * samples_per_symbol);

        delete[] sync_chirp;
        return preamble_signal;
    }

    // Modulate symbols into chirp signals
    Sample* modulate_symbols(const std::vector<int>& symbols, int& length) {
        int samples_per_symbol = 2 * std::pow(2, sf);
        length = symbols.size() * samples_per_symbol;
        Sample* result = new Sample[length];

        for (size_t i = 0; i < symbols.size(); ++i) {
            Sample* chirp_signal = chirp(true, sf, bw, fs, symbols[i], cfo, samples_per_symbol);
            std::copy(chirp_signal, chirp_signal + samples_per_symbol, result + i * samples_per_symbol);
            delete[] chirp_signal;
        }
        return result;
    }

    // Estimate CFO
    double estimate_cfo(const Sample* signal, int length) {
        int samples_per_symbol = 2 * std::pow(2, sf);
        int preamble_samples = preamble_len * samples_per_symbol;

        std::vector<std::complex<float>> mixed_signal(preamble_samples);
        Sample* ref_chirp = chirp(true, sf, bw, fs, 0, 0, samples_per_symbol);

        for (int i = 0; i < preamble_samples; ++i) {
            float I = signal[i].I * ref_chirp[i % samples_per_symbol].I + signal[i].Q * ref_chirp[i % samples_per_symbol].Q;
            float Q = signal[i].Q * ref_chirp[i % samples_per_symbol].I - signal[i].I * ref_chirp[i % samples_per_symbol].Q;
            mixed_signal[i] = std::complex<float>(I, Q);
        }
        delete[] ref_chirp;

        fftwf_complex* in = reinterpret_cast<fftwf_complex*>(mixed_signal.data());
        fftwf_complex* out = fftwf_alloc_complex(preamble_samples);
        fftwf_plan plan = fftwf_plan_dft_1d(preamble_samples, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

        fftwf_execute(plan);

        int peak_index = 0;
        float max_magnitude = 0;
        for (int i = 0; i < preamble_samples; ++i) {
            float magnitude = std::sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]);
            if (magnitude > max_magnitude) {
                max_magnitude = magnitude;
                peak_index = i;
            }
        }

        fftwf_destroy_plan(plan);
        fftwf_free(out);

        double frequency_offset = (peak_index - preamble_samples / 2) * (fs / preamble_samples);
        return frequency_offset;
    }

    // Recover symbol timing
    int recover_symbol_timing(const Sample* signal, int length) {
        int samples_per_symbol = 2 * std::pow(2, sf);
        int preamble_samples = preamble_len * samples_per_symbol;

        std::vector<std::complex<float>> mixed_signal(preamble_samples);
        Sample* ref_chirp = chirp(true, sf, bw, fs, sync_word, 0, samples_per_symbol);

        for (int i = 0; i < preamble_samples; ++i) {
            float I = signal[i].I * ref_chirp[i % samples_per_symbol].I + signal[i].Q * ref_chirp[i % samples_per_symbol].Q;
            float Q = signal[i].Q * ref_chirp[i % samples_per_symbol].I - signal[i].I * ref_chirp[i % samples_per_symbol].Q;
            mixed_signal[i] = std::complex<float>(I, Q);
        }
        delete[] ref_chirp;

        fftwf_complex* in = reinterpret_cast<fftwf_complex*>(mixed_signal.data());
        fftwf_complex* out = fftwf_alloc_complex(preamble_samples);
        fftwf_plan plan = fftwf_plan_dft_1d(preamble_samples, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

        fftwf_execute(plan);

        int peak_index = 0;
        float max_magnitude = 0;
        for (int i = 0; i < preamble_samples; ++i) {
            float magnitude = std::sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]);
            if (magnitude > max_magnitude) {
                max_magnitude = magnitude;
                peak_index = i;
            }
        }

        fftwf_destroy_plan(plan);
        fftwf_free(out);

        return peak_index;
    }

    // Demodulate symbols
    std::pair<std::vector<int>, int> demodulate_symbols(const Sample* signal, int symbol_count, bool apply_cfo, int start = 0) {
        int samples_per_symbol = 2 * std::pow(2, sf);
        std::vector<int> symbols;

        for (int i = 0; i < symbol_count; ++i) {
            int symbol = fft_dechirp(signal, start + i * samples_per_symbol, apply_cfo);
            symbols.push_back(symbol);
        }
        return {symbols, start + symbol_count * samples_per_symbol};
    }

    int fft_dechirp(const Sample* signal, int start_idx, bool apply_cfo) {
        int samples_per_symbol = 2 * std::pow(2, sf);
        std::vector<std::complex<float>> mixed_signal(samples_per_symbol);

        Sample* ref_chirp = chirp(true, sf, bw, fs, 0, apply_cfo ? cfo : 0, samples_per_symbol);

        for (int i = 0; i < samples_per_symbol; ++i) {
            float I = signal[start_idx + i].I * ref_chirp[i].I + signal[start_idx + i].Q * ref_chirp[i].Q;
            float Q = signal[start_idx + i].Q * ref_chirp[i].I - signal[start_idx + i].I * ref_chirp[i].Q;
            mixed_signal[i] = std::complex<float>(I, Q);
        }
        delete[] ref_chirp;

        fftwf_complex* in = reinterpret_cast<fftwf_complex*>(mixed_signal.data());
        fftwf_complex* out = fftwf_alloc_complex(samples_per_symbol);
        fftwf_plan plan = fftwf_plan_dft_1d(samples_per_symbol, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

        fftwf_execute(plan);

        float max_magnitude = 0;
        int peak_index = 0;
        for (int i = 0; i < samples_per_symbol; ++i) {
            float magnitude = std::sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]);
            if (magnitude > max_magnitude) {
                max_magnitude = magnitude;
                peak_index = i;
            }
        }

        fftwf_destroy_plan(plan);
        fftwf_free(out);

        return peak_index;
    }

    Sample* chirp(bool is_up, int sf, double bw, double fs, double freq_offset, double cfo, int samples) {
        int N = (samples > 0) ? samples : static_cast<int>(2 * std::pow(2, sf));
        Sample* chirp_signal = new Sample[N];
        double T = N / fs;
        double k = bw / T;
        double sign = is_up ? 1 : -1;

        for (int n = 0; n < N; ++n) {
            double t = n / fs;
            double phase = 2 * M_PI * (freq_offset * t + sign * 0.5 * k * t * t) + cfo * t;
            chirp_signal[n].I = cos(phase);
            chirp_signal[n].Q = sin(phase);
        }
        return chirp_signal;
    }

    // Interleaving and encoding functions
    std::vector<int> diagonal_interleave(const std::vector<int>& data) {
        int M = sf;
        int N = (data.size() / M) * M;
        std::vector<int> interleaved_data(N);

        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < M; ++j) {
                int idx = j * M + ((i + j) % M);
                if (idx < N) {
                    interleaved_data[idx] = data[i * M + j];
                }
            }
        }
        return interleaved_data;
    }

    std::vector<int> diagonal_deinterleave(const std::vector<int>& data) {
        int M = sf;
        int N = (data.size() / M) * M;
        std::vector<int> deinterleaved_data(N);

        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < M; ++j) {
                int idx = j * M + ((i + j) % M);
                if (idx < N) {
                    deinterleaved_data[i * M + j] = data[idx];
                }
            }
        }
        return deinterleaved_data;
    }

    std::vector<int> gray_code(const std::vector<int>& data) {
        std::vector<int> gray_coded_data;
        for (int value : data) {
            gray_coded_data.push_back(value ^ (value >> 1));
        }
        return gray_coded_data;
    }

    std::vector<int> gray_decode(const std::vector<int>& data) {
        std::vector<int> gray_decoded_data;
        for (int value : data) {
            int decoded = value;
            int mask = value >> 1;
            while (mask != 0) {
                decoded ^= mask;
                mask >>= 1;
            }
            gray_decoded_data.push_back(decoded);
        }
        return gray_decoded_data;
    }

    std::vector<int> hamming_encode(const std::vector<int>& data) {
        std::vector<int> encoded_data;
        for (int byte : data) {
            int d1 = (byte >> 3) & 1;
            int d2 = (byte >> 2) & 1;
            int d3 = (byte >> 1) & 1;
            int d4 = byte & 1;

            int p1 = d1 ^ d2 ^ d4;
            int p2 = d1 ^ d3 ^ d4;
            int p3 = d2 ^ d3 ^ d4;

            int encoded_byte = (p1 << 6) | (p2 << 5) | (d1 << 4) | (p3 << 3) | (d2 << 2) | (d3 << 1) | d4;
            encoded_data.push_back(encoded_byte);
        }
        return encoded_data;
    }

    std::vector<int> hamming_decode(const std::vector<int>& data) {
        std::vector<int> decoded_data;
        for (int encoded_byte : data) {
            int p1 = (encoded_byte >> 6) & 1;
            int p2 = (encoded_byte >> 5) & 1;
            int d1 = (encoded_byte >> 4) & 1;
            int p3 = (encoded_byte >> 3) & 1;
            int d2 = (encoded_byte >> 2) & 1;
            int d3 = (encoded_byte >> 1) & 1;
            int d4 = encoded_byte & 1;

            int s1 = p1 ^ d1 ^ d2 ^ d4;
            int s2 = p2 ^ d1 ^ d3 ^ d4;
            int s3 = p3 ^ d2 ^ d3 ^ d4;
            int syndrome = (s1 << 2) | (s2 << 1) | s3;

            if (syndrome != 0) {
                encoded_byte ^= (1 << (7 - syndrome));
            }

            int decoded_byte = (d1 << 3) | (d2 << 2) | (d3 << 1) | d4;
            decoded_data.push_back(decoded_byte);
        }
        return decoded_data;
    }

    bool compute_crc(const std::vector<int>& data) {
        uint16_t crc = 0xFFFF;
        for (int byte : data) {
            crc ^= static_cast<uint8_t>(byte) << 8;
            for (int i = 0; i < 8; ++i) {
                crc = (crc & 0x8000) ? (crc << 1) ^ 0x1021 : crc << 1;
            }
        }
        return crc == 0;
    }
};

#endif // LORAPHY_H
