#ifndef LORAPHY_H
#define LORAPHY_H

#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include <fftw3.h>

class LoRaPHY {
public:
    LoRaPHY(double rf_freq, int sf, double bw, double fs, uint8_t sync_word = 0x34)
        : rf_freq(rf_freq), sf(sf), bw(bw), fs(fs), sync_word(sync_word), cfo(0.0), symbol_timing_offset(0), crc_enabled(true), has_header(true), whitening_sequence(generate_whitening_sequence(255)) {}

    // Modulation function with Sample* output
    Sample* modulate(const std::vector<int>& symbols, int& length) {
        // Generate the preamble with sync word at the end
        Sample* preamble_signal = generate_preamble(sync_word);
        int preamble_length = preamble_len * (2 * std::pow(2, sf));

        // Create and encode the header
        std::vector<int> header = create_header(symbols.size());
        header = encode(whiten(header));

        // Modulate the header and data
        Sample* header_signal = modulate_symbols(header, length);
        int data_length;
        Sample* data_signal = modulate_symbols(encode(whiten(symbols)), data_length);

        // Combine preamble, header, and data into one signal
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

    // Demodulation function with sync word handling and timing recovery
    std::pair<std::vector<int>, bool> demodulate(const Sample* signal, int length) {
        cfo = estimate_cfo(signal, length);            // Estimate CFO from preamble
        symbol_timing_offset = recover_symbol_timing(signal, length);  // Correct timing offset

        // Demodulate header and data
        auto [header_symbols, data_start] = demodulate_symbols(signal, 8, true);
        auto [header, valid_header] = decode(dewhiten(header_symbols));

        // In demodulate
        std::cout << "Demodulated Header: ";
        for (int val : header_symbols) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
        
        if (!valid_header) {
            return {{}, false};
        }

        int payload_length = header[0];
        auto [data_symbols, _] = demodulate_symbols(signal, payload_length, false, data_start);
        auto [data, valid_crc] = decode(dewhiten(data_symbols));

        return {data, valid_crc};
    }

    // Encoding function that applies Hamming encoding, Gray coding, and interleaving
    std::vector<int> encode(const std::vector<int>& data) {
        std::vector<int> encoded_data = data;
        encoded_data = hamming_encode(encoded_data);
        encoded_data = gray_code(encoded_data);
        encoded_data = diagonal_interleave(encoded_data);
        return encoded_data;
    }

    // Decoding function that reverses interleaving, Gray coding, and Hamming decoding
    std::pair<std::vector<int>, bool> decode(const std::vector<int>& symbols) {
        std::vector<int> decoded_data = diagonal_deinterleave(symbols);
        decoded_data = gray_code(decoded_data);  // Decode Gray-coded symbols
        decoded_data = hamming_decode(decoded_data);  // Decode Hamming
        bool checksum_valid = crc_enabled ? compute_crc(decoded_data) : true;
        return std::make_pair(decoded_data, checksum_valid);
    }
    
    // Header creation function
    std::vector<int> create_header(int payload_length) {
        std::vector<int> header(8, 0); 
        header[0] = payload_length; 
        header[1] = 0x04;  
        header[2] = crc_enabled ? 1 : 0;
        
        // In create_header
        std::cout << "Created Header: ";
        for (int val : header) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
        return header;
    }

private:
    // Properties
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

    // Generate whitening sequence
    static std::vector<int> generate_whitening_sequence(int length) {
        std::vector<int> sequence(length);
        uint8_t reg = 0xFF;
        for (int i = 0; i < length; ++i) {
            sequence[i] = reg;
            reg = (reg << 1) ^ ((reg & 0x80) ? 0xE1 : 0);
        }
        return sequence;
    }

    // Data whitening
    std::vector<int> whiten(const std::vector<int>& data) {
        std::vector<int> whitened_data(data.size());
        for (size_t i = 0; i < data.size(); ++i) {
            whitened_data[i] = data[i] ^ whitening_sequence[i % whitening_sequence.size()];
        }
        return whitened_data;
    }

    // Data dewhitening
    std::vector<int> dewhiten(const std::vector<int>& data) {
        return whiten(data);
    }

    // Diagonal interleave
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

    // Diagonal deinterleave
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

    // Hamming (7,4) encode
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

    // Hamming (7,4) decode
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

    // Gray coding
    std::vector<int> gray_code(const std::vector<int>& data) {
        std::vector<int> gray_coded_data;
        for (int value : data) {
            gray_coded_data.push_back(value ^ (value >> 1));
        }
        return gray_coded_data;
    }

    // CRC-16 computation (LoRa standard polynomial: x^16 + x^12 + x^5 + 1)
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

    // Generate preamble with sync word
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

    // Symbol timing recovery
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

    // Generate chirp signal
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

    // Modulate symbols into a chirp signal
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

    // Demodulate symbols from a signal
    std::pair<std::vector<int>, int> demodulate_symbols(const Sample* signal, int symbol_count, bool apply_cfo, int start = 0) {
        int samples_per_symbol = 2 * std::pow(2, sf);
        std::vector<int> symbols;

        for (int i = 0; i < symbol_count; ++i) {
            int symbol = fft_dechirp(signal, start + i * samples_per_symbol, apply_cfo);
            symbols.push_back(symbol);
        }
        return {symbols, start + symbol_count * samples_per_symbol};
    }

    // Estimate CFO from preamble
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
};

#endif // LORAPHY_H
