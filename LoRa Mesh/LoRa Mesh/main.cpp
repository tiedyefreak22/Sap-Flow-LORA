//
//  main.cpp
//  LoRa Mesh
//
//  Created by Kevin Hardin on 9/30/24.
//

#include <iostream>
#include "RTLSDR.h"
#include "AlignTrack.h"
#include "LoRaPHY.h"
#include <cstdlib>
#include <vector>
#include <complex>
#include <ctime>
#include <cstring> // For strlen

int load_data = 0;
int SDR = 0;
int sim_data = 1;

// LoRa Parameters
int SF = 8; // Spreading factor
int M=2^SF; // no. of samples in one symbol
int BW = 125e3 ; // Bandwidth
int Fs = 1e6;  // sampling freq
int fc = 915e6 ;  // carrier center frequency
float noise_sigma = 0.0;
const char* message = "Hello world!";
int numPreambleSymbols = 8; // Standard LoRa preamble symbol count
double symbolRate = BW / (2^SF); // Symbol rate
int shift = 0;

// LoRa Setup
LoRaPHY phy(fc, SF, BW, Fs);

int main(int argc, const char * argv[]) {
    try {
        //std::vector<Sample> samples;
        Sample* samples;
        int modulated_length;
        
        if (load_data) {
            // Load data from a .bin file
            const char* filename = "/Users/kevinhardin/Documents/GitHub/Becca-Sap-Flow-LORA/output.bin";
            size_t sampleCount = 0;

            samples = readBinaryFile(filename, &sampleCount);
//            if (samples == NULL) {
//                return 1; // Error reading file
//            }

            // Print the first few samples
            for (size_t i = 0; i < sampleCount && i < 10; ++i) {
                printf("Sample %zu: I = %f, Q = %f\n", i, samples[i].I, samples[i].Q);
            }
        } else if (sim_data) {
            // Encode payload
            printf("[encode] message: %s\n", message);
            
            std::vector<int> preencoded_data;

            // Convert each character to an integer and apply encoding steps
            for (size_t i = 0; i < sizeof(message); ++i) {
                preencoded_data.push_back(static_cast<int>(message[i])); // Convert char to int
            }

            // Directly encode, modulate, demodulate, and decode without transformations
            std::vector<int> encoded_data = phy.encode(preencoded_data);
            int modulated_length;
            samples = phy.modulate(encoded_data, modulated_length);
            
//            for (size_t i = 0; i < sizeof(samples); i++) {
//                printf("%f, %f\n", samples[i].I, samples[i].Q);
//            }
            
        } else if (SDR) {
            // Capture data from RTL-SDR
            init_rtl_sdr();
            samples = receive_rtl_sdr(dev);
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//        // Convert raw data to Signal format and take FFT
//        [received_fft] = fft(received_signal);
//        new_received_fft=received_fft;
//        
//        new_received_signal = ifft(new_received_fft(:,1));
//
        // Step 3: Demodulate the received signal
        // Sample test without whitening and error correction
        // Demodulate the received signal
        auto [demodulated_symbols, header_valid] = phy.demodulate(samples, modulated_length);

        if (!header_valid) {
            std::cerr << "Header is invalid, transmission may be corrupted." << std::endl;
        } else {
            auto [decoded_data, crc_valid] = phy.decode(demodulated_symbols);
            // Print received data if CRC passes
            if (crc_valid) {
                for (int byte : decoded_data) {
                    std::cout << static_cast<char>(byte);
                }
            } else {
                std::cerr << "CRC check failed, data may be corrupted." << std::endl;
            }
        }

        // Clean up dynamically allocated memory
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // std::vector<Signal> overlappingSignals;
        // Fill overlappingSignals with processed data
        
        // Apply AlignTrack deconfliction algorithm
        // std::vector<Signal> deconflictedSignals = alignTrack(overlappingSignals);

        // Decode each deconflicted signal
//        for (const auto& signal : deconflictedSignals) {
//            decodeLoRaSymbols(signal);
//        }

    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
        return 1;
    }

    return 0;

    
}
