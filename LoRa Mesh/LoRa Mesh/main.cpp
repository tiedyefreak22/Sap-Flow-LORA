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
        std::vector<Sample> samples;
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
            
            std::vector<int> symbols = phy.encode(message);

            // Baseband Modulation
            samples = phy.modulate(symbols);
                        
            // Print the first few samples
            for (int i = 0; i < 10; ++i) {
                printf("Sample %d: I = %f, Q = %f\n", i, samples[i].I, samples[i].Q);
            }
            
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
        // Demodulation
        int num_symbols = 13; // Adjust based on the message length and encoding

        // Estimate CFO
        double cfo = phy.estimateCFO(samples, M, 8, Fs, BW);

        // Demodulate the signal
        std::vector<int> demodulated_symbols = phy.demodulate(samples, M, num_symbols, Fs, BW, cfo);

        for (int i = 0; i < sizeof(demodulated_symbols); i++) {
            printf("%d\n", demodulated_symbols[i]);
        }
        // Decode the demodulated symbols
        std::vector<uint8_t> decoded_data = phy.decode(demodulated_symbols);

        // Convert decoded data to a character array
        char* decoded_message = phy.convertToCharArray(decoded_data);

        // Print the decoded message
        std::cout << "Decoded message: " << decoded_message << std::endl;

        // Free dynamically allocated memory
        delete[] decoded_message;
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
