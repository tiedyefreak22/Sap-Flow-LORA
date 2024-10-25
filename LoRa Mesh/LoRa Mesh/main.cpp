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
        Sample* samples = nullptr;
        if (load_data) {
            // Load data from a .bin file
            const char* filename = "/Users/kevinhardin/Documents/GitHub/Becca-Sap-Flow-LORA/output.bin";
            size_t sampleCount = 0;

            samples = readBinaryFile(filename, &sampleCount);
            if (samples == NULL) {
                return 1; // Error reading file
            }

            // Print the first few samples
            for (size_t i = 0; i < sampleCount && i < 10; ++i) {
                printf("Sample %zu: I = %f, Q = %f\n", i, samples[i].I, samples[i].Q);
            }
        } else if (sim_data) {
            // Encode payload
            printf("[encode] message: %s\n", message);
            
            // Convert to std::string
            std::string str(message);

            // Create a vector to store the integers
            std::vector<int> intVector;

            // Iterate over each character in the string
            for (char c : str) {
                // Convert each character to an integer and add to the vector
                intVector.push_back(c - '0');
            }
            
            std::vector<int> symbols = phy.encode(intVector);

            // Baseband Modulation
            Sample* signalIQ1 = phy.modulate(symbols);
            std::srand(std::time(nullptr)); // use current time as seed for random generator
            int random_value = (1 + (std::rand() % 10000));
//            // Create a vector with a single complex number
//            std::vector<std::complex<float>> complexVector;
//            complexVector.push_back(std::complex<float>(static_cast<float>(noise), 0.0f));
//    
//            std::vector<std::complex<float>> samples(signalIQ1.size());
            Sample temp[sizeof(signalIQ1)];
            for (int i = 0; i < sizeof(signalIQ1); i++) {
                Sample noise = {static_cast<float>(noise_sigma * ((static_cast<float>(random_value) / 10000.0) * 3 - 3)), 0};
                temp[i] = add(signalIQ1[i], noise);
            }
            
            samples = temp;
            
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
        std::vector<std::complex<float>> new_samples;
        for (int i = 0; i < sizeof(samples); i++){
            new_samples[i] = {samples[i].I, samples[i].Q};
        }
        std::vector<int> symbols_d = phy.demodulate(new_samples);
        for (int i = 0; i < sizeof(symbols_d); i++) {
            printf("%c", symbols_d[i]);
        }
        // Decoding
        auto [data, checksum] = phy.decode(symbols_d);
        char* str = nullptr;
        for (int i = 0; i < sizeof(data) - 2; i++) {
            str[i] = data[i];
            printf("%c", str[i]);
        }
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
