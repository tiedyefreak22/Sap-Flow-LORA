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
int sim_data = 1;

// LoRa Parameters
int SF = 8; // Spreading factor
int M=2^SF; // no. of samples in one symbol
int BW = 125e3 ; // Bandwidth
int Fs = 1e6;  // sampling freq
int fc = 915e6 ;  // carrier center frequency
float noise_sigma = 0.0;
std::string message = "Hello world!";
std::vector<int> message_symbols(message.begin(), message.end());
int numPreambleSymbols = 8; // Standard LoRa preamble symbol count
double symbolRate = BW / (2^SF); // Symbol rate
int shift = 0;
uint8_t sync_word = 0x34;          // Default sync word

// LoRa Setup
LoRaPHY phy(fc, SF, BW, Fs);
//phy.setPreambleLength(numPreambleSymbols); // Set preamble length (if applicable)

int main(int argc, const char * argv[]) {
    try {
        //std::vector<Sample> samples;
        Sample* samples;
        int modulated_length = 0;
        
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
            // Modulate the encoded data into a complex sample array
            int modulated_length;
            Sample* modulated_signal = phy.modulate(message_symbols, modulated_length);

            // Debug output to visualize the encoded and modulated signal
            std::cout << "[main] Encoded data: ";
            for (int value : message_symbols) {
                std::cout << value << " ";
            }
            std::cout << std::endl;

            std::cout << "[main] Modulated signal (first few samples): ";
            for (int i = 0; i < std::min(modulated_length, 10); ++i) {
                std::cout << "(" << modulated_signal[i].I << ", " << modulated_signal[i].Q << ") ";
            }
            std::cout << std::endl;
            
        } else if (~sim_data and ~load_data) {
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
        // Demodulate and decode the received signal
        auto [demodulated_data, success] = phy.demodulate(samples, modulated_length);

        // Check if decoding was successful
        if (success) {
            // Convert the decoded data back to a string
            std::string decoded_message(demodulated_data.begin(), demodulated_data.end());
            std::cout << "[main] Decoded message: " << decoded_message << std::endl;
        } else {
            std::cerr << "[main] Decoding failed. Transmission may be corrupted." << std::endl;
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
