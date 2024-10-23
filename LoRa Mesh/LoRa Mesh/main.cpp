//
//  main.cpp
//  LoRa Mesh
//
//  Created by Kevin Hardin on 9/30/24.
//

#include <iostream>
#include "RTLSDR.h"
#include "AlignTrack.h"
#include "LoRa_Demodulation.h"
#include "LoRaPHY.h"

int load_data = 0;

int main(int argc, const char * argv[]) {
    try {
        uint8_t *rtlData;

        if (load_data) {
            // Load data from a .bin file
            const char* filename = "/Users/kevinhardin/Documents/GitHub/Becca-Sap-Flow-LORA/output.bin";
            size_t sampleCount = 0;

            Sample* samples = readBinaryFile(filename, &sampleCount);
            if (samples == NULL) {
                return 1; // Error reading file
            }

            // Print the first few samples
            for (size_t i = 0; i < sampleCount && i < 10; ++i) {
                printf("Sample %zu: I = %f, Q = %f\n", i, samples[i].I, samples[i].Q);
            }
        } else {
            // Capture data from RTL-SDR
            init_rtl_sdr();
            rtlData = receive_rtl_sdr(dev);
        }

        // Convert raw data to Signal format
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
