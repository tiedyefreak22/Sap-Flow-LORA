//
//  AlignTrack.h
//  LoRa Mesh
//
//  Created by Kevin Hardin on 9/30/24.
//

#ifndef AlignTrack_h
#define AlignTrack_h

#include <complex>
#include <vector>
#include <iostream>
#include <fftw3.h>

//#include <rtl-sdr.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <fftw3.h>
//#include <math.h>
//#include <iostream>
//#include <cmath>
//#include <fstream>

// 2. Preprocessing: Time-Frequency Analysis
fftw_complex *analyze_time_frequency(uint8_t *buffer) {
    int N = SAMPLES;
    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_plan p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Prepare the input data (convert raw I/Q samples into complex numbers)
    // The buffer contains len bytes of raw I/Q samples
    // Each sample is 2 bytes: I (in-phase) and Q (quadrature)
    for (int i = 0; i < N; i++) {
        in[i][0] = buffer[2*i];   // I (real)
        in[i][1] = buffer[2*i+1]; // Q (imaginary)
    }
    fftw_execute(p);  // Execute the FFT to transform to frequency domain
    fftw_destroy_plan(p);
    return out;
}

// 3. Signal separation and inverse FFT back to time domain, then demodulate
fftw_complex *separate_collisions(fftw_complex *frequency_data) {
    // Step 3a: Identify different signal components based on frequency peak magnitudes
    for (uint32_t peak_index = 0; peak_index < SAMPLES; peak_index++) {
        double magnitude = sqrt(frequency_data[peak_index][0] * frequency_data[peak_index][0] + frequency_data[peak_index][1] * frequency_data[peak_index][1]);
        if (magnitude > THRESHOLD) {
            // printf("Peak detected at index %d with magnitude %f\n", peak_index, magnitude);
            
            // Step 3b: Separate each signal component based on time-frequency shifts
            // You can use cross-correlation to determine timing offsets.
            int i;
//            double frequency_shift = peak_index * (DEFAULT_SAMPLE_RATE / SAMPLES);
//            double time_offset = peak_index * (1.0 / SAMPLES); // Calculate frequency shift
            double frequency_shift = 0.0;
            double time_offset = 0.0;

            // printf("Separating signal at frequency shift: %f\n", frequency_shift);
            
            // Step 3c: Align the signals by compensating for time and frequency shifts
            // Dechirping by applying conjugate of ideal chirp
            for (i = 0; i < SAMPLES; i++) {
                frequency_data[i][0] *= cos(-2 * M_PI * frequency_shift * i);  // Real part
                frequency_data[i][1] *= sin(-2 * M_PI * frequency_shift * i);  // Imaginary part
            }
        }
    }
    // Apply inverse FFT after dechirping to obtain symbol values
    fftw_complex *demodulated_output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * SAMPLES);
    fftw_plan p = fftw_plan_dft_1d(SAMPLES, frequency_data, demodulated_output, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    // Process `out` which now contains demodulated symbols
    printf("Signal separated and demodulated successfully\n");

    fftw_destroy_plan(p);
    return demodulated_output;
}

// Function to adjust the phase of a signal
void adjustPhase(std::vector<std::complex<float>>& signal, float phaseOffset) {
    for (size_t i = 0; i < signal.size(); ++i) {
        float phase = std::arg(signal[i]);
        float magnitude = std::abs(signal[i]);
        signal[i] = std::polar(magnitude, phase + phaseOffset);
    }
}

// Function to calculate the phase offset from the cross-correlation peak
float calculatePhaseOffset(const std::vector<std::complex<float>>& crossCorr) {
    size_t peakIndex = 0;
    float maxMagnitude = 0.0;

    for (size_t i = 0; i < crossCorr.size(); ++i) {
        float magnitude = std::abs(crossCorr[i]);
        if (magnitude > maxMagnitude) {
            maxMagnitude = magnitude;
            peakIndex = i;
        }
    }

    return std::arg(crossCorr[peakIndex]);
}

// Function to compute the energy of a signal
float computeEnergy(const std::vector<std::complex<float>>& signal) {
    float energy = 0.0f;
    for (const auto& sample : signal) {
        energy += std::norm(sample);
    }
    return energy;
}

// Function to maximize energy between up to three signals
//std::vector<std::complex<float>> maximizeEnergy(const std::vector<std::vector<std::complex<float>>>& signals) {
//    std::vector<std::complex<float>> bestSignal;
//    float maxEnergy = 0.0;
//
//    for (size_t i = 0; i < signals.size(); ++i) {
//        auto currentSignal = signals[i];
//
//        // Adjust phase by comparing with other signals
//        for (size_t j = 0; j < signals.size(); ++j) {
//            if (i != j) {
//                std::vector<std::complex<float>> crossCorr = crossCorrelate(signals[i], signals[j]);
//                float phaseOffset = calculatePhaseOffset(crossCorr);
//                adjustPhase(currentSignal, phaseOffset);
//            }
//        }
//
//        // Compute energy and keep track of the best signal
//        float energy = computeEnergy(currentSignal);
//        if (energy > maxEnergy) {
//            maxEnergy = energy;
//            bestSignal = currentSignal;
//        }
//    }
//
//    return bestSignal;
//}

// AlignTrack deconfliction process for up to three overlapping signals
//std::vector<std::vector<std::complex<float>>> alignTrack(const std::vector<std::vector<std::complex<float>>>& overlappingSignals) {
//    std::vector<std::vector<std::complex<float>>> deconflictedSignals;
//
//    for (size_t i = 0; i < overlappingSignals.size(); ++i) {
//        auto signal = overlappingSignals[i];
//
//        // Time synchronization using cross-correlation
//        for (size_t j = 0; j < i; ++j) {
//            std::vector<std::complex<float>> crossCorr = crossCorrelate(overlappingSignals[i], overlappingSignals[j]);
//            float phaseOffset = calculatePhaseOffset(crossCorr);
//            adjustPhase(signal, phaseOffset);
//        }
//
//        // Maximize energy using iterative adjustments for all signals
//        signal = maximizeEnergy(overlappingSignals);
//
//        // Store the deconflicted signal
//        deconflictedSignals.push_back(signal);
//    }
//
//    return deconflictedSignals;
//}

#endif /* AlignTrack_h */
