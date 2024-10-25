//
//  RTLSDR.h
//  LoRa Mesh
//
//  Created by Kevin Hardin on 9/30/24.
//

#ifndef RTLSDR_h
#define RTLSDR_h

#include <rtl-sdr.h>
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>
#include <iostream>
#include <cmath>
#include <fstream>

// Define a structure to hold I/Q samples
typedef struct {
    float I; // In-phase component
    float Q; // Quadrature component
} Sample;

Sample    add    (Sample l, Sample r)    {return ({Sample z = {l.I + r.I, l.Q + r.Q}; z;});}

#define SAMPLES 2048  // Number of samples to capture in each chunk
#define THRESHOLD 1000 // Magnitude threshold for peak detection
double NOISE_THRESHOLD = 45.0; // Threshold for detecting meaningful signals (adjust according to your environment)

#define DEFAULT_SAMPLE_RATE 1000000  // 1.0 MSPS
#define DEFAULT_FREQUENCY  915000000 // 915 MHz (LoRA)
//#define DEFAULT_SAMPLE_RATE 2048000  // 2.048 MSPS
//#define DEFAULT_FREQUENCY  1090000000 // 1.09 GHz (ADS-B)
//#define DEFAULT_SAMPLE_RATE 2048000  // 2.048 MSPS
//#define DEFAULT_FREQUENCY  94100000 // 94.1 MHz (Local Radio)

rtlsdr_dev_t *dev;

// Function to save the RTL-SDR data to a file
void save_data_to_file(const char* filename, double* freq, double* magnitude, int N) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    for (int i = 0; i < N; ++i) {
        outfile << freq[i] << " " << magnitude[i] << std::endl;
    }

    outfile.close();
    std::cout << "Data saved successfully to " << filename << std::endl;
}

// 1a. Initialize RTL-SDR
rtlsdr_dev_t *init_rtl_sdr() {
    int device_index = 0; // Index for the RTL-SDR device (use 0 if there's only one)
    // Open the RTL-SDR device
    if (rtlsdr_open(&dev, device_index) != 0) {
        fprintf(stderr, "Error opening RTL-SDR device\n");
        return NULL;
    }

    // Set center frequency
    if (rtlsdr_set_center_freq(dev, DEFAULT_FREQUENCY) != 0) {
        fprintf(stderr, "Error setting frequency\n");
        rtlsdr_close(dev);
        return NULL;
    }
    printf("Tuned to %u Hz.\n", DEFAULT_FREQUENCY);

    // Set sample rate
    if (rtlsdr_set_sample_rate(dev, DEFAULT_SAMPLE_RATE) != 0) {
        fprintf(stderr, "Error setting sample rate\n");
        rtlsdr_close(dev);
        return NULL;
    }
    printf("Sample rate set to %u Hz.\n", DEFAULT_SAMPLE_RATE);

    // Set the tuner gain mode to automatic
    if (rtlsdr_set_tuner_gain_mode(dev, 0) != 0) {
        fprintf(stderr, "Error setting gain mode\n");
        rtlsdr_close(dev);
        return NULL;
    }
    printf("Tuner gain set to auto.\n");
    return dev;
}

// 1b. RTL-SDR Signal Capture
void callback(unsigned char *buf, uint32_t len, void *ctx) {
    // This function is called each time a block of data is available
    double power_sum = 0;
    for (uint32_t i = 0; i < len; i++) {
        power_sum += (buf[i] - 127) * (buf[i] - 127);  // Normalize and square the sample values
    }

    double rms_power = sqrt(power_sum / len);  // RMS power

    if (rms_power > NOISE_THRESHOLD) {
        printf("Received %u bytes of data with power: %f\n", len, rms_power);

        // Append data to file
        // save_data_to_file("/Users/kevinhardin/Documents/rtl_sdr_output.dat");
    }
}

std::vector<Sample> receive_rtl_sdr(rtlsdr_dev_t *dev) {
    // Reset buffer before starting
    if (rtlsdr_reset_buffer(dev) != 0) {
        fprintf(stderr, "Error resetting buffer\n");
        rtlsdr_close(dev);
        //return NULL;
    }

    // Read samples (blocking mode, synchronous)
    std::vector<Sample>* buffer = (std::vector<Sample> *)malloc(SAMPLES * 2 * sizeof(uint8_t));
    if (buffer == NULL) {
        fprintf(stderr, "Memory allocation error\n");
        rtlsdr_close(dev);
        //return NULL;
    }

    int r = rtlsdr_read_async(dev, *callback, NULL, 0, 0);
    if (r < 0) {
        fprintf(stderr, "Error starting asynchronous read\n");
        rtlsdr_close(dev);
        //return NULL;
    }
    return *buffer;
}

// Function to read a binary file containing I/Q samples
std::vector<Sample> readBinaryFile(const char* filename, size_t* sampleCount) {
    // Open the binary file
    FILE* file = fopen(filename, "rb");
    if (!file) {
        perror("Could not open file");
        //return NULL;
    }

    // Move the file pointer to the end of the file
    fseek(file, 0, SEEK_END);
    long fileSize = ftell(file);
    rewind(file); // Move the file pointer back to the beginning

    // Calculate the number of samples (each sample has 2 floats)
    *sampleCount = fileSize / (sizeof(float) * 2);
    
    // Allocate memory for samples
    std::vector<Sample>* samples = (std::vector<Sample>*)malloc(*sampleCount * sizeof(Sample));
    if (!samples) {
        perror("Memory allocation failed");
        fclose(file);
        //return NULL;
    }

    // Read the samples from the file
    size_t readCount = fread(samples, sizeof(Sample), *sampleCount, file);
    if (readCount != *sampleCount) {
        perror("Error reading samples from file");
        free(samples);
        fclose(file);
        //return NULL;
    }

    fclose(file); // Close the file
    return *samples; // Return the pointer to the samples
}

#endif /* RTLSDR_h */
