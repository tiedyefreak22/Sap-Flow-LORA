/* 
 * AlignTrack Lora Mesh
 * Author: Kevin Hardin
 * Date: 9/18/24
 */

#include <rtl-sdr.h>
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>
#include <iostream>
#include <cmath>
#include <fstream>

#define SAMPLES 2048  // Number of samples to capture in each chunk
//#define DEFAULT_SAMPLE_RATE 1000000  // 1.0 MSPS
//#define DEFAULT_FREQUENCY  868000000 // 868 MHz (LoRA)
#define DEFAULT_SAMPLE_RATE 2048000  // 2.048 MSPS
#define DEFAULT_FREQUENCY  1090000000 // 1.09 GHz (ADS-B)
//#define DEFAULT_SAMPLE_RATE 2048000  // 2.048 MSPS
//#define DEFAULT_FREQUENCY  94100000 // 94.1 MHz (Local Radio)

rtlsdr_dev_t *dev = NULL;

// Function to save the FFT result to a file with error handling
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

// 3. Signal Separation
// Pseudo-code
// void separate_collisions(fftw_complex *frequency_data, int N) {
//     // Step 3a: Identify different signal components based on frequency peaks
//     for (int i = 0; i < N; i++) {
//         if (is_peak(frequency_data[i])) {
//             // Step 3b: Separate each signal component based on time-frequency shifts
//             // You can use cross-correlation to determine timing offsets.
//             double time_offset = calculate_time_offset(frequency_data[i]);
//             double frequency_shift = calculate_frequency_shift(frequency_data[i]);

//             // Step 3c: Align the signals by compensating for time and frequency shifts
//             align_signals(frequency_data[i], time_offset, frequency_shift);
//         }
//     }
// }

// 4. LoRa Demodulation
// void lora_demodulate(fftw_complex *aligned_signal) {
    // Step 4a: Correlate the aligned signal with known LoRa chirps
    // Step 4b: Decode the symbols based on the frequency shifts
    // Step 4c: Apply FEC (Forward Error Correction) to correct any errors
// }

// 5. Error Correction & Decoding
// AlignTrack also uses error correction mechanisms like FEC (Forward Error Correction)
// to ensure that any errors caused by collisions are corrected. LoRa includes built-in
// FEC using Hamming codes or Reed-Solomon codes, so you’ll implement or use a library
// for FEC.

// 1b. RTL-SDR Signal Capture
void callback(unsigned char *buf, uint32_t len, void *ctx) {
    // This function is called each time a block of data is available
    double power_sum = 0;
    for (uint32_t i = 0; i < len; i++) {
        power_sum += (buf[i] - 127) * (buf[i] - 127);  // Normalize and square the sample values
    }

    double rms_power = sqrt(power_sum / len);  // RMS power

    // Threshold for detecting meaningful signals (adjust according to your environment)
    double noise_threshold = 3.0;

    if (rms_power > noise_threshold) {
        printf("Received %u bytes of data with power: %f\n", len, rms_power);

        // Process and report the data
        
        // Perform time-frequency analysis
        fftw_complex *frequency_data = analyze_time_frequency(buf);


        // Calculate magnitude and frequency values
        double magnitude[SAMPLES];
        double freq[SAMPLES];

        for (int i = 0; i < SAMPLES / 2; ++i) {
            double real = frequency_data[i][0];
            double imag = frequency_data[i][1];
            freq[i] = i * DEFAULT_SAMPLE_RATE / SAMPLES;  // Frequency axis
            magnitude[i] = sqrt(real*real + imag*imag);  // Magnitude
            printf("Real: %f, Imag: %f, Mag: %f\n", real, imag, magnitude[i]);
        }

        // Save data to file
        save_data_to_file("/Users/kevinhardin/Documents/fft_output.dat", freq, magnitude, SAMPLES / 2);

        // Separate the colliding signals
        // separate_collisions(frequency_data, SAMPLES);

        // Demodulate the separated signals
        // lora_demodulate(frequency_data);
        
        // Error Correction & Decoding
    }
}

uint8_t *receive_rtl_sdr(rtlsdr_dev_t *dev) {
    // Reset buffer before starting
    if (rtlsdr_reset_buffer(dev) != 0) {
        fprintf(stderr, "Error resetting buffer\n");
        rtlsdr_close(dev);
        return NULL;
    }

    // Read samples (blocking mode, synchronous)
    uint8_t *buffer = (uint8_t *)malloc(SAMPLES * 2 * sizeof(uint8_t));
    if (buffer == NULL) {
        fprintf(stderr, "Memory allocation error\n");
        rtlsdr_close(dev);
        return NULL;
    }

    int r = rtlsdr_read_async(dev, *callback, NULL, 0, 0);
    if (r < 0) {
        fprintf(stderr, "Error starting asynchronous read\n");
        rtlsdr_close(dev);
        return NULL;
    }
    return buffer;
}

// setup() runs once, when the device is first turned on
void setup() {
    init_rtl_sdr();
}

// loop() runs over and over again, as quickly as it can execute.
void loop() {
    uint8_t *buffer = receive_rtl_sdr(dev);
    rtlsdr_close(dev);
}

int main() {
    setup();
    loop();
}
