// LoRaPHY.cpp
#include "LoRaPHY.h"
#include <iostream>
#include <vector>
#include <fftw3.h>

//#ifdef __cplusplus
//extern "C" {
//#endif
//
//#include <liquid/liquid.h>  // liquid-dsp header
//
//#ifdef __cplusplus
//}
//#endif

LoRaPHY::LoRaPHY(double rf, int spreading_factor, double bandwidth, double sampling_rate)
    : rf_freq(rf), sf(spreading_factor), bw(bandwidth), fs(sampling_rate) {
    if (rtlsdr_open(&dev, 0) < 0) {
        std::cerr << "Failed to open RTL-SDR device" << std::endl;
        exit(1);
    }
}

LoRaPHY::~LoRaPHY() {
    rtlsdr_close(dev);
}

std::vector<double> LoRaPHY::modulate(const std::vector<int>& symbols) {
    // Modulation logic
    std::vector<double> signal;
    return signal;
}

std::vector<int> LoRaPHY::demodulate(const std::vector<double>& signal) {
    // Demodulation logic
    std::vector<int> symbols;
    return symbols;
}

std::vector<int> LoRaPHY::encode(const std::vector<int>& data) {
    // Encoding logic
    return data;
}

std::vector<int> LoRaPHY::decode(const std::vector<int>& symbols) {
    // Decoding logic
    return symbols;
}

void LoRaPHY::receive() {
    uint8_t buffer[16384];
    int n_read;
    rtlsdr_read_sync(dev, buffer, sizeof(buffer), &n_read);
    
    std::vector<double> iq_data(buffer, buffer + n_read);
    std::vector<int> symbols = demodulate(iq_data);
    std::vector<int> data = decode(symbols);

    for (int d : data) {
        std::cout << d << " ";
    }
    std::cout << std::endl;
}

void LoRaPHY::transmit(const std::vector<int>& data) {
    std::vector<int> encoded_data = encode(data);
    std::vector<double> signal = modulate(encoded_data);

    // Transmission logic
}
