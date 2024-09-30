// main.cpp
#include <iostream>
#include "LoRaPHY.h"
#include "Aligntrack.h"

int main() {
    // Parameters for LoRa PHY layer
    double rf_freq = 915e6;    // Carrier frequency
    int sf = 7;                // Spreading factor
    double bw = 125e3;         // Bandwidth
    double fs = 1e6;           // Sampling rate

    init_rtl_sdr();

    uint8_t *buffer = receive_rtl_sdr(dev);
    rtlsdr_close(dev);

    // Initialize LoRaPHY object
    LoRaPHY phy(rf_freq, sf, bw, fs);

    // Sample data to transmit
    std::vector<int> data_to_send = {1, 2, 3, 4, 5};

    // Transmit data
    phy.transmit(data_to_send);

    // Receive and process data
    phy.receive();

    return 0;
}
