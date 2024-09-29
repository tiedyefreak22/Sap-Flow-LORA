// LoRaPHY.h
#ifndef LORAPHY_H
#define LORAPHY_H

#include <vector>
#include <complex>
#include <rtl-sdr.h>

class LoRaPHY {
public:
    double rf_freq;
    int sf;
    double bw;
    double fs;
    rtlsdr_dev_t* dev;

    LoRaPHY(double rf, int spreading_factor, double bandwidth, double sampling_rate);
    ~LoRaPHY();

    std::vector<double> modulate(const std::vector<int>& symbols);
    std::vector<int> demodulate(const std::vector<double>& signal);
    std::vector<int> encode(const std::vector<int>& data);
    std::vector<int> decode(const std::vector<int>& symbols);
    void receive();
    void transmit(const std::vector<int>& data);
};

#endif // LORAPHY_H
