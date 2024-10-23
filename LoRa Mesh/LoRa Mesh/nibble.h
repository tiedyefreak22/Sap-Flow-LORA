//
//  nibble.h
//  LoRa Mesh
//
//  Created by Kevin Hardin on 10/23/24.
//

#ifndef nibble_h
#define nibble_h

class nibble {
private:
    unsigned char value;

public:
    // Constructor to initialize the nibble
    nibble(unsigned int val = 0) {
        set(val);
    }

    // Set method to ensure value fits within 4 bits (0-15)
    void set(unsigned int val) {
        value = val & 0x0F;  // Mask the value to only keep the lower 4 bits
    }

    // Get the current nibble value
    unsigned char get() const {
        return value;
    }

    // Overload assignment operator for assigning unsigned int to nibble
    nibble& operator=(unsigned int val) {
        set(val);
        return *this;
    }

    // Overload the addition operator
    nibble operator+(const nibble& other) const {
        return nibble((this->value + other.value) & 0x0F);  // Keep only lower 4 bits
    }
};

#endif /* nibble_h */
