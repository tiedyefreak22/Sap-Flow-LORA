import galois
import re
from crcmod import mkCrcFun

def parse_polynomial_to_coeffs(poly_str):
    # Remove whitespace
    poly_str = poly_str.replace(" ", "")

    # Get non-numeric characters
    non_numeric_char = set(''.join([char for char in poly_str if char.isalpha()]))

    if len(non_numeric_char) > 1:
        raise ValueError("Polynomial must only contain one variable")

    non_numeric_char = non_numeric_char.pop()
    
    # Pattern to match: coefficient, x, power
    # Matches: -3x^2, 2x, 5, -x
    pattern = rf'([+-]?\d*)\*?{non_numeric_char}\^?(\d*)'
    
    # Find all terms
    terms = re.findall(pattern, poly_str)
    
    # Handle constant term (if any)
    const_pattern = rf'(?<![{non_numeric_char}\^])[+-]?\d+(?![{non_numeric_char}\^])'
    constant = re.findall(const_pattern, poly_str)
    
    # Determine the degree to know the list size
    max_power = 0
    parsed_terms = []
    
    for coeff, power in terms:
        if power == '':
            power = 1# if non_numeric_char in coeff or not coeff else 0
        else:
            power = int(power)
        
        if coeff == '' or coeff == '+':
            coeff = 1
        elif coeff == '-':
            coeff = -1
        else:
            coeff = int(coeff)

        parsed_terms.append((coeff, power))
        if power > max_power:
            max_power = power

    # Create coefficients list filled with zeros
    coeffs = [0] * (max_power + 1)
    for c, p in parsed_terms:
        coeffs[p] += c
    
    # Add constant if exists
    if constant:
        # Simple check, might need better logic for complex constant placement
        coeffs[0] = int(constant[-1])
        
    return coeffs

def coeffs_to_hex(coeffs):
    val = 0
    for i, coeff in enumerate(coeffs):
        val += coeff * (2**i)
    return hex(val)

def CRCGenerator(
    msg,
    Polynomial = None,
    InitialConditions = 0,
    DirectMethod = False,
    FinalXOR = 0,
    ChecksumsPerFrame = 1,
):
    if Polynomial is None:
        Polynomial = 'z^16 + z^12 + z^5 + 1'
    else:
        Polynomial = Polynomial

    coeffs = parse_polynomial_to_coeffs(Polynomial)
    hex_poly = coeffs_to_hex(coeffs)

    # return calculate_crc(msg, value)
    return mkCrcFun(hex_poly, initCrc = InitialConditions, xorOut=FinalXOR, rev=DirectMethod)

def bitget(A, bit):
    """
    Replicates the MATLAB bitget(A, bit) function in Python.
    bit_pos is 1-based, consistent with MATLAB.
    """
    # MATLAB's bit position 1 corresponds to Python's bit position 0.
    # The bitwise AND operation (A & (1 << (bit_pos - 1))) checks if the specific bit is set.
    # If the result is non-zero, the bit is 1; otherwise, it is 0.

    if bit > 0:
        return (A >> (bit - 1)) & 1
    else:
        raise ValueError("bitget: function uses matlab 1 indexing")

def calculate_crc(data: bytes, polynomial = 0x1021, init_value = 0xFFFF):
    """
    Calculates 16-bit CRC using a given polynomial.
    """
    crc = init_value
    for byte in data:
        # XOR byte into the MSB of the CRC
        crc ^= (byte << 8)
        
        # Process each bit
        for _ in range(8):
            if (crc & 0x8000):
                crc = (crc << 1) ^ polynomial
            else:
                crc = (crc << 1)
            # Maintain 16-bit
            crc &= 0xFFFF
            
    return hex(crc)

def gf(x):
    # Matlab's default gf() galois field is 2
    GF = galois.GF(2)
    return GF(x)