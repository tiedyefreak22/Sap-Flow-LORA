from numpy import floor, ceil, sqrt, prod, log10, log, sin, cos, pi, exp, tan, arcsin, arccos, arctan, zeros, ones, pi, sinc, isscalar, real, imag, arctan2, isinf, inf, mod, isnan, unique, linspace
import numpy as np
import math
from scipy.interpolate import interp1d, Akima1DInterpolator, PchipInterpolator
from scipy.constants import c

def isfloat(value):
    try:
        float(value)
        return True
    except (ValueError, TypeError):
        # Catches cases where conversion fails, e.g., "abc" or a non-string input
        return False

def deg2rad(angleInDegrees):
    return np.radians(angleInDegrees)

def rad2deg(angleInRadians):
    return np.degrees(angleInRadians)

def asin(array):
    return arcsin(array)

def acos(array):
    return arccos(array)

def atan(array):
    return arctan(array)

def asind(array):
    return rad2deg(asin(deg2rad(array)))
                   
def sind(array):
    return rad2deg(sin(deg2rad(array)))

def acosd(array):
    return rad2deg(acos(deg2rad(array)))
                   
def cosd(array):
    return rad2deg(cos(deg2rad(array)))

def atand(array):
    return rad2deg(atan(deg2rad(array)))
                   
def tand(array):
    return rad2deg(tan(deg2rad(array)))

def atan2(Y, X):
    return arctan2(Y, X)

def numel(array, axis = None):
    if axis == None:
        return np.size(array)
    else:
        return np.size(array, axis = axis - 1)

def size(array, axis = None):
    if axis == None:
        return np.shape(array)
    else:
        return np.size(array, axis = axis - 1)

def nan(sz1, sz2 = None, sz3 = None):
    if (sz2 is None) and (sz3 is None):
        return np.full(sz1, np.nan)
    elif (sz3 is None):
        return np.full((sz1, sz2), np.nan)
    else:
        return np.full((sz1, sz2, sz3), np.nan)

def angle(h):
    '''
    %ANGLE  Phase angle.
    %   ANGLE(H) returns the phase angles, in radians, of a matrix with
    %   complex elements.  
    %
    %   Class support for input X:
    %      float: double, single
    %
    %   See also ABS, UNWRAP.
    
    %   Copyright 1984-2010 The MathWorks, Inc. 
    '''
    p = atan2(imag(h), real(h))
    
    return p

def singcolvec(vector): # replaces "(:)" in Matlab
    if isscalar(vector):
        return [vector]
    else:
        return np.reshape(np.ravel(vector, order = "F"), (np.size(np.ravel(vector, order = "F")), 1))

def singrowvec(vector): # replaces "(:).'" in Matlab
    if isscalar(vector):
        return [vector]
    else:
        return np.ravel(vector, order = "F")

def single(value):
    return np.float32(value)

def eps(x = 1.0, datatype = "double", like = None):
    '''
       EPS('double') is the same as EPS, or EPS(1.0).
       EPS('single') is the same as EPS(single(1.0)), or single(2^-23).

    EPS  Spacing of floating point numbers.
       EPS(X) is the positive distance from ABS(X) to the next larger in
       magnitude floating point number of the same precision as X.
       X may be either double precision or single precision.
       For all X, EPS(X) is equal to EPS(ABS(X)).
    
       EPS, with no arguments, is the distance from 1.0 to the next larger double
       precision number, that is EPS with no arguments returns 2^(-52).
    
       EPS(CLASSNAME) returns the positive distance from 1.0 to the next 
       floating point number in the precision of CLASSNAME.
    
       EPS('double') is the same as EPS, or EPS(1.0).
       EPS('single') is the same as EPS(single(1.0)), or single(2^-23).
    
       EPS('like', Y) returns the positive distance from 1.0 to the next
       floating point number in the precision of the numeric variable Y, with
       the same data type, sparsity, and complexity (real or complex) as the
       numeric variable Y.
    
       Except for numbers whose absolute value is smaller than REALMIN,
       if 2^E <= ABS(X) < 2^(E+1), then
          EPS(X) returns 2^(E-23) if ISA(X,'single')
          EPS(X) returns 2^(E-52) if ISA(X,'double')
    
       For all X of class double such that ABS(X) <= REALMIN, EPS(X)
       returns 2^(-1074).   Similarly, for all X of class single such that
       ABS(X) <= REALMIN('single'), EPS(X) returns 2^(-149).
    
       Replace expressions of the form
          if Y < EPS * ABS(X)
       with
          if Y < EPS(X)
    '''
    if datatype.lower() == "double":
        return math.nextafter(x, inf)
    elif datatype.lower() == "single":
        return math.nextafter(single(x), inf)
    elif like is not None:
        match prototype.dtype:
            case np.float64:
                return math.nextafter(x, inf)
            case np.int64:
                return math.nextafter(x, inf)
            case np.float32:
                return math.nextafter(single(x), inf)
            case np.int32:
                return math.nextafter(single(x), inf)
            case _:
                raise ValueError("Value is not an single/double int nor float!")
    else:
        raise ValueError("Something is wrong with eps arguments")

def rem(arrayA, arrayB):
    return np.remainder(arrayA, arrayB)

def interp1(x, v, xq, method = "linear", extrap = False):
    match method.lower():
        case "linear":
            if extrap:
                F = interp1d(x, v, "linear", "extrapolate")
            else:
                F = interp1d(x, v, "linear")
        case "nearest":
            if extrap:
                F = interp1d(x, v, "nearest", "extrapolate")
            else:
                F = interp1d(x, v, "nearest")
        case "next":
            if extrap:
                F = interp1d(x, v, "next", "extrapolate")
            else:
                F = interp1d(x, v, "next")
        case "previous":
            if extrap:
                F = interp1d(x, v, "previous", "extrapolate")
            else:
                F = interp1d(x, v, "previous")
        case "pchip":
            F = PchipInterpolator(x, v, extrapolate = extrap)
        case "cubic":
            if extrap:
                F = interp1d(x, v, "cubic", "extrapolate")
            else:
                F = interp1d(x, v, "cubic")
        case "v5cubic":
            if extrap:
                F = interp1d(x, v, "cubic", "extrapolate")
            else:
                F = interp1d(x, v, "cubic")
        case "makima":
            F = Akima1DInterpolator(x, v, method = "makima", extrapolate = extrap)
        case "spline":
            if extrap:
                F = interp1d(x, v, "slinear", "extrapolate")
            else:
                F = interp1d(x, v, "slinear")
        case _:
            raise ValueError("Please enter one of the following methods: linear, nearest, next, previous, pchip, cubic, v5cubic, makima, or spline.")

    return F(xq)

def matlabmin(A, B = None, dim = None):
    if (B is None) and (dim is None):
        return np.min(A), np.argmin(A)
    elif (B == []) and isinstance(dim, str) and (dim.lower == "all"):
        return np.min(A), np.argmin(A)
    elif (B == []) and (dim is not None):
        match (dim, A.shape):
            case (0, _):
                raise ValueError("This function uses matlab's 1 indexing for input arguments, not python's 0 indexing (not applicable for output arguments)")
            case (dim, ()) if dim > 0:
                return A, 0
            case (1, (size,)):
                return A, ones(np.shape(A), dtype = "int")
            case (2, (size,)):
                return np.min(A, axis = dim - 2), np.argmin(A)
            case (3, (size,)):
                return A, ones(np.shape(A), dtype = "int")
            case (1, (rows, columns)):
                return np.min(A, axis = dim - 1), np.argmin(A, axis = dim - 1)
            case (2, (rows, columns)):
                return singcolvec(np.min(A, axis = dim - 1)), singcolvec(np.argmin(A, axis = dim - 1))
            case (3, (rows, columns)):
                return A, ones(np.shape(A), dtype = "int")
            case (1, (dim1, dim2, dim3)):
                return np.min(A, axis = dim - 1), np.argmin(A, axis = dim - 1)
            case (2, (dim1, dim2, dim3)):
                return singcolvec(np.min(A, axis = dim - 1)), singcolvec(np.argmin(A, axis = dim - 1))
            case (3, (dim1, dim2, dim3)):
                return np.minimum(A[:, :, 0], A[:, :, 1]), np.argmin(A, axis = dim - 1)
            case _:
                raise ValueError("This function cannot handle dimensions above 3")
    elif (B is not None) and (dim is None):
        return np.minimum(A, B)
    else:
        raise ValueError("Incompatible matlabmin arguments")

def matlabmax(A, B = None, dim = None):
    if (B is None) and (dim is None):
        return np.max(A), np.argmax(A)
    elif (B == []) and isinstance(dim, str) and (dim.lower == "all"):
        return np.max(A), np.argmax(A)
    elif (B == []) and (dim is not None):
        match (dim, A.shape):
            case (0, _):
                raise ValueError("This function uses matlab's 1 indexing for input arguments, not python's 0 indexing (not applicable for output arguments)")
            case (dim, ()) if dim > 0:
                return A, 0
            case (1, (size,)):
                return A, ones(np.shape(A), dtype = "int")
            case (2, (size,)):
                return np.max(A, axis = dim - 2), np.argmax(A)
            case (3, (size,)):
                return A, ones(np.shape(A), dtype = "int")
            case (1, (rows, columns)):
                return np.max(A, axis = dim - 1), np.argmax(A, axis = dim - 1)
            case (2, (rows, columns)):
                return singcolvec(np.max(A, axis = dim - 1)), singcolvec(np.argmax(A, axis = dim - 1))
            case (3, (rows, columns)):
                return A, ones(np.shape(A), dtype = "int")
            case (1, (dim1, dim2, dim3)):
                return np.max(A, axis = dim - 1), np.argmax(A, axis = dim - 1)
            case (2, (dim1, dim2, dim3)):
                return singcolvec(np.max(A, axis = dim - 1)), singcolvec(np.argmax(A, axis = dim - 1))
            case (3, (dim1, dim2, dim3)):
                return np.maximum(A[:, :, 0], A[:, :, 1]), np.argmax(A, axis = dim - 1)
            case _:
                raise ValueError("This function cannot handle dimensions above 3")
    elif (B is not None) and (dim is None):
        return np.maximum(A, B)
    else:
        raise ValueError("Incompatible matlabmax arguments")

def mag2db(y):
    '''
    MAG2DB  Magnitude to dB conversion.
    
       YDB = MAG2DB(Y) converts magnitude data Y into dB values.
       Negative values of Y are mapped to NaN.
    
       See also DB2MAG.
    
       Copyright 1986-2021 The MathWorks, Inc.
    '''
    y[y < 0] = np.nan
    ydb = 20 * log10(y)

    return ydb
    
def mag2dbScalar(y):
    # mag2db for scalar y.
    if y < 0:
       ydb = np.full_like(y, np.nan)
    else:
       ydb = 20 * log10(y)

    return ydb

def freq2wavelen(freq, c = c):
    '''
    freq2wavelen Convert frequency to wavelength
       WAVELEN = freq2wavelen(FREQ) returns the wavelength in m assuming the
       speed of light in a vacuum.
    
       FREQ is an M-length vector of positive frequencies in Hz. The output
       WAVELEN is an M-length vector of wavelengths in m. The default value of
       the speed of light in this case is obtained from
       physconst('LightSpeed').
    
       WAVELEN = freq2wavelen(FREQ,C) specifies the wave propagation speed C
       as a positive scalar in meters/second.
    
       [WAVELEN,C] = freq2wavelen(...) also returns the wave propagation speed
       in m/s, which was used for the calculation.
    
        Example:
          Calculate the wavelength corresponding to a frequency of 3 GHz. 
       lambda = freq2wavelen(3e9)
    
       See also phased, physconst, wavelen2freq.
    
       Copyright 2020 The MathWorks, Inc.
    '''
    # Calculate wavelength 
    wavelen = c / freq
    return wavelen, c

def wavelen2freq(freq, c = c):
    '''
    %wavelen2freq Convert wavelength to frequency
    %   FREQ = wavelen2freq(WAVELEN) returns the frequency in Hz assuming the
    %   speed of light in a vacuum.
    %
    %   WAVELEN is an M-length vector of positive wavelengths in m. The output
    %   FREQ is an M-length vector of frequencies in Hz. The default value of
    %   the speed of light in this case is obtained from
    %   physconst('LightSpeed').
    %
    %   FREQ = wavelen2freq(WAVELEN,C) specifies the wave propagation speed C
    %   as a positive scalar in meters/second.
    %
    %   [FREQ,C] = wavelen2freq(...) also returns the wave propagation speed in
    %   m/s, which was used for the calculation.
    %
    %   % Example:
    %   %   Calculate the frequency corresponding to a wavelength of 2 meters. 
    %   freq = wavelen2freq(2)
    %
    %   See also phased, physconst, freq2wavelen.
    
    %   Copyright 2020 The MathWorks, Inc.
    '''    
    # Calculate frequency
    freq = c / wavelen
    
    return freq, c

def sign(nums):
    if np.size(nums) == 1:
        if nums > 0:
            return 1
        elif nums < 0:
            return -1
        else:
            return 0
    else:
        output = []
        for num in nums:
            if num > 0:
                output.append(1)
            elif num < 0:
                output.append(-1)
            else:
                output.append(0)
        return output

def pow2db(y):
    '''
    Convert power/RCS to decibels
    (both work because they are relative to 1W or 1m)

    Args:
        y (np.array): N-element input vector array [W or m]

    Returns:
        np.array: Power measurement [dBW or dBsm]
    '''
    return 10 * log10(y)

def db2pow(y):
    '''
    Convert decibels to power/RCS
    (both work because they are relative to 1W or 1m)

    Args:
        y (np.array): N-element input vector array [dBW or dBsm]

    Returns:
        np.array: Power [W or m]
    '''
    return 10 ** (y / 10)

def num2str(array, precision = None, formatSpec = None):
    if isscalar(array):
        array = [array]
    if (precision is not None) and (formatSpec is not None):
        raise ValueError("Please only specify either precision or formatSpec")
    elif precision is not None:
        return [str(round(num, precision)) for num in array]
    elif formatSpec is not None:
        return [formatSpec % num for num in array]
    else:
        return [str(num) for num in array]

def strcmpi(array1, array2):
    shape1 = np.shape(array1)
    if isinstance(array1, list):
        lc_array1 = [i.lower() for i in np.ravel(array1)]
    elif isinstance(array1, str):
        lc_array1 = [array1.lower()]
    else:
        raise ValueError("Arguments must either be strings or arrays of strings")

    shape2 = np.shape(array2)
    raise ValueError("Input array shapes much currently be the same")
    
    if isinstance(array2, list):
        lc_array2 = [i.lower() for i in np.ravel(array2)]
    elif isinstance(array2, str):
        lc_array2 = [array2.lower()]
    else:
        raise ValueError("Arguments must either be strings or arrays of strings")

    output_array = []
    [output_array.append(i == j) for i, j in  zip(lc_array1, lc_array2)]
    return np.array(output_array).reshape((np.maximum(shape1, shape2)))