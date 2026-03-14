from numpy import floor, ceil, sqrt, prod, log10, log, sin, cos, pi, exp, tan, arcsin, arccos, arctan, zeros, ones, pi, sinc, isscalar, real, imag, arctan2, isinf, inf, mod, isnan, linspace, isfinite, cumsum, fliplr, copysign, nan, logspace, hypot
from numpy.random import rand
from scipy.special import iv
import scipy.special as sps
from scipy.integrate import quad
import scipy.integrate as spi
import numpy as np
import math
from scipy.interpolate import interp1d, Akima1DInterpolator, PchipInterpolator
from scipy.constants import c
import sys
import matplotlib.pyplot as plt
import matplotlib
from numbers import Number

# Get the smallest positive normalized float
realmin = sys.float_info.min

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
    output = asin(deg2rad(array))
    return np.where(abs(output - np.round(output)) < sqrt(eps(1)), np.round(output).astype(np.int_), output)
    
def sind(array):
    output = sin(deg2rad(array))
    return np.where(abs(output - np.round(output)) < sqrt(eps(1)), np.round(output).astype(np.int_), output)
    
def acosd(array):
    output = acos(deg2rad(array))
    return np.where(abs(output - np.round(output)) < sqrt(eps(1)), np.round(output).astype(np.int_), output)
    
def cosd(array):
    output = cos(deg2rad(array))
    return np.where(abs(output - np.round(output)) < sqrt(eps(1)), np.round(output).astype(np.int_), output)

def atand(array):
    output = atan(deg2rad(array))
    return np.where(abs(output - np.round(output)) < sqrt(eps(1)), np.round(output).astype(np.int_), output)
    
def tand(array):
    output = tan(deg2rad(array))
    return np.where(abs(output - np.round(output)) < sqrt(eps(1)), np.round(output).astype(np.int_), output)
    
def atan2(Y, X):
    return arctan2(Y, X)

# def numel(array, axis = None):
#     if axis == None:
#         return np.size(array)
#     else:
#         return np.size(array, axis = axis - 1)

def numel(array):
    return np.size(np.ravel(array))

def size(array, axis = None):
    ndim = np.ndim(array)
    if numel(axis) == 1:
        match (ndim, axis):
            case (_, 0):
                raise ValueError("This function uses matlab's 1 indexing for input arguments, not python's 0 indexing (not applicable for output arguments)")
            case (0, None):
                return tuple(np.ravel(ones((1, 2), dtype = np.int_)))
            case (0, _):
                return unlevel(ones(1, dtype = np.int_))
            case (1, None):
                return tuple((1, np.size(array)))
            case (1, 1):
                return unlevel(ones(1, dtype = np.int_))
            case (1, 2):
                return np.size(array)
            case (1, _):
                return unlevel(ones(1, dtype = np.int_))
            case (2, None):
                if (np.ndim(array) == 1) or ((np.ndim(array) == 2) and (np.shape(array)[1] == 1)): # single column
                    return tuple((np.size(array), 1))
                else:
                    return np.shape(array)
            case (2, 1):
                if (np.ndim(array) == 1) or ((np.ndim(array) == 2) and (np.shape(array)[1] == 1)): # single column
                    return np.size(array)
                else:
                    return np.size(array, axis = axis - 1)
            case (2, 2):
                if (np.ndim(array) == 1) or ((np.ndim(array) == 2) and (np.shape(array)[1] == 1)): # single column
                    return unlevel(ones(1, dtype = np.int_))
                else:
                    return np.size(array, axis = axis - 1)
            case (2, _):
                return unlevel(ones(1, dtype = np.int_))
            case (3, None):
                return np.shape(array)
            case (3, (1 | 2 | 3)):
                return np.size(array, axis = axis - 1)
            case (3, _):
                return unlevel(ones(1, dtype = np.int_))
            case (_, _):
                raise ValueError("This function can only go up to 3 dimensions!")

def Nan(sz1, sz2 = None, sz3 = None):
    if (sz2 is None) and (sz3 is None):
        return np.full(sz1, np.nan)
    elif (sz3 is None):
        return np.full((sz1, sz2), np.nan)
    else:
        return np.full((sz1, sz2, sz3), np.nan)

def Inf(sz1, sz2 = None, sz3 = None):
    if (sz2 is None) and (sz3 is None):
        return np.full(sz1, inf)
    elif (sz3 is None):
        return np.full((sz1, sz2), inf)
    else:
        return np.full((sz1, sz2, sz3), inf)

def unlevel(obj):
    while isinstance(obj, (list, np.ndarray)) and len(obj) == 1:
        obj = obj[0]
    return obj

def numpyify(x):
    if isinstance(x, str):
        return x
    elif isinstance(x, list):
        return np.asarray(unlevel(x))
    elif isinstance(x, np.ndarray):
        return np.asarray(unlevel(x.tolist()))
    elif np.ndim(x) == 0:
        return x
    else:
        print("numpyify error")

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
    if numel(x) == 1:
        if datatype.lower() == "double":
            return math.nextafter(x, inf) - x
        elif datatype.lower() == "single":
            return math.nextafter(single(x), inf) - x
        elif like is not None:
            match prototype.dtype:
                case np.float64:
                    return math.nextafter(x, inf) - x
                case np.int64:
                    return math.nextafter(x, inf) - x
                case np.float32:
                    return math.nextafter(single(x), inf) - x
                case np.int32:
                    return math.nextafter(single(x), inf) - x
                case _:
                    raise ValueError("Value is not an single/double int nor float!")
        else:
            raise ValueError("Something is wrong with eps arguments")
    else:
        output = []
        for value in x:
            if datatype.lower() == "double":
                output.append(math.nextafter(value, inf) - value)
            elif datatype.lower() == "single":
                output.append(math.nextafter(single(value), inf) - value)
            elif like is not None:
                match prototype.dtype:
                    case np.float64:
                        output.append(math.nextafter(value, inf) - value)
                    case np.int64:
                        output.append(math.nextafter(value, inf) - value)
                    case np.float32:
                        output.append(math.nextafter(single(value), inf) - value)
                    case np.int32:
                        output.append(math.nextafter(single(value), inf) - value)
                    case _:
                        raise ValueError("Value is not an single/double int nor float!")
            else:
                raise ValueError("Something is wrong with eps arguments")
        return np.array(output)

def rem(arrayA, arrayB):
    return np.remainder(arrayA, arrayB)

def interp1(x, v_array, xq, method = "linear", extrap = None):
    '''
    Sample points, specified as a row or column vector of real numbers. The values in x must be distinct. The length of x must conform to one of the following requirements:
        If v is a vector, then length(x) must equal length(v).
        If v is an array, then length(x) must equal size(1, v).

    Sample values, specified as a vector, matrix, or array of real or complex numbers. If v is a matrix or an array, then each row contains a separate set of 1-D values.
    '''
    def interp1_helper(x, v, xq, method = "linear", extrap = None):
        # sort both x and v based on x
        paired_sorted = sorted(zip(x, v))
        x, v = zip(*paired_sorted)
        x = list(x)
        v = list(v)
        
        # check if x is strictly increasing; if it isn't, throw an error
        if not all(x[i] < x[i + 1] for i in range(len(x) - 1)):
            raise ValueError("interp1: x must be distinct values")
            
        match method.lower():
            case "linear":
                if extrap: # extrapolate if extrap is True
                    F = interp1d(x, v, "linear", fill_value = "extrapolate", bounds_error = False)
                else: # if extrapolate is false or None, don't extrapolate
                    F = interp1d(x, v, "linear", bounds_error = False)
            case "nearest":
                if extrap:
                    F = interp1d(x, v, "nearest", fill_value = "extrapolate", bounds_error = False)
                else:
                    F = interp1d(x, v, "nearest", bounds_error = False)
            case "next":
                if extrap:
                    F = interp1d(x, v, "next", fill_value = "extrapolate", bounds_error = False)
                else:
                    F = interp1d(x, v, "next", bounds_error = False)
            case "previous":
                if extrap:
                    F = interp1d(x, v, "previous", fill_value = "extrapolate", bounds_error = False)
                else:
                    F = interp1d(x, v, "previous", bounds_error = False)
            case "pchip": # extrap by default; returns nans by default if not extrapolating
                if extrap is not None: # if extrapolate is specified, do that
                    F = PchipInterpolator(x, v, extrapolate = extrap)
                else: # if extrapolate is not specified, do it by default
                    F = PchipInterpolator(x, v, extrapolate = True)
            case "cubic":
                if extrap:
                    F = interp1d(x, v, "cubic", fill_value = "extrapolate", bounds_error = False)
                else:
                    F = interp1d(x, v, "cubic", bounds_error = False)
            case "v5cubic":
                if extrap:
                    F = interp1d(x, v, "cubic", fill_value = "extrapolate", bounds_error = False)
                else:
                    F = interp1d(x, v, "cubic", bounds_error = False)
            case "makima": # extrap by default; returns nans by default if not extrapolating
                if extrap is not None:
                    F = Akima1DInterpolator(x, v, method = "makima", extrapolate = extrap)
                else:
                    F = Akima1DInterpolator(x, v, method = "makima", extrapolate = True)
            case "spline": # extrap by default
                if extrap is not None:
                    if extrap:
                        F = interp1d(x, v, "slinear", fill_value = "extrapolate", bounds_error = False)
                    else:
                        F = interp1d(x, v, "slinear", bounds_error = False)
                else:
                    F = interp1d(x, v, "slinear", fill_value = "extrapolate", bounds_error = False)
            case _:
                raise ValueError("Please enter one of the following methods: linear, nearest, next, previous, pchip, cubic, v5cubic, makima, or spline.")
    
        return F(xq)

    # interp1_helper is the guts of the interp function; the actual interp1 function just handles looping through v if required
    output_array = []
    if np.ndim(v_array) > 1:
        for i in range(np.shape(v_array)[0]):
            output_array.append(interp1_helper(singrowvec(x), v_array[i, :], xq, method = method, extrap = extrap))
        return output_array
    else:
        return interp1_helper(singrowvec(x), v_array, xq, method = method, extrap = extrap)

def matlabmin(A, B = None, dim = None):
    if (B is None) and (dim is None):
        return np.nanmin(A), np.nanargmin(A)
    elif (B == []) and isinstance(dim, str) and (dim.lower() == "all"):
        return np.nanmin(A), np.nanargmin(A)
    elif (B == []) and (dim is not None):
        match (dim, A.shape):
            case (0, _):
                raise ValueError("This function uses matlab's 1 indexing for input arguments, not python's 0 indexing (not applicable for output arguments)")
            case (dim, ()) if dim > 0:
                return A, 0
            case (1, (size,)):
                return A, ones(np.shape(A), dtype = "int")
            case (2, (size,)):
                return np.nanmin(A, axis = dim - 2), np.nanargmin(A)
            case (3, (size,)):
                return A, ones(np.shape(A), dtype = "int")
            case (1, (rows, columns)):
                return np.nanmin(A, axis = dim - 1), np.nanargmin(A, axis = dim - 1)
            case (2, (rows, columns)):
                return singcolvec(np.nanmin(A, axis = dim - 1)), singcolvec(np.nanargmin(A, axis = dim - 1))
            case (3, (rows, columns)):
                return A, ones(np.shape(A), dtype = "int")
            case (1, (dim1, dim2, dim3)):
                return np.nanmin(A, axis = dim - 1), np.nanargmin(A, axis = dim - 1)
            case (2, (dim1, dim2, dim3)):
                return singcolvec(np.nanmin(A, axis = dim - 1)), singcolvec(np.nanargmin(A, axis = dim - 1))
            case (3, (dim1, dim2, dim3)):
                return np.fmin(A[:, :, 0], A[:, :, 1]), np.nanargmin(A, axis = dim - 1)
                # return np.minimum(A[:, :, 0], A[:, :, 1]), np.nanargmin(A, axis = dim - 1)
            case _:
                raise ValueError("This function cannot handle dimensions above 3")
    elif (B is not None) and (dim is None):
        return np.fmin(A, B)
        # return np.minimum(A, B)
    else:
        raise ValueError("Incompatible matlabmin arguments")

def matlabmax(A, B = None, dim = None):
    if (B is None) and (dim is None):
        return np.nanmax(A), np.nanargmax(A)
    elif (B == []) and isinstance(dim, str) and (dim.lower() == "all"):
        return np.nanmax(A), np.nanargmax(A)
    elif (B == []) and (dim is not None):
        match (dim, A.shape):
            case (0, _):
                raise ValueError("This function uses matlab's 1 indexing for input arguments, not python's 0 indexing (not applicable for output arguments)")
            case (dim, ()) if dim > 0:
                return A, 0
            case (1, (size,)):
                return A, ones(np.shape(A), dtype = "int")
            case (2, (size,)):
                return np.nanmax(A, axis = dim - 2), np.nanargmax(A)
            case (3, (size,)):
                return A, ones(np.shape(A), dtype = "int")
            case (1, (rows, columns)):
                return np.nanmax(A, axis = dim - 1), np.nanargmax(A, axis = dim - 1)
            case (2, (rows, columns)):
                return singcolvec(np.nanmax(A, axis = dim - 1)), singcolvec(np.nanargmax(A, axis = dim - 1))
            case (3, (rows, columns)):
                return A, ones(np.shape(A), dtype = "int")
            case (1, (dim1, dim2, dim3)):
                return np.nanmax(A, axis = dim - 1), np.nanargmax(A, axis = dim - 1)
            case (2, (dim1, dim2, dim3)):
                return singcolvec(np.nanmax(A, axis = dim - 1)), singcolvec(np.nanargmax(A, axis = dim - 1))
            case (3, (dim1, dim2, dim3)):
                return np.fmax(A[:, :, 0], A[:, :, 1]), np.nanargmax(A, axis = dim - 1)
                # return np.maximum(A[:, :, 0], A[:, :, 1]), np.nanargmax(A, axis = dim - 1)
            case _:
                raise ValueError("This function cannot handle dimensions above 3")
    elif (B is not None) and (dim is None):
        return np.fmax(A, B)
        # return np.maximum(A, B)
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

def db2mag(ydb):
    '''
    %DB2MAG  dB to magnitude conversion.
    %
    %   Y = DB2MAG(YDB) computes Y such that YDB = 20*log10(Y).
    
    %   Copyright 1986-2007 The MathWorks, Inc.
    %   $Revision $
    '''
    y = 10 ** (ydb / 20)
    return y

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
    if isinstance(array1, (list, np.ndarray)):
        lc_array1 = [i.lower() for i in np.ravel(array1)]
    elif isinstance(array1, str):
        lc_array1 = [array1.lower()]
    else:
        raise ValueError(f"Arguments must either be strings or arrays of strings; not {type(array1)}.")

    shape2 = np.shape(array2)
    if isinstance(array2, (list, np.ndarray)):
        lc_array2 = [i.lower() for i in np.ravel(array2)]
    elif isinstance(array2, str):
        lc_array2 = [array2.lower()]
    else:
        raise ValueError(f"Arguments must either be strings or arrays of strings; not {type(array2)}.")

    if shape1 != shape2:
        if numel(array1) == 1:
            lc_array1 = np.ravel(np.full(shape2, unlevel(array1), dtype = object))
            shape1 = shape2
        elif numel(array2) == 1:
            lc_array2 = np.ravel(np.full(shape1, unlevel(array2), dtype = object))
            shape2 = shape1
        else:
            raise ValueError("Input array shapes must currently be the same or one must only be length 1")

    output_array = []
    [output_array.append(i == j) for i, j in  zip(lc_array1, lc_array2)]
    
    return np.array(output_array).reshape((np.maximum(shape1, shape2)))

def strcmp(array1, array2):
    shape1 = np.shape(array1)
    if isinstance(array1, (list, np.ndarray)):
        lc_array1 = [i for i in np.ravel(array1)]
    elif isinstance(array1, str):
        lc_array1 = [array1]
    else:
        raise ValueError(f"Arguments must either be strings or arrays of strings; not {type(array1)}.")

    shape2 = np.shape(array2)
    if isinstance(array2, (list, np.ndarray)):
        lc_array2 = [i for i in np.ravel(array2)]
    elif isinstance(array2, str):
        lc_array2 = [array2]
    else:
        raise ValueError(f"Arguments must either be strings or arrays of strings; not {type(array2)}.")

    if shape1 != shape2:
        if numel(array1) == 1:
            lc_array1 = np.ravel(np.full(shape2, unlevel(array1), dtype = object))
            shape1 = shape2
        elif numel(array2) == 1:
            lc_array2 = np.ravel(np.full(shape1, unlevel(array2), dtype = object))
            shape2 = shape1
        else:
            raise ValueError("Input array shapes must currently be the same or one must only be length 1")

    output_array = []
    [output_array.append(i == j) for i, j in  zip(lc_array1, lc_array2)]
    
    return np.array(output_array).reshape((np.maximum(shape1, shape2)))

def isempty(array):
    if np.size(array) > 0:
        return False
    else:
        return True

def isequal(array1, array2):
    return np.array_equal(array1, array2)
        
def repmat(a, m, n = 1):
    return np.tile(a, (m, n))

def bsxfun(op_string, array1, array2):
    # if shapes of arrays don't match, explicitly expand them so they do
    if not (np.shape(array1) == np.shape(array2)):
        array1, array2 = broadcast_arrays(array1, array2)
        
    match op_string.lower():
        case 'plus':
            return np.array(array1) + np.array(array2)
        case 'minus':
            return np.array(array1) - np.array(array2)
        case 'times':
            return np.array(array1) * np.array(array2)
        case 'rdivide':
            return np.array(array1) / np.array(array2)
        case 'ldivide':
            return np.array(array2) / np.array(array1)
        case 'power':
            return np.array(array1) ** np.array(array2)
        case 'eq':
            return np.array(array1) == np.array(array2)
        case 'ne':
            return np.array(array1) != np.array(array2)
        case 'gt':
            return np.array(array1) > np.array(array2)
        case 'ge':
            return np.array(array1) >= np.array(array2)
        case 'lt':
            return np.array(array1) < np.array(array2)
        case 'le':
            return np.array(array1) <= np.array(array2)
        case 'and':
            return np.array(array1) & np.array(array2)
        case 'or':
            return np.array(array1) | np.array(array2)
        case 'xor':
            return np.logical_xor(np.array(array1), np.array(array2))
        case 'bitand':
            return np.bitwise_and(np.array(array1), np.array(array2))
        case 'bitor':
            return np.bitwise_or(np.array(array1), np.array(array2))
        case 'bitxor':
            return np.bitwise_xor(np.array(array1), np.array(array2))
        case 'max':
            return np.maximum(np.array(array1), np.array(array2))
        case 'min':
            return np.minimum(np.array(array1), np.array(array2))
        case 'mod':
            return mod(np.array(array1), np.array(array2))
        case 'rem':
            return rem(np.array(array1), np.array(array2))
        case 'atan2':
            return atan2(np.array(array1), np.array(array2))
        case 'atan2d':
            return np.degrees(atan2(np.array(array1), np.array(array2)))
        case 'hypot':
            return np.hypot(np.array(array1), np.array(array2))

def diff(array, n = 1, dim = None):
    if dim is None:
        if np.ndim(array) == 1:
            return np.diff(array, n = n)
        elif np.ndim(array) == 2:
            return np.diff(array, n = n, axis = 0)
    else:
        return np.diff(array, n = n, axis = dim - 1)

def length(array):
    # L = length(X) returns the length of the largest array dimension in X. For vectors, the length is simply the number of elements. For arrays with more dimensions, the length is max(size(X)). The length of an empty array is zero.
    if np.ndim(array) == 1:
        return np.size(array)
    else:
        return max(size(array))

def broadcast_arrays(array1, array2):
    if (np.size(array1) == 0) or (np.size(array2) == 0):
        raise ValueError("broadcast_arrays: Input arrays can't be empty!")

    if np.shape(array1) == () and np.shape(array2) != ():
        return np.full_like(array2, array1), array2
    elif np.shape(array2) == () and np.shape(array1) != ():
        return array1, np.full_like(array1, array2)
    elif np.shape(array1) == () and  np.shape(array2) == ():
        return array1, array2

    # make sure array have row and column dimensions
    if np.ndim(array1) == 1:
        array1 = array1.reshape(-1, 1)
    if np.ndim(array2) == 1:
        array2 = array2.reshape(1, -1)

    # transpose arrays (if required) so max dim is in second place
    array1Tflag = False
    array2Tflag = False
    if np.argmax(np.shape(array1)) == 0:
        array1 = array1.T
        array1Tflag = True
    if np.argmax(np.shape(array2)) == 0:
        array2 = array2.T
        array2Tflag = True

    max_dims = [np.shape(array1)[1], np.shape(array2)[1]]
    min_dims = [np.shape(array1)[0], np.shape(array2)[0]]
    if max_dims[0] == max_dims[1]:
        common_shape = (max(min_dims), max(max_dims))
    else:
        common_shape = (min(max_dims), max(max_dims))

    # print("Array1 shape: %s" % str(np.shape(array1)))
    # print("Array2 shape: %s" % str(np.shape(array2)))
    # print("Broadcast shape: %s" % str(common_shape))

    # try broadcast; transpose if failed
    try:
        array1 = np.broadcast_to(array1, common_shape)
    except:
        array1 = np.broadcast_to(array1.T, common_shape)
        if array1Tflag:
            array1Tflag = False
        else:
            array1Tflag = True
    try:
        array2 = np.broadcast_to(array2, common_shape)
    except:
        array2 = np.broadcast_to(array2.T, common_shape)
        if array2Tflag:
            array2Tflag = False
        else:
            array2Tflag = True

    # print("Array 1 transposed: %s" % str(array1Tflag))
    # print("Array 2 transposed: %s" % str(array1Tflag))

    if array1.shape != common_shape:
        array1.T
    if array2.shape != common_shape:
        array2.T

    if array1Tflag and array2Tflag:
        array1 = array1.T
        array2 = array2.T
        
    return array1, array2

def matlabmean(array, axis = 1):
    if axis == 0:
        raise ValueError("matlabmean function uses matlab 1 indexing")
    elif axis == 1:
        return np.mean(array, axis = axis - 1)
    elif axis == 2:
        return singcolvec(np.mean(array, axis = axis - 1))
    else:
        raise ValueError("matlabmean function only supports axis <= 2")

def matlabstd(array, axis = 1, ddof = 1):
    if axis == 0:
        raise ValueError("matlabstd function uses matlab 1 indexing")
    elif axis == 1:
        return np.std(array, axis = axis - 1, ddof = ddof)
    elif axis == 2:
        return singcolvec(np.std(array, axis = axis - 1, ddof = ddof))
    else:
        raise ValueError("matlabstd function only supports axis <= 2")

def matlabsum(array, axis = 1):
    if axis == 0:
        raise ValueError("matlabsum function uses matlab 1 indexing")
    elif axis == 1:
        return np.sum(array, axis = axis - 1)
    elif axis == 2:
        return singcolvec(np.sum(array, axis = axis - 1))
    else:
        raise ValueError("matlabsum function only supports axis <= 2")

def unique(array):
    return np.unique(array, return_index = True)

def gammaincinv(y, a):
    return sps.gammaincinv(a, y)

def gammainc(x, a):
    return sps.gammainc(a, x)

def besseli(nu, Z, scale = False):
    if not scale:
        return iv(nu, Z)
    else:
        return iv(nu, Z) * exp(-abs(real(Z)))

def is_single_column(arr):
    # Checks if a numpy array is a single column (1D or 2D).
    # Note: non-array single values are not counted as single columns; must be in array
    if (np.ndim(arr) == 1 and np.size(arr) == 1) or (np.ndim(arr) == 2 and np.shape(arr)[1] == 1):
        return True
    else:
        return False

def is_single_row(arr):
    # Checks if a numpy array is a single row (1D or 2D).
    # Note: non-array single values are not counted as single rows; must be in array
    if (np.ndim(arr) == 1) or (np.ndim(arr) == 2 and np.shape(arr)[0] == 1):
        return True
    else:
        return False

def find(X, n = None, direction = "first"):
    if direction.lower() == "first":
        if n is None:
            return np.nonzero(np.ravel(np.asarray(X).T))[0]
        else:
            return np.nonzero(np.ravel(np.asarray(X).T))[0][:n]
    elif direction.lower() == "last":
        if n is None:
            return np.nonzero(np.ravel(np.asarray(X).T))[0]
        else:
            return np.nonzero(np.ravel(np.asarray(X).T))[0][:n]
    else:
        raise ValueError("find function received invalid direction argument")

def sort(array, axis = 1):
    if axis == 0:
        raise ValueError("sort function uses matlab 1 indexing")
    if np.ndim(array) == 0:
        return array
    else:
        return np.sort(array, axis = axis - 1)

def magic(n):
    n = int(n)
    if n < 3:
        raise ValueError("Size must be at least 3")
    if n % 2 == 1:
        p = np.arange(1, n+1)
        return n*np.mod(p[:, None] + p - (n+3)//2, n) + np.mod(p[:, None] + 2*p-2, n) + 1
    elif n % 4 == 0:
        J = np.mod(np.arange(1, n+1), 4) // 2
        K = J[:, None] == J
        M = np.arange(1, n*n+1, n)[:, None] + np.arange(n)
        M[K] = n*n + 1 - M[K]
    else:
        p = n//2
        M = magic(p)
        M = np.block([[M, M+2*p*p], [M+3*p*p, M+p*p]])
        i = np.arange(p)
        k = (n-2)//4
        j = np.concatenate((np.arange(k), np.arange(n-k+1, n)))
        M[np.ix_(np.concatenate((i, i+p)), j)] = M[np.ix_(np.concatenate((i+p, i)), j)]
        M[np.ix_([k, k+p], [0, k])] = M[np.ix_([k+p, k], [0, k])]
    return M

def arrayfun(func, A, UniformOutput = True, ErrorHandler = None):
    """Applies the given function to the value."""
    return unlevel(list(map(func, A)))

def cellfun(func, A, UniformOutput = True, ErrorHandler = None):
    """Applies the given function to the value."""
    return unlevel(list(map(func, A)))

def quadgk(
    fun, # Integrand
    a, # Upper intregration limit
    b, # Upper intregration limit
    AbsTol = None, # Absolute error tolerance
    RelTol = None, # Relative error tolerance
    Waypoints = None, # Integration waypoints
    MaxIntervalCount = 650, # Maximum number of intervals allowed
):

    if AbsTol is None:
        AbsTol = 1e-10 # Matlab double precision tolerance; python is double precision by default

    if RelTol is None:
        RelTol = 1e-6 # Matlab double precision tolerance; python is double precision by default

    if (Waypoints is not None) and any(isinstance(x, complex) for x in Waypoints):
        raise ValueError("quadgk cannot yet handle complex waypoints")

    q, errbnd = quad(fun, a, b, limit = int(MaxIntervalCount), epsabs = AbsTol, epsrel = RelTol, points = Waypoints)

    return q, errbnd

def cumtrapz(var1, var2 = None, dim = 1):
    if dim == 0:
        raise ValueError("cumtrapz function uses matlab 1 indexing")
    if var2 is None:
        return spi.cumulative_trapezoid(var1, axis = dim - 1, initial = 0)
    else:
        return spi.cumulative_trapezoid(var2, var1, axis = dim - 1, initial = 0)

def intersect(set1, set2):
    return list(set(set1).intersection(set2))

def fix(array):
    return np.trunc(array).astype(int)

def matlabunion(set1, set2):
    return list(set(set1).union(set2))

def line(*args, **kwargs):
    # Make kwargs keys lowercase for easier matching
    normalized_kwargs = {k.lower(): v for k, v in kwargs.items()}

    # If a Matplotlib axis object is specified first, react accordingly
    if isinstance(args[0], matplotlib.axes.Axes):
        ax = args[0]
        x = args[1]
        y = args[2]
        if ax.name == "3d":
            z = args[3]
        else:
            z = None
    else:
        ax = None
        x = args[0]
        y = args[1]
        if isinstance(args[2], (list, np.ndarray)):
            z = args[2]
        else:
            z = None

    if "color" in normalized_kwargs.keys():
        color = normalized_kwargs["color"]
    else:
        color = "b"
        
    if "linestyle" in normalized_kwargs.keys():
        linestyle = normalized_kwargs["linestyle"]
    else:
        linestyle = "-"

    if "linewidth" in normalized_kwargs.keys():
        linewidth = normalized_kwargs["linewidth"]
    else:
        linewidth = 0.5
        
    if "marker" in normalized_kwargs.keys():
        marker = normalized_kwargs["marker"]
    else:
        marker = None
        
    if "markersize" in normalized_kwargs.keys():
        markersize = normalized_kwargs["markersize"]
    else:
        markersize = 6
            
    # Translate Matlab marker type into Matplotlib marker type, if present in keyword dict
    if marker is not None:
        match normalized_kwargs["marker"]:
            case "o": # Circle
                marker = "o"
            case "+": # Plus sign
                marker = "+"
            case "*": #	Asterisk
                marker = "*"
            case ".": #	Point
                marker = "."
            case "x": #	Cross
                marker = "x"
            case "_": #	Horizontal line
                marker = "_"
            case "|": #	Vertical line
                raise ValueError("line: that marker choice is not available in Matplotlib")
            case "square": # Square
                marker = "s"
            case "diamond": # Diamond
                marker = "d"
            case "^": #	Upward-pointing triangle
                marker = "^"
            case "v": #	Downward-pointing triangle
                marker = "v"
            case ">": #	Right-pointing triangle
                marker = ">"
            case "<": #	Left-pointing triangle
                marker = "<"
            case "pentagram": #	Pentagram
                raise ValueError("line: that marker choice is not available in Matplotlib")
            case "hexagram": # Hexagram
                raise ValueError("line: that marker choice is not available in Matplotlib")
            case "none": # No markers
                pass
            case None: # No markers
                pass
            case _:
                raise ValueError("line: that marker choice is not available in Matplotlib")

    if ax is None:
        if z is None:
            ax = plt.gca()
        else:
            if plt.gca().name != "3d":
                plt.clf()
                fig = plt.figure()
                ax = fig.add_subplot(projection='3d')

    if z is None:
        if marker is None:
            ax.plot(x, y, color = color, linestyle = linestyle, linewidth = linewidth)
        else:
            ax.plot(x, y, color = color, linestyle = linestyle, linewidth = linewidth, marker = marker, markersize = markersize)
    else:
        if marker is None:
            ax.plot(x, y, z, color = color, linestyle = linestyle, linewidth = linewidth)
        else:
            ax.plot(x, y, z, color = color, linestyle = linestyle, linewidth = linewidth, marker = marker, markersize = markersize)

def text(*args, **kwargs):
    # Make kwargs keys lowercase for easier matching
    normalized_kwargs = {k.lower(): v for k, v in kwargs.items()}

    # If a Matplotlib axis object is specified first, react accordingly
    if isinstance(args[0], matplotlib.axes.Axes):
        ax = args[0]
        x = args[1]
        y = args[2]
        if ax.name == "3d":
            z = args[3]
            txt = args[4]
        else:
            z = None
            txt = args[3]
    else:
        ax = None
        x = args[0]
        y = args[1]
        if isinstance(args[2], Number):
            z = args[2]
            txt = args[3]
        else:
            z = None
            txt = args[2]

    if "fontsize" in normalized_kwargs.keys():
        fontsize = normalized_kwargs["fontsize"]
    else:
        fontsize = 10 # maybe 9
        
    if "fontweight" in normalized_kwargs.keys():
        fontweight = normalized_kwargs["fontweight"]
    else:
        fontweight = "normal"

    if "fontname" in normalized_kwargs.keys():
        fontname = normalized_kwargs["fontname"]
    else:
        fontname = "Arial"
        
    if "color" in normalized_kwargs.keys():
        color = normalized_kwargs["color"]
    else:
        color = "b"
        
    if "horizontalalignment" in normalized_kwargs.keys():
        horizontalalignment = normalized_kwargs["horizontalalignment"]
    else:
        horizontalalignment = "left"

    if "verticalalignment" in normalized_kwargs.keys():
        verticalalignment = normalized_kwargs["verticalalignment"]
    else:
        verticalalignment = "baseline"

    if "position" in normalized_kwargs.keys():
        position = normalized_kwargs["position"]
    else:
        position = [0, 0]
        
    if "units" in normalized_kwargs.keys():
        raise ValueError("text: 'units' keyword not available")

    usetex = False
    if "interpreter" in normalized_kwargs.keys():
        if normalized_kwargs["interpreter"].lower() == "tex":
            usetex = True
        else:
            raise ValueError("text: only interpreter available is 'tex'")

    if ax is None:
        if z is None:
            ax = plt.gca()
        else:
            if plt.gca().name != "3d":
                plt.clf()
                fig = plt.figure()
                ax = fig.add_subplot(projection='3d')

    txt_array = []
    if z is None:
        for x_in, y_in, item in zip(x, y, txt):
            txt_array.append(ax.text(
                x_in,
                y_in,
                item,
                fontsize = fontsize,
                fontweight = fontweight,
                fontname = fontname,
                color = color,
                horizontalalignment = horizontalalignment,
                verticalalignment = verticalalignment,
                # position = position,
                usetex = usetex,
            ))
    else: # I don't believe Matplotlib will handle z coordinates
        for x_in, y_in, z_in, item in zip(inv_trans.transform((x, y, z)), txt):
            txt_array.append(ax.text(
                x_in,
                y_in,
                z_in,
                item,
                fontsize = fontsize,
                fontweight = fontweight,
                fontname = fontname,
                color = color,
                horizontalalignment = horizontalalignment, 
                verticalalignment = verticalalignment,
                # position = position,
                usetex = usetex,
            ))

    return txt_array

def Extent(text):
    fig = plt.gcf()
    ax = plt.gca()
    
    # You must draw the canvas to get a valid renderer and extent
    fig.canvas.draw()
    
    # Get the renderer instance
    renderer = fig.canvas.get_renderer()
    
    # Get the bounding box in window coordinates (pixels)
    bb_display = text.get_window_extent(renderer=renderer)
    
    # Transform the bounding box to data coordinates
    # This allows you to compare extents of different text objects in the same data space
    transf = ax.transData.inverted()
    bb_datacoords = bb_display.transformed(transf)
    
    return [bb_datacoords.x0, bb_datacoords.y0, bb_datacoords.x1, bb_datacoords.y1]

# def broadcast_sum(array1, array2):
#     if is_single_row(array1) and is_single_column(array2):
#         new_array1 = np.tile(array1, (np.size(array2), 1))
#         new_array2 = np.tile(array2, np.size(array1))
#         return new_array1 + new_array2
#     elif is_single_row(array2) and is_single_column(array1):
#         new_array2 = np.tile(array2, (np.size(array1), 1))
#         new_array1 = np.tile(array1, np.size(array2))
#         return new_array1 + new_array2
#     else:
#         raise ValueError("broadcast_sum requires one row vector and one column vector")

# def broadcast_sum(array1, array2, array3 = None, array4 = None):
#     array_of_arrays = [array1, array2]
    
#     # Check to make sure all arrays are either single columns or single rows
#     if not (is_single_row(array1) or is_single_column(array1)):
#         raise ValueError("broadcast_sum requires only row or column vectors")
#     elif not (is_single_row(array2) or is_single_column(array2)):
#         raise ValueError("broadcast_sum requires only row or column vectors")
#     elif array3 is not None:
#         if not (is_single_row(array3) or is_single_column(array3)):
#             raise ValueError("broadcast_sum requires only row or column vectors")
#         else:
#             array_of_arrays.append(array3)
#     elif array4 is not None:
#         if not (is_single_row(array4) or is_single_column(array4)):
#             raise ValueError("broadcast_sum requires only row or column vectors")
#         else:
#             array_of_arrays.append(array4)

#     # Check to make sure no vectors are already the same shape; if they are, add them
#     num_of_lists = sum(1 for array in array_of_arrays if isinstance(array, (list, np.ndarray)))
#     size_of_lists = [np.size(i) for i in array_of_arrays]
#     unique_lists, unique_idx = unique(size_of_lists)
#     size_of_unique_lists = [size_of_lists[i] for i in unique_idx]
#     new_array_of_arrays = []
#     if np.size(unique_idx) != num_of_lists:
#         missing_lists = [number for number in range(0, num_of_lists) if number not in unique_idx]
#         size_of_missing_lists = [np.size(i) for i in missing_lists]
#         for missing_list, size_of_missing_list in zip(missing_lists, size_of_missing_lists):
#             condition = [i == size_of_missing_list for i in size_of_lists]
#             new_array_of_arrays.append(sumif(array, condition))

# def sumif(array, condition, sum_range = None):
#     # sum array if single array is given, based on condition
#     # sum arrays element-wise if multiple arrays are given, based on condition
#     if sum_range is not None:
#         if isinstance(array, (list, np.ndarray)) and all(isinstance(item, (list, np.ndarray)) for item in array):        
#             return np.add(array[condition][sum_range])
#         else:
#             return sum(array[condition][sum_range])
#     else:
#         if isinstance(array, (list, np.ndarray)) and all(isinstance(item, (list, np.ndarray)) for item in array):
#             return np.add(*[i for idx, i in enumerate(array) if condition[idx]])
#         else:
#             return sum(array[condition])

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

def gf(array):
    return array

def CRCGenerator(name, value, msg):
    return calculate_crc(msg, value)

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