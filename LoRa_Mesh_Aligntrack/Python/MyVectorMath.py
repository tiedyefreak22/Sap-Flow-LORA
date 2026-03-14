import numpy as np
from numpy.typing import NDArray

def rotate_vector(vector: NDArray, gamma: float, beta: float, alpha: float):
    """
    Intrinsically rotates a 3D vector.

    Args:
        vector (np.array): The N x 3 input vector e.g. [x, y, z]
        gamma  (float): roll angle of rotation in degrees, about x-axis
        beta   (float): pitch angle of rotation in degrees, about y-axis
        alpha  (float): yaw angle of rotation in degrees, about z-axis

    Returns:
        np.array: The rotated vector.
    """
    yaw = np.array(
        [[cos(np.radians(alpha)), -sin(np.radians(alpha)), 0],
         [sin(np.radians(alpha)),  cos(np.radians(alpha)), 0],
         [0,                       0,                      1]]
    )
    pitch = np.array(
        [[cos(np.radians(beta)),   0,  sin(np.radians(beta))],
         [0,                       1,                      0],
         [-sin(np.radians(beta)),  0,  cos(np.radians(beta))]]
    )
    roll = np.array(
        [[1,                       0,                      0],
         [0, cos(np.radians(gamma)), -sin(np.radians(gamma))],
         [0, sin(np.radians(gamma)),  cos(np.radians(gamma))]]
    )
    
    R = yaw @ pitch @ roll # Rotation matrix

    # np.dot handles transposition of vector to columns implicitly
    return np.dot(R, vector)

def translate_vector(vector: NDArray, vx: float, vy: float, vz: float):
    """
    Translates a 3D vector.

    Args:
        vector (np.array): The N x 3 input vector e.g. [x, y, z]
        vx     (float): x translation value
        vy     (float): y translation value
        vz     (float): z translation value

    Returns:
        np.array: The translated vector.
    """
    # Create a column of ones (N, 1)
    ones = ones((vector.shape[0], 1))
    
    # Append the column to the points array to get (N, 4) homogeneous coordinates
    vector_4d_homogeneous = np.hstack((vector, ones))

    T = np.array( # Translation matrix
        [[1, 0, 0, vx],
         [0, 1, 0, vy],
         [0, 0, 1, vz],
         [0, 0, 0,  1]]
    )

    new_vector = np.dot(T, vector_4d_homogeneous.T).T
    return np.delete(new_vector, -1, axis=1)

def mag_vec(vector: NDArray):
    """
    Gets magnitude of vector.

    Args:
        vector (np.array): N-element input vector e.g. [x, y, z]

    Returns:
        float: The magnitude of the vector.
    """
    running_sum = 0
    for component_idx in range(np.size(vector)):
        running_sum += vector[component_idx] ** 2
    return sqrt(running_sum)

def norm_vec(vector: NDArray):
    """
    Normalizes vector by it's magnitude.

    Args:
        vector (np.array): N-element input vector e.g. [x, y, z]

    Returns:
        np.array: The normalized vector.
    """
    return vector / mag_vec(vector)