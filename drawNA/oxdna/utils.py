import numpy as np
import scipy.linalg as la

# Debesh's oxDNA constants
SHIFT_BASE = 0.5
SHIFT_ACROSS = 0.56 - (SHIFT_BASE * 0.9)
SHIFT_ROUND = -0.105

def get_rotation_matrix(axis : np.ndarray, theta : float) -> np.ndarray:
    """Returns a rotation matrix for rotating an angle theta
    (rad) around an axis vector
    """
    return la.expm(np.cross(np.eye(3), axis/la.norm(axis)*theta))
    
def next_5p(
    position: np.ndarray,
    a1: np.ndarray, 
    a3: np.ndarray,
    angle: float,
    rise: float
):
    """Returns the pos_com for a new nucleotide in the 5' direction
    """
    # shift round in a1 direction
    rotation_matrix = get_rotation_matrix(a3, angle)
    a1 = np.dot(
            rotation_matrix,
            a1,
        )
    
    # shift up in a3 direction
    new_base = position + rise * a3
    new_pos = new_base - a1 * SHIFT_BASE
    new_pos += SHIFT_ROUND * np.cross(a3, a1)
    return new_pos

def get_box() -> np.ndarray:
    pass

def round_to_multiple(n, mo=0.34, decimal_places=2):
    """
    Function rounds to the nearest multiple of value given
    Returns output to (default) 2 decimal places

    Arguments:

        n --- value (integer or float) to round  
        mo --- "multiple_of" is the value (integer or float) which 
        we want to round to a multiple of  
        decimal_places --- no. of decimals to return
    """
    a = (n // mo) * mo  # Smaller multiple
    b = a + mo  # Larger multiple
    closest_multiple = b if n - a > b - n else a  # Return of closest of two
    return round(closest_multiple, decimal_places)

# the three commented lines would define quantities that are not necessary
def quat_to_exyz (myquat):
    sqw = myquat[0] * myquat[0];
    sqx = myquat[1] * myquat[1];
    sqy = myquat[2] * myquat[2];
    sqz = myquat[3] * myquat[3];

    invs = 1 / (sqx + sqy + sqz + sqw)
    m00 = (sqx - sqy - sqz + sqw) * invs ;
    #m11 = (-sqx + sqy - sqz + sqw) * invs ;
    m22 = (-sqx - sqy + sqz + sqw) * invs ;
    
    tmp1 = myquat[1] * myquat[2];
    tmp2 = myquat[3] * myquat[0];
    m10 = 2.0 * (tmp1 + tmp2) * invs ;
    #m01 = 2.0 * (tmp1 - tmp2) * invs ;

    tmp1 = myquat[1] * myquat[3];
    tmp2 = myquat[2] * myquat[0];
    m20 = 2.0 * (tmp1 - tmp2) * invs ;
    m02 = 2.0 * (tmp1 + tmp2) * invs ;
    tmp1 = myquat[2] * myquat[3];
    tmp2 = myquat[1] * myquat[0];
    #m21 = 2.0 * (tmp1 + tmp2) * invs ;
    m12 = 2.0 * (tmp1 - tmp2) * invs ; 

    mya1 = np.array([m00, m10, m20])
    mya3 = np.array([m02, m12, m22])

    return mya1, mya3