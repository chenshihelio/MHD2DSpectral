import numpy as np
import struct

def read_parallel_info(filename="./diags/parallel_info.dat"):
    """Read parallel information from the given binary file.
    """
    with open(filename, 'rb') as file:
        file.read(4)
        npe = int(struct.unpack("f", file.read(4))[0])
        nvar = int(struct.unpack("f", file.read(4))[0])
    return npe, nvar

def read_EBM(filename="./diags/EBM_info.dat"):
    # read the EBM_info.dat
    # which includes time,radius,Ur

    file_EBM = np.array(np.loadtxt(filename))

    if len(file_EBM.shape) == 1:
        file_EBM = np.reshape(file_EBM, (int(len(file_EBM) / 3), 3))

    t_EBM = file_EBM[:, 0]
    radius = file_EBM[:, 1]
    Ur_EBM = file_EBM[:, 2]

    return t_EBM, radius, Ur_EBM

def read_grid(filename='./diags/grid.dat'):
    """
    Read grid data from the given binary file.

    Parameters:
    filename (str): Path to the grid data file.

    Returns:
    np.ndarray: xgrid values.
    np.ndarray: ygrid values.
    """
    with open(filename, 'rb') as file:
        file.seek(4)
        nx = int(struct.unpack("f", file.read(4))[0])
        ny = int(struct.unpack("f", file.read(4))[0])
        file.read(8)
        xgrid = np.fromfile(file, dtype='float32', count=nx)
        ygrid = np.fromfile(file, dtype='float32', count=ny)
    return xgrid, ygrid

def read_uu(filename, nx, ny, nvar):
    """
    Read uu data from the given binary file.

    Parameters:
    filename (str): Path to the uu data file.
    nx (int): Number of x grid points.
    ny (int): Number of y grid points.
    nvar (int): Number of variables.

    Returns:
    float: Time value.
    np.ndarray: uu data array.
    """
    with open(filename, 'rb') as file_uu:
        file_uu.seek(4)
        t = struct.unpack("f", file_uu.read(4))[0]
        file_uu.read(4)
        
        uu = np.fromfile(file_uu, dtype='float64', count=nvar*ny*nx).reshape((nvar, ny, nx))
        #change the order of the array
        uu = np.moveaxis(uu, [0, 1, 2], [2, 1, 0])
    return t, uu