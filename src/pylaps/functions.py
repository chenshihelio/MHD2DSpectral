import numpy as np
from numpy import cos, sin, sqrt

def calc_deriv(arr, x, dim=0):
    # dim = 0(x) or 1(y)
    # arr: 2D array
    # x: 1D array

    nx = len(x)
    if arr.shape[dim] != nx:
        print("Shape of arr does not match len(x)!!")
        return None

    dx = x[1] - x[0]
    kx = np.fft.fftfreq(nx, d=dx)

    # fft along the
    arr_fftx = np.fft.fft(arr, axis=dim)
    kxarr = np.zeros(arr_fftx.shape, dtype=complex)

    for i in range(arr.shape[(dim + 1) % 2]):
        if dim == 0:
            kxarr[:, i] = 2 * np.pi * 1j * kx * arr_fftx[:, i]
        elif dim == 1:
            kxarr[i, :] = 2 * np.pi * 1j * kx * arr_fftx[i, :]

    darr_dx = np.real(np.fft.ifft(kxarr, axis=dim))

    return darr_dx


def current_density(Bx, By, Bz, x, y):

    # calculate derivatives
    # dBx_dy
    dBx_dy = calc_deriv(Bx, y, dim=1)

    # dBx_dz = 0
    dBx_dz = 0

    # dBy_dx
    dBy_dx = calc_deriv(By, x, dim=0)

    # dBy_dz = 0
    dBy_dz = 0

    # dBz_dx
    dBz_dx = calc_deriv(Bz, x, dim=0)

    # dBz_dy
    dBz_dy = calc_deriv(Bz, y, dim=1)

    # calculate current density
    Jx = dBz_dy - dBy_dz
    Jy = dBx_dz - dBz_dx
    Jz = dBy_dx - dBx_dy

    return [Jx, Jy, Jz]

def VpVg_fs(cs, ca, theta=None):
    """Calculates the phase velocities and group velocities for fast and slow waves."""
    if theta is None:
        theta = np.arange(0, 360, 1) * np.pi / 180

    cos_angle = cos(theta)
    sin_angle = sin(theta)

    cm = sqrt(cs**2 + ca**2)

    cn4 = cm**4 - 4 * cs**2 * ca**2 * cos_angle**2
    cn2 = sqrt(cn4)

    tmp2 = cm**2 - cn2

    Cp_slow = sqrt(0.5 * (tmp2))
    Cp_fast = sqrt(0.5 * (cm**2 + cn2))

    Cg_slow_perp = (
        sin_angle * Cp_slow * (1 - cs**2 * ca**2 / Cp_slow**2 / cn2 * cos_angle**2)
    )
    Cg_slow_para = (
        cos_angle * Cp_slow * (1 + cs**2 * ca**2 / Cp_slow**2 / cn2 * sin_angle**2)
    )

    Cg_fast_perp = (
        sin_angle * Cp_fast * (1 + cs**2 * ca**2 / Cp_fast**2 / cn2 * cos_angle**2)
    )
    Cg_fast_para = (
        cos_angle * Cp_fast * (1 - cs**2 * ca**2 / Cp_fast**2 / cn2 * sin_angle**2)
    )

    return {
        "Vps": Cp_slow,
        "Vpf": Cp_fast,
        "Vgs_perp": Cg_slow_perp,
        "Vgs_para": Cg_slow_para,
        "Vgf_perp": Cg_fast_perp,
        "Vgf_para": Cg_fast_para,
    }