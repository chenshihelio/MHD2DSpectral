import numpy as np
import struct
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.gridspec import GridSpec
import glob
import sys
from matplotlib.colors import TwoSlopeNorm,SymLogNorm


def read_parallel_info(filename='./output/parallel_info.dat'):
    # read parallel information ----------------
    file_prl = open(filename,'rb')

    skip = struct.unpack("f",file_prl.read(4))

    npe = int((struct.unpack("f",file_prl.read(4)))[0])
    nvar = int((struct.unpack("f",file_prl.read(4)))[0])

    skip = struct.unpack("f",file_prl.read(4))

    file_prl.close()
    return npe, nvar


def read_grid(filename='./output/grid.dat'):
    # read nx, ny, nz and grid-------------------------
    file_grid = open(filename, 'rb')

    skip = struct.unpack("f",file_grid.read(4))

    nx = int((struct.unpack("f",file_grid.read(4)))[0])
    ny = int((struct.unpack("f",file_grid.read(4)))[0])

    skip = struct.unpack("f",file_grid.read(4))

    xgrid = np.zeros(nx)
    ygrid = np.zeros(ny)

    skip = struct.unpack("f",file_grid.read(4))

    for i in range(nx):
        xgrid[i] = (struct.unpack("f",file_grid.read(4)))[0]
    for i in range(ny):
        ygrid[i] = (struct.unpack("f",file_grid.read(4)))[0]

    skip = struct.unpack("f",file_grid.read(4))

    file_grid.close()
    return xgrid,ygrid

def read_EBM(filename='./output/EBM_info.dat'):
    #read the EBM_info.dat
    #which includes time,radius,Ur

    file_EBM = np.array(np.loadtxt(filename))

    if len(file_EBM.shape)==1:
        file_EBM = np.reshape(file_EBM,(int(len(file_EBM)/3),3))

    t_EBM = file_EBM[:,0]
    radius = file_EBM[:,1]
    Ur_EBM = file_EBM[:,2]

    return t_EBM,radius,Ur_EBM

def read_uu(filename,nx,ny,nvar):
    file_uu = open(filename, 'rb')
    uu = np.zeros([nx,ny,nvar])

    skip = (struct.unpack("f",file_uu.read(4)))[0]
    t = (struct.unpack("f",file_uu.read(4)))[0] 
    skip = (struct.unpack("f",file_uu.read(4)))[0]

    for ivar in range(nvar):
            for iy in range(ny):
                for ix in range(nx):
                    uu[ix,iy,ivar] = (struct.unpack(\
                        "d",file_uu.read(8)))[0] 

    file_uu.close()
    
    return t, uu



def calc_deriv(arr,x,dim=0):
    # dim = 0(x) or 1(y)
    # arr: 2D array
    # x: 1D array

    nx = len(x)
    if arr.shape[dim]!=nx:
        print('Shape of arr does not match len(x)!!')
        return None


    dx = x[1] - x[0]
    kx = np.fft.fftfreq(nx,d=dx)

    # fft along the 
    arr_fftx = np.fft.fft(arr,axis=dim)
    kxarr = np.zeros(arr_fftx.shape,dtype=complex)
    
    for i in range(arr.shape[(dim+1)%2]):
        if dim==0:
            kxarr[:,i] = 2*np.pi* 1j * kx * arr_fftx[:,i]
        elif dim==1:
            kxarr[i,:] = 2*np.pi* 1j * kx * arr_fftx[i,:]


    darr_dx = np.real(np.fft.ifft(kxarr,axis=dim))
    arr_fftx = None 
    kxarr = None

    return darr_dx



def current_density(Bx,By,Bz,x,y):
    
    # calculate derivatives
    # dBx_dy
    dBx_dy = calc_deriv(Bx,y,dim=1)

    # dBx_dz = 0
    dBx_dz = 0

    # dBy_dx
    dBy_dx = calc_deriv(By,x,dim=0)

    # dBy_dz = 0
    dBy_dz = 0

    # dBz_dx
    dBz_dx = calc_deriv(Bz,x,dim=0)

    # dBz_dy 
    dBz_dy = calc_deriv(Bz,y,dim=1)

    # calculate current density
    Jx = dBz_dy - dBy_dz 
    Jy = dBx_dz - dBz_dx 
    Jz = dBy_dx - dBx_dy

    return [Jx,Jy,Jz]




def calc_VpVg_fastandslow(Cs,Ca,theta):
    cos_angle = np.cos(theta)
    sin_angle = np.sin(theta)

    Cm = np.sqrt(Cs**2 + Ca**2)

    Cn4 = Cm**4 - 4 * Cs**2 * Ca**2 * cos_angle**2
    if Cn4<0:
        Cn4 = 0

    Cn2 = np.sqrt(Cn4)

    tmp2 = Cm**2 - Cn2
    if tmp2<0:
        tmp2 = 0

    Cp_slow = np.sqrt(0.5 * (tmp2) )
    Cp_fast = np.sqrt(0.5 * (Cm**2 + Cn2) )

    Cg_slow_perp = sin_angle * Cp_slow * (1 - Cs**2 * Ca**2 / Cp_slow**2 
        / Cn2 * cos_angle**2)
    Cg_slow_para = cos_angle * Cp_slow * (1 + Cs**2 * Ca**2 / Cp_slow**2 
        / Cn2 * sin_angle**2)

    Cg_fast_perp = sin_angle * Cp_fast * (1 + Cs**2 * Ca**2 / Cp_fast**2 
        / Cn2 * cos_angle**2)
    Cg_fast_para = cos_angle * Cp_fast * (1 - Cs**2 * Ca**2 / Cp_fast**2 
        / Cn2 * sin_angle**2)

    return {'Vps':Cp_slow,'Vpf':Cp_fast,
        'Vgs_perp':Cg_slow_perp,'Vgs_para':Cg_slow_para,
        'Vgf_perp':Cg_fast_perp,'Vgf_para':Cg_fast_para}

        

def func_for_maximum_Vgsperp_costheta(h,cs = 0.5, ca = 1.0):
    cs2 = cs*cs
    ca2 = ca*ca 
    cm2 = cs2 + ca2 
    cm4 = cm2 * cm2 
    cn4 = cm4 - 4*cs2*ca2*h
    if cn4<0:
        cn4 = 0
    
    cn2 = np.sqrt(cn4)

    vp2 = 0.5 * (cm2 - cn2)


    r = cs2 * ca2 / vp2 / cn2 

    return r * (1-h) * (1-r*h) - (1-r*h) - 2 * (1-h)\
        *(r - h * r*r * ( 2 -  cm2 / cn2))



#-------------------------------------
def main():
    print('You are running the functions module... Nothing will happen...')
    input('Press any key to exit...')
    return 1

if __name__ == '__main__':
    main()