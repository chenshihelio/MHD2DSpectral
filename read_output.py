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
    return t, uu



#-------------------------------------


xgrid, ygrid = read_grid()
nx,ny= len(xgrid), len(ygrid)
print('nx, ny =', nx, ny)

XX,YY = np.meshgrid(xgrid,ygrid,indexing='ij')

kx = np.fft.rfftfreq(nx, xgrid[1] - xgrid[0])
ky = np.fft.rfftfreq(ny, ygrid[1] - ygrid[0])

npe,nvar = read_parallel_info()                                 
print('npe, nvar = ', npe, nvar)

t_EBM, radius, Ur = read_EBM()

files = sorted(glob.glob('./output/out*dat'))
nout = len(files)

time = np.zeros(nout)



for nt in range(nout):
    t, uu = read_uu(files[nt],nx,ny,nvar)

    print('t = {:.3f}'.format(t))

    time[nt] = t 

    rho = uu[:,:,0]

    ux = uu[:,:,1]
    uy = uu[:,:,2]
    uz = uu[:,:,3]

    Bx = uu[:,:,4]
    By = uu[:,:,5]
    Bz = uu[:,:,6]

    p = uu[:,:,7]

    Babs = np.sqrt(Bx**2 + By**2 + Bz**2)

    rho1 = rho - np.average(rho)
    p1 = p - np.average(p)
    Babs1 = Babs - np.average(Babs)

    
    # plot a 1D cut along x direction
    ind_y = int(ny/2)
    fig = plt.figure(figsize=[10,10])

    sub = fig.add_subplot(411)
    sub.plot(xgrid,rho1[:,ind_y],label=r'$\delta \rho$')
    sub.plot(xgrid,p1[:,ind_y],color='C1',label=r'$\delta P$')
    sub.plot(xgrid,Babs1[:,ind_y],color='k',label=r'$\delta |B|$')
    sub.legend()
    
    sub = fig.add_subplot(412,sharex=sub)
    sub.plot(xgrid,ux[:,ind_y],label=r'$u_x$')
    sub.plot(xgrid,Bx[:,ind_y],label=r'$b_x$')
    sub.legend()

    sub = fig.add_subplot(413,sharex=sub)
    sub.plot(xgrid,uy[:,ind_y],label=r'$u_y$')
    sub.plot(xgrid,By[:,ind_y],label=r'$b_y$')
    sub.legend()

    sub = fig.add_subplot(414,sharex=sub)
    sub.plot(xgrid,uz[:,ind_y],label=r'$u_z$')
    sub.plot(xgrid,Bz[:,ind_y],label=r'$b_z$')
    sub.legend()

    sub.set_xlabel(r'$x$')
    
    fig.suptitle(r'$t={:.2f}$'.format(t))
    fig.tight_layout(rect=[0,0,1,0.95])

    fig.savefig('./figure/cut_along_x_{:03d}.png'.format(nt))
    plt.close(fig)


    # plot a 1D cut along y direction
    ind_x = int(nx/2)
    fig = plt.figure(figsize=[10,10])

    sub = fig.add_subplot(411)
    sub.plot(ygrid,rho1[ind_x,:],label=r'$\delta \rho$')
    sub.plot(ygrid,p1[ind_x,:],color='C1',label=r'$\delta P$')
    sub.plot(ygrid,Babs1[ind_x,:],color='k',label=r'$\delta |B|$')
    sub.legend()

    sub = fig.add_subplot(412,sharex=sub)
    sub.plot(ygrid,ux[ind_x,:],label=r'$u_x$')
    sub.plot(ygrid,Bx[ind_x,:],label=r'$b_x$')
    sub.legend()

    sub = fig.add_subplot(413,sharex=sub)
    sub.plot(ygrid,uy[ind_x,:],label=r'$u_y$')
    sub.plot(ygrid,By[ind_x,:],label=r'$b_y$')
    sub.legend()

    sub = fig.add_subplot(414,sharex=sub)
    sub.plot(ygrid,uz[ind_x,:],label=r'$u_z$')
    sub.plot(ygrid,Bz[ind_x,:],label=r'$b_z$')
    sub.legend()

    fig.savefig('./figure/cut_along_y_{:03d}.png'.format(nt))
    plt.close(fig)
    

    # plot 2D color-coded maps
    figsize=[8,7]
    # rho---
    vlim = [] # [0.75,1.2]
    if_div_cmap = True
    fig = plt.figure(figsize=figsize)
    sub = fig.add_subplot(111)
    if if_div_cmap:
        pm = sub.pcolormesh(XX,YY,rho1,shading='gouraud',cmap='RdBu',norm=TwoSlopeNorm(vcenter=0))
    else:
        if len(vlim)!=2:
            pm = sub.pcolormesh(XX,YY,rho,shading='gouraud',cmap=plt.cm.plasma)
        else:
            pm = sub.pcolormesh(XX,YY,rho,shading='gouraud',cmap=plt.cm.plasma,
                vmin=vlim[0],vmax=vlim[1])
    
    sub.set_aspect('equal')
    cb = fig.colorbar(pm,shrink=0.9,aspect=15)
    sub.set_xlabel(r'$x$',fontsize=12)
    sub.set_ylabel(r'$y$',fontsize=12)
    sub.set_title(r'$\delta \rho$'+' $t={:.2f}$'.format(t))
    fig.savefig('./figure/rho1_{:03d}.png'.format(nt))
    plt.close(fig)


    # P---
    vlim = []# [0.15,0.3]
    if_div_cmap = True
    fig = plt.figure(figsize=figsize)
    sub = fig.add_subplot(111)
    if if_div_cmap:
        pm = sub.pcolormesh(XX,YY,p1,shading='gouraud',cmap='RdBu',norm=TwoSlopeNorm(vcenter=0))
    else:
        if len(vlim)!=2:
            pm = sub.pcolormesh(XX,YY,p,shading='gouraud',cmap=plt.cm.plasma)
        else:
            pm = sub.pcolormesh(XX,YY,p,shading='gouraud',cmap=plt.cm.plasma,
                vmin=vlim[0],vmax=vlim[1])
    
    sub.set_aspect('equal')
    cb = fig.colorbar(pm,shrink=0.9,aspect=15)
    sub.set_xlabel(r'$x$',fontsize=12)
    sub.set_ylabel(r'$y$',fontsize=12)
    sub.set_title(r'$\delta P$'+' $t={:.2f}$'.format(t))
    fig.savefig('./figure/P1_{:03d}.png'.format(nt))
    plt.close(fig)


    # Babs---
    vlim = [] # [0.75,1.2]
    if_div_cmap = True
    fig = plt.figure(figsize=figsize)
    sub = fig.add_subplot(111)
    if if_div_cmap:
        pm = sub.pcolormesh(XX,YY,Babs1,shading='gouraud',cmap='RdBu',norm=TwoSlopeNorm(vcenter=0))
    else:
        if len(vlim)!=2:
            pm = sub.pcolormesh(XX,YY,Babs,shading='gouraud',cmap=plt.cm.plasma)
        else:
            pm = sub.pcolormesh(XX,YY,Babs,shading='gouraud',cmap=plt.cm.plasma,
                vmin=vlim[0],vmax=vlim[1])
    
    sub.set_aspect('equal')
    cb = fig.colorbar(pm,shrink=0.9,aspect=15)
    sub.set_xlabel(r'$x$',fontsize=12)
    sub.set_ylabel(r'$y$',fontsize=12)
    sub.set_title(r'$\delta |B|$'+' $t={:.2f}$'.format(t))
    fig.savefig('./figure/Babs1_{:03d}.png'.format(nt))
    plt.close(fig)



    # Bz---
    vlim = [] # [0.15,0.3]
    fig = plt.figure(figsize=figsize)
    sub = fig.add_subplot(111)
    if len(vlim)!=2:
        pm = sub.pcolormesh(XX,YY,Bz,shading='gouraud',cmap=plt.cm.RdBu,norm=TwoSlopeNorm(vcenter=0))
    else:
        pm = sub.pcolormesh(XX,YY,Bz,shading='gouraud',cmap=plt.cm.RdBu,
            vmin=vlim[0],vmax=vlim[1])
    
    sub.set_aspect('equal')
    cb = fig.colorbar(pm,shrink=0.9,aspect=15)
    sub.set_xlabel(r'$x$',fontsize=12)
    sub.set_ylabel(r'$y$',fontsize=12)
    sub.set_title(r'$B_z$'+' $t={:.2f}$'.format(t))
    fig.savefig('./figure/Bz_{:03d}.png'.format(nt))
    plt.close(fig)