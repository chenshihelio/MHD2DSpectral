import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.gridspec import GridSpec
import glob
import sys
from matplotlib.colors import TwoSlopeNorm,SymLogNorm

sys.path.append('../')
from functions import *


#-------------------------------------
xgrid, ygrid = read_grid()
nx,ny= len(xgrid), len(ygrid)
print('nx, ny =', nx, ny)

XX,YY = np.meshgrid(xgrid,ygrid,indexing='ij')

kx = np.fft.rfftfreq(nx, xgrid[1] - xgrid[0])
ky = np.fft.rfftfreq(ny, ygrid[1] - ygrid[0])

npe,nvar = read_parallel_info()                                 
print('npe, nvar = ', npe, nvar)

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



    

    # plot 2D color-coded maps
    figsize=[15,5]
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(2,3)

    sub = fig.add_subplot(gs[0,0])
    pm = sub.pcolormesh(XX,YY,rho,shading='gouraud',cmap=plt.cm.plasma)
    cb = fig.colorbar(pm,shrink=0.9,aspect=15)
    sub.set_xlabel(r'$x$',fontsize=12)
    sub.set_ylabel(r'$y$',fontsize=12)
    sub.set_title(r'$\rho$'+' $t={:.2f}$'.format(t))

    sub = fig.add_subplot(gs[1,0])
    pm = sub.pcolormesh(XX,YY,p,shading='gouraud',cmap=plt.cm.plasma)
    cb = fig.colorbar(pm,shrink=0.9,aspect=15)
    sub.set_xlabel(r'$x$',fontsize=12)
    sub.set_ylabel(r'$y$',fontsize=12)
    sub.set_title(r'$P$'+' $t={:.2f}$'.format(t))

    sub = fig.add_subplot(gs[0,1])
    pm = sub.pcolormesh(XX,YY,ux,shading='gouraud',cmap=plt.cm.RdBu_r,
        vmin=-np.max(np.abs(ux)),vmax=np.max(np.abs(ux)))
    cb = fig.colorbar(pm,shrink=0.9,aspect=15)
    sub.set_xlabel(r'$x$',fontsize=12)
    sub.set_ylabel(r'$y$',fontsize=12)
    sub.set_title(r'$U_x$'+' $t={:.2f}$'.format(t))

    sub = fig.add_subplot(gs[1,1])
    pm = sub.pcolormesh(XX,YY,uy,shading='gouraud',cmap=plt.cm.RdBu_r,
        vmin=-np.max(np.abs(uy)),vmax=np.max(np.abs(uy)))
    cb = fig.colorbar(pm,shrink=0.9,aspect=15)
    sub.set_xlabel(r'$x$',fontsize=12)
    sub.set_ylabel(r'$y$',fontsize=12)
    sub.set_title(r'$U_y$'+' $t={:.2f}$'.format(t))

    sub = fig.add_subplot(gs[0,2])
    pm = sub.pcolormesh(XX,YY,Bx,shading='gouraud',cmap=plt.cm.RdBu_r,
        vmin=-np.max(np.abs(Bx)),vmax=np.max(np.abs(Bx)))
    cb = fig.colorbar(pm,shrink=0.9,aspect=15)
    sub.set_xlabel(r'$x$',fontsize=12)
    sub.set_ylabel(r'$y$',fontsize=12)
    sub.set_title(r'$B_x$'+' $t={:.2f}$'.format(t))

    sub = fig.add_subplot(gs[1,2])
    pm = sub.pcolormesh(XX,YY,By,shading='gouraud',cmap=plt.cm.RdBu_r,
        vmin=-np.max(np.abs(By)),vmax=np.max(np.abs(By)))
    cb = fig.colorbar(pm,shrink=0.9,aspect=15)
    sub.set_xlabel(r'$x$',fontsize=12)
    sub.set_ylabel(r'$y$',fontsize=12)
    sub.set_title(r'$B_y$'+' $t={:.2f}$'.format(t))

    gs.tight_layout(fig)

    fig.savefig('./figure/{:03d}.png'.format(nt))
    plt.close(fig)
