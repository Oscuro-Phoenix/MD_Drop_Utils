from numba import jit
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

@jit
def get_rho(cell_size,xlim,ylim,zlim,x0,y0,z0,padding,dx,dz):
    
    p = cell_size/2
    q = padding
    Nx = int(1/dx)
    Nz = int(1/dz)
    #Populate Density Matrix for each cell
    rho = np.zeros(((xlim[1]-xlim[0]+2*q)*Nx,(zlim[1]-zlim[0]+2*q)*Nz))

    for x in range(-q*Nx,(xlim[1]-xlim[0]+q)*Nx):
        for z in range(-q*Nz,(zlim[1]-zlim[0]+q)*Nz):
            y_temp = y0[(x0>x/Nx+xlim[0]-p) & (x0<x/Nx+xlim[0]+p) & (z0>z/Nz+zlim[0]-p) & (z0<z/Nz+zlim[0]+p)]
            if y_temp.size != 0 and y_temp.max() != y_temp.min():
                rho[x+q*Nx][z+q*Nz] = y_temp.shape[0]/(y_temp.max()-y_temp.min())
    return rho


def plot_snapshot(X, Y, Z, timestep):
    plt.scatter(X[:,timestep],Z[:,timestep])


def plot_rho(rho,padding,xlim,zlim,dx,dz):
    q = padding
    Nx = int(1/dx)
    Nz = int(1/dz)
    z_ax = np.linspace((zlim[0]-q)*Nz, (zlim[1]+q-1)*Nz,  (zlim[1]+2*q-zlim[0])*Nz)
    x_ax = np.linspace((xlim[0]-q)*Nx, (xlim[1]+q-1)*Nx,  (xlim[1]+2*q-xlim[0])*Nx)
    x_ax, z_ax = np.meshgrid(x_ax,z_ax)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    surf = ax.plot_surface(z_ax, x_ax, rho.T,cmap='viridis')
    ax.set_xlabel('Z')
    ax.set_ylabel('X')
    ax.set_title(r'Level Plot of $ \rho $')
    fig.colorbar(surf)
    ax.view_init(-90, 0)
    fig
