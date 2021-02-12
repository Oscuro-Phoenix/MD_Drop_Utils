from numba import jit
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

@jit
def get_rho(cell_size,xlim,ylim,zlim,x0,y0,z0,dx,dz):
    
    p = cell_size/2
    q = 0
    Nx = int(1/dx)
    Nz = int(1/dz)
    #Populate Density Matrix for each cell
    rho = np.zeros(((xlim[1]-xlim[0])*Nx,(zlim[1]-zlim[0])*Nz))

    for x in range(0,(xlim[1]-xlim[0])*Nx):
        for z in range(0,(zlim[1]-zlim[0])*Nz):
            y_temp = y0[(x0>x/Nx+xlim[0]-p) & (x0<x/Nx+xlim[0]+p) & (z0>z/Nz+zlim[0]-p) & (z0<z/Nz+zlim[0]+p)]
            if y_temp.size != 0 and y_temp.max() != y_temp.min():
                rho[x][z] = y_temp.shape[0]/(y_temp.max()-y_temp.min())
    return rho


def plot_snapshot(X, Y, Z, timestep):
    plt.scatter(X[:,timestep],Z[:,timestep])


def plot_rho(rho,xlim,zlim,dx,dz):
    q = 0
    Nx = int(1/dx)
    Nz = int(1/dz)
    z_ax = np.linspace((zlim[0])*Nz, (zlim[1]-1)*Nz,  (zlim[1]-zlim[0])*Nz)
    x_ax = np.linspace((xlim[0])*Nx, (xlim[1]-1)*Nx,  (xlim[1]-xlim[0])*Nx)
    x_ax, z_ax = np.meshgrid(x_ax,z_ax)
    fig = plt.figure(figsize=(20,10))
    ax = plt.axes(projection='3d')
    surf = ax.plot_surface(z_ax, x_ax, rho.T,cmap='viridis')
    ax.set_xlabel('Z')
    ax.set_ylabel('X')
    ax.set_title(r'Level Plot of $ \rho $')
    fig.colorbar(surf)
    ax.view_init(-90, 0)
    fig
