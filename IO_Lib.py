import numpy as np

#Use mask and apply it over all arrays
def get_clean_coord(X,Y,Z,timestep, xlim, ylim, zlim):
    i = timestep
    X_clean = np.ma.masked_where(X[:,i]>xlim[1], X[:,i]) 
    Y_clean = np.ma.masked_where(X[:,i]>xlim[1], Y[:,i])
    Z_clean = np.ma.masked_where(X[:,i]>xlim[1], Z[:,i])
    
    X_clean = np.ma.masked_where(X_clean[:]<xlim[0], X_clean[:]) 
    Y_clean = np.ma.masked_where(X_clean[:]<xlim[0], Y_clean[:])
    Z_clean = np.ma.masked_where(X_clean[:]<xlim[0], Z_clean[:])
    
    X_clean = np.ma.masked_where(Y_clean[:]>ylim[1], X_clean[:]) 
    Y_clean = np.ma.masked_where(Y_clean[:]>ylim[1], Y_clean[:]) 
    Z_clean = np.ma.masked_where(Y_clean[:]>ylim[1], Z_clean[:])
    
    X_clean = np.ma.masked_where(Y_clean[:]<ylim[0], X_clean[:]) 
    Y_clean = np.ma.masked_where(Y_clean[:]<ylim[0], Y_clean[:]) 
    Z_clean = np.ma.masked_where(Y_clean[:]<ylim[0], Z_clean[:])
    
    
    X_clean = np.ma.masked_where(Z_clean[:]>zlim[1], X_clean[:]) 
    Y_clean = np.ma.masked_where(Z_clean[:]>zlim[1], Y_clean[:])
    Z_clean = np.ma.masked_where(Z_clean[:]>zlim[1], Z_clean[:])
    
    X_clean = np.ma.masked_where(Z_clean[:]<zlim[0], X_clean[:]) 
    Y_clean = np.ma.masked_where(Z_clean[:]<zlim[0], Y_clean[:])
    Z_clean = np.ma.masked_where(Z_clean[:]<zlim[0], Z_clean[:])
    
    return X_clean, Y_clean, Z_clean


# Get coordinates of atom having name : atomname from file named : filename and store in a local txt file

def gen_txt(path, filename, atomname):
    pdb = open(path+filename+'.pdb')
    trimmed_pdb = open('Coord_'+filename+'.txt',"w")
    for line in pdb:
        if line[0] == 'A' and line[13] == atomname:
            trimmed_pdb.write(line[31:45]+' '+line[47:55] +'\n')
    pdb.close()
    trimmed_pdb.close()
    
def get_coordinates(filename, N_molecules, N_timesteps):
    trimmed_pdb = open('Coord_'+filename+'.txt')
    x = []
    for line in trimmed_pdb:
        x.append(line)
    
    X = np.zeros((N_molecules,N_timesteps))
    Y = np.zeros((N_molecules,N_timesteps))
    Z = np.zeros((N_molecules,N_timesteps))
    for j in range(0,N_molecules*N_timesteps):
        c = [float(i) for i in x[j].split()]
        i1 = j%N_molecules
        j1 = int(j/N_molecules)
        X[i1][j1] = c[0]
        Y[i1][j1] = c[1]
        Z[i1][j1] = c[2]
        
    return X,Y,Z
