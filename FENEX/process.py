# import pandas as pd
import numpy as np
import os
def read_covdata_MC(filename: str) -> 'tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray]':
    """
    Reads the covariance data from out of a simulation in a spacific format
    ----------
    filename:str
        File name of the covariance data
    Returns
    -------
    f1new : float
        The next integration point f1 in the coexistence line
    
    f: np.ndarray [2]
        Integration points f1 and f2
    
    free_energy: np.ndarray [iphase,npoints]
        Free energy of the initial integartion point  [:,0]
        the rest of the array is zeros
    
    z: np.ndarray [2]
        Conjugate varibales (z1,z2) 
    
    cov: np.ndarray [3]
        Covariances [cov(z1,Z1),cov(z2,Z2),cov(z1,Z2)] 
    
    stats: np.ndarray [2]
        Additional stats of the simulation
        stats[0] : acceptance probability of the perturbation
        stats[1] : average potential energy
    """
    lines_in_file = open(filename, 'r').readlines()

    f = np.zeros(2)
    z = np.zeros(2)
    cov = np.zeros(3)
    stats = np.zeros(2)

    _, _, f[0], f[1] = lines_in_file[0].split()

    important_info = lines_in_file[-5:]

    for i in range(2):
        z[i], stats[i] = ( float(x) for x in important_info[i].split() )

    cov[0] = float(important_info[3])
    cov[1] = float(important_info[2])
    cov[2] = float(important_info[4])

    return f,z,cov,stats

def read_test_system() -> 'tuple[int,float,np.ndarray,np.ndarray,np.ndarray,np.ndarray,np.ndarray]':
    """
    Reads input file of FENEX method
    ----------
    filename:str
        File name of the covariance data
    Returns
    -------
    Npoints : int
        Number of simulaterd integartion points

    f1new : float
        The next integration point f1 in the coexistence line
    
    f: np.ndarray [2,Npoints,2]
        Integration points f1 and f2
    
    free_energy: np.ndarray [npoints,iphase]
        Free energy of the initial integartion point  [:,0]
        the rest of the array is zeros
    
    z: np.ndarray [2,Npoints,2]
        Conjugate varibales (z1,z2) 
    
    cov: np.ndarray [3,Npoints,2]
        Covariances [cov(z1,Z1),cov(z2,Z2),cov(z1,Z2)] 
    
    stats: np.ndarray [2,Npoint,2]
        Additional stats of the simulation
        stats[0] : acceptance probability of the perturbation
        stats[1] : average potential energy
    """

    dir_path = os.path.dirname(os.path.realpath(__file__))
    filename = os.path.join(dir_path, "data", "simulation_data.dat")

    with open(filename, 'r') as file:
        fe = file.readline().split()[:2]
        Npoints = int(file.readline().split()[:1][0])
        free_energy = np.zeros((Npoints,2))
        f = np.zeros((2,Npoints,2))
        z = np.zeros((2,Npoints,2))
        stats = np.zeros((2,Npoints,2))
        cov = np.zeros((3,Npoints,2))
        free_energy[0,:] = fe
        for i in range(Npoints):
            for j in range(2):
                filephase=file.readline().split()
                # Recommend finding another way to find filepath - still depends on script
                # being executed in a particular place
                f[:,i,j],z[:,i,j],cov[:,i,j],stats[:,i,j] = read_covdata_MC(filephase[0])
        f1new=float(file.readline().split()[0])
        
    return Npoints,f1new,f,free_energy, z, cov,stats


