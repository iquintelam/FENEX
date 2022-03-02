from functions import *
import pandas as pd
import numpy as np
import linecache

def read_covdata_MC(filename: str) -> 'tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray]':
    lines_in_file = open(filename, 'r').readlines()
    nbins = len(lines_in_file) -4
    f = np.zeros(2)
    z = np.zeros(2)
    cov = np.zeros(3)
    stats = np.zeros(2)
    _,_,f[0],f[1] = open(filename, 'r').readline().split()
    for i in range(2):
        z[i] , stats[i] = linecache.getline(filename,nbins+i).split()
    cov[0]  =  linecache.getline(filename,nbins+3)
    cov[1]  =  linecache.getline(filename,nbins+2)
    cov[2]  =  linecache.getline(filename,nbins+4)
    return f,z,cov,stats

def read_input_MC(filename: str):
    
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
                f[:,i,j],z[:,i,j],cov[:,i,j],stats[:,i,j] = read_covdata_MC(filephase[0])
        f1new=float(file.readline().split()[0])
        
    return Npoints,f1new,f,free_energy, z, cov,stats


