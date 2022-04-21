"""Provide functions used to estimate coexistence points"""
from ctypes import c_char
import numpy as np
from scipy import optimize

def delta_f(f1new: float,f2new: float,f: np.ndarray) -> np.ndarray:
    """
    Calculate the difference between next and current integration points

    Parameters
    ----------
    f1new : float
        The next integration point f1 in the coexistence line
    
    f2new : float
        The next integration point f2 in the coexistence line
    
    f: np.ndarray [2,iphase]
        Current integration points f1 and f2
    Returns
    -------
    df : [2,iphase]
        Difference between next and current integration points
    """
    df = np.zeros((2,2))
    df[0,:] = f1new - f[0,:]
    df[1,:] = f2new - f[1,:]
    return df

def poly_coefficients(df: np.ndarray,z: np.ndarray,cov: np.ndarray) -> np.ndarray:
    """
    Calculate the coefficients in the free energy polynomial 

    Parameters
    ----------
    df : [2,iphase]
        Difference between next and current integration points
    
    z: np.ndarray [2,iphase]
        Conjugate varibales (z1,z2) of currrent point (f1,f2) for both I and II phases

    cov: np.ndarray [3,iphase]
        Covariances [cov(z1,Z1),cov(z2,Z2),cov(z1,Z2)] of current point for both I and II phases

    Returns
    -------
    df : [6,2]
        Coefficients in the free energy polynomial
    """
    coef = np.zeros((6,2))
    coef[0,:] = z[0,:]*df[0,:]
    coef[1,:] = z[1,:]*df[1,:]
    coef[2,:] = cov[0,:]*df[0,:]**2
    coef[3,:] = cov[1,:]*df[1,:]**2
    coef[4,:] = cov[2,:]*df[0,:]*df[1,:]
    coef[5,:] = cov[2,:]*df[0,:]
    return coef

def first_guess_newton(free_energy: np.ndarray,z: np.ndarray,f1new: float,f: np.ndarray)->float:
    """
    Calculate the first guess for the new point f2 in the coexistence line

    Parameters
    ----------

    free_energy: np.ndarray [iphase]
        Free energy of phases I and II in the current point 
    
    z: np.ndarray [2,iphase]
        Conjugate varibales (z1,z2) of currrent point (f1,f2) for both I and II phases

    f1new : float
        The next integration point f1 in the coexistence line
    
    f: np.ndarray [2,iphase]
        Current integration points f1 and f2
    
    Returns
    -------
    f2new0 : float
        First guess for the new point f2 in the coexistence line
    
    """
    f20_a = free_energy[0] - free_energy[1] +z[0,0]*(f1new-f[0,0]) - z[0,1]*(f1new-f[0,1])
    f20_b = -z[1,0]*f[1,0]+z[1,1]*f[1,1]
    f2new0 = (f20_a+f20_b)/(z[1,1] - z[1,0])

    return f2new0

def f_first_point(f2new: float,free_energy:np.ndarray,z:np.ndarray,cov:np.ndarray,f1new:float,f:np.ndarray)->float:
    """
    Free energy difference function of the next point to be optimized if the Number of points is equal to 1
    or one is refining the estimation

    Parameters
    ----------

    f1new : float
        The next integration point f1 in the coexistence line
    
    f2new : float
        The next integration point f2 in the coexistence line
    
    free_energy: np.ndarray [iphase]
        Free energy of phases I and II in the current point 
    
    z: np.ndarray [2,iphase]
        Conjugate varibales (z1,z2) of currrent point (f1,f2) for both I and II phases
    
    cov: np.ndarray [3,iphase]
        Covariances [cov(z1,Z1),cov(z2,Z2),cov(z1,Z2)] of current point for both I and II phases
    
    f: np.ndarray [2,iphase]
        Current integration points f1 and f2

     Returns
    -------
    fun[1] - fun[2] : float
        Free energy difference in the next point   
    """
    fun = np.zeros(2)    
    df = delta_f(f1new,f2new,f)
    coef=poly_coefficients(df,z,cov)
    fun   = (free_energy+coef[0,:] +
             coef[1,:]-0.5*(coef[2,:] +
             coef[3,:]) - coef[4,:])
    return fun[0] - fun[1] 

def df_first_point(f2new: float,free_energy:np.ndarray,z: np.ndarray,cov: np.ndarray,f1new: float,f: np.ndarray) -> float:

    """
    Calculate the partial derivative of Free energy difference function (f_first_point) with respect 
    to the new coexistence point (f2)
    
    Parameters
    ----------

    f1new : float
        The next integration point f1 in the coexistence line
    
    f2new : float
        The next integration point f2 in the coexistence line
    
    z: np.ndarray [2,iphase]
        Conjugate varibales (z1,z2) of currrent point (f1,f2) for both I and II phases
    
    cov: np.ndarray [3,iphase]
        Covariances [cov(z1,Z1),cov(z2,Z2),cov(z1,Z2)] of current point for both I and II phases
    
    f: np.ndarray [2,iphase]
        Current integration points f1 and f2

     Returns
    -------
    dfun[1] - dfun[2] : float
        Partial derivative of Free energy difference function  
    """
    dfun = np.zeros(2)
    df = delta_f(f1new,f2new,f)
    dfun = z[1,:] -cov[1,:]*df[1,:] - cov[2,:]*df[0,:]
    return dfun[0] - dfun[1]

def refine_coex(int_1) -> 'tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray]':
    """
    Refines the near coexistence data

    Parameters
    ----------
    Attributes of class int_1:
    f1new : float
        The next integration point f1 in the coexistence line
    
    f2new : float
        The next integration point f2 in the coexistence line
    
    z: np.ndarray [2,iphase]
        Conjugate varibales (z1,z2) of currrent point (f1,f2) for both I and II phases
    
    cov: np.ndarray [3,iphase]
        Covariances [cov(z1,Z1),cov(z2,Z2),cov(z1,Z2)] of current point for both I and II phases
    
    f: np.ndarray [2,iphase]
        Current integration points f1 and f2
    ene: np.ndarray [iphase] stats[1,:,:]
        Average potential energy for both I and II phases

     Returns
    -------
    zsat : np.ndarray
        Refined conjugate varibales (z1,z2) of currrent point (f1,f2) for both I and II phases
    f2sat : np.ndarray
        Refined point (f2) for both I and II phases
    enesat : np.ndarray
        Estimation of saturated average potential energy for both I and II phases
    """
    ene = int_1.stats[1,:,:]
    Npoints = np.shape(int_1.f)[1]
    enesat = np.zeros((Npoints,2))
    zsat = np.zeros((2,Npoints,2))
    f2sat = np.zeros(Npoints)
    free_energy_zsat = np.zeros((Npoints,2))
    for i in range(Npoints):
        f2new0 = int_1.f[1,i,0] 
        f1sat  = int_1.f[0,i,0]
        f2sat[i] = optimize.newton(f_first_point, f2new0,fprime=df_first_point,
                              args=(int_1.free_energy[i,:],int_1.z[:,i,:],int_1.cov[:,i,:],f1sat,int_1.f[:,i,:]),
                              maxiter=500)
        
        df = delta_f(f1sat,f2sat[i],int_1.f[:,i,:])
        zsat[0,i,:] = int_1.z[0,i,:] -int_1.cov[0,i,:]*df[0,:] - int_1.cov[2,i,:]*df[1,:]
        zsat[1,i,:] = int_1.z[1,i,:] -int_1.cov[1,i,:]*df[1,:] - int_1.cov[2,i,:]*df[0,:]
        coef=poly_coefficients(df,zsat[:,i,:],int_1.cov[:,i,:])
        #redundancy 
     
        free_energy_zsat[i,:]   = (int_1.free_energy[i,:]+coef[0,:] +
                 coef[1,:]-0.5*(coef[2,:] +
                 coef[3,:]) - coef[4,:])
        
        enesat[i,:] = ene[i,:] - int_1.cov[2,i,:]*df[1,:]
        
        sat_prop = dict()
        sat_prop["z_sat"] =zsat
        sat_prop["free_energy_sat"] = free_energy_zsat
        sat_prop["ene_sat"] = enesat
        sat_prop["f2_sat"] = f2sat
    return sat_prop#zsat,free_energy_zsat,enesat,f2sat

def cal_free_energy(f_a: np.ndarray,f_b: np.ndarray,z_a: np.ndarray,z_b: np.ndarray,cov_a: np.ndarray,cov_b: np.ndarray,free_energy_a: np.ndarray) -> np.ndarray:
    """
    Calculate the free-energy differences between the previous
    and current stats simulated for each phase, and hence estimate
    the free energy of current state

    Parameters
    ----------

    f_b : np.ndarray
        Current integration points f1 and f2
    
    f_b : np.ndarray
        Previous integration points f1 and f2
    
    z_b: np.ndarray [2,iphase]
        Conjugate varibales (z1,z2) of currrent point (f1,f2) for both I and II phases
    
    z_a: np.ndarray [2,iphase]
        Conjugate varibales (z1,z2) of previous point (f1,f2) for both I and II phases
    
    cov_b: np.ndarray [3,iphase]
        Covariances [cov(z1,Z1),cov(z2,Z2),cov(z1,Z2)] of current point for both I and II phases

    cov_a: np.ndarray [3,iphase]
        Covariances [cov(z1,Z1),cov(z2,Z2),cov(z1,Z2)] of previous point for both I and II phases
    
    free_energy_a: np.ndarray [iphase]
        Free energy of phases I and II in the previous point

     Returns
    -------
    free_energy_b : np.ndarray
        Free energy of phases I and II in the current point


    """
    df = delta_f(f_b[0,:],f_b[1,:],f_a)
    coef_a=poly_coefficients(df,z_a,cov_a)
    coef_b=poly_coefficients(df,z_b,cov_b)
    free_energy_b = (free_energy_a+0.5*(coef_a[0,:]+coef_b[0,:]+
    coef_a[1,:]+coef_b[1,:])+(coef_b[2,:]-coef_a[2,:]+coef_b[3,:]-coef_a[3,:])/12
    +(coef_b[4,:]-coef_a[4,:])/6) 
    
    return free_energy_b

def calculate_first_point(f1new: float,f: np.ndarray,free_energy: np.ndarray, z: np.ndarray, cov: np.ndarray) -> float:
    """
    Calculate the first next point in the integration in the Nf1f2 ensemble or 
    Refines the  next point in the integration in the Nf1f2 ensemble
    Parameters
    ----------
    f1new : float
        The next integration point f1 in the coexistence line
    
    f: np.ndarray [2,iphase]
        Current integration points f1 and f2
    
    free_energy: np.ndarray [iphase]
        Free energy of the phases  
    
    z: np.ndarray [2,iphase]
        Conjugate varibales (z1,z2) of previous point (f1,f2) for both I and II phases
    
    cov: np.ndarray [3,iphase]
        Covariances [cov(z1,Z1),cov(z2,Z2),cov(z1,Z2)] of previous point for both I and II phases


    Returns
    -------
    f2new : float
        The next property 2 point in the integration or
        Refined current point(f2)
    """

    f2new0 = first_guess_newton(free_energy[0,:],z[:,0,:],f1new,f[:,0,:])
    f2new = optimize.newton(f_first_point, f2new0,fprime=df_first_point,
                          args=(free_energy[0,:],z[:,0,:],cov[:,0,:],f1new,f[:,0,:]),
                          maxiter=500)
    return f2new

def f_next_point(f2new: float,free_energyb: np.ndarray,za: np.ndarray,zb: np.ndarray,covb: np.ndarray,f1new: float,fa: np.ndarray,fb: np.ndarray) -> float:

    """
    Free energy difference function of the next point to be optimized if the Number of points is greater than 1

    REF : J. Chem. Phys. 146, 134508 (2017)


    Parameters
    ----------

    f1new : float
        The next integration point f1 in the coexistence line
    
    f2new : float
        The next integration point f2 in the coexistence line
    
    free_energyb: np.ndarray [iphase]
        Free energy of phases I and II in the current point 
    
    zb: np.ndarray [2,iphase]
        Conjugate varibales (z1,z2) of currrent point (f1,f2) for both I and II phases
    
    za: np.ndarray [2,iphase]
        Conjugate varibales (z1,z2) of previous point (f1,f2) for both I and II phases
    
    covb: np.ndarray [3,iphase]
        Covariances [cov(z1,Z1),cov(z2,Z2),cov(z1,Z2)] of current point for both I and II phases
    
    f_b : np.ndarray [2,iphase]
        Current integration points f1 and f2
    
    f_b : np.ndarray [2,iphase]
        Previous integration points f1 and f2

     Returns
    -------
    fun[1] - fun[2] : float
        Free energy difference in the next point   
    """

    fun = np.zeros(2)
    dfab = delta_f(fb[0,:],fb[1,:],fa)
    
    omega1 = (za[0,:]-zb[0,:]-covb[0,:]*dfab[0]
              -covb[2,:]*dfab[1])/(3*dfab[0]**2)
    omega2 = (za[1,:]-zb[1,:]-covb[1,:]*dfab[1]
              -covb[2,:]*dfab[0])/(3*dfab[1]**2)
    df = delta_f(f1new,f2new,fb)
    coef=poly_coefficients(df,zb,covb)
    fun   = (free_energyb+coef[0,:] +
             coef[1,:]-0.5*(coef[2,:] +
             coef[3,:]) - coef[4,:] + 
             omega1*df[0]**3+omega2*df[1]**3)
    return fun[0] - fun[1]

def df_next_point(f2new: float,free_energyb: np.ndarray,za: np.ndarray,zb: np.ndarray,covb: np.ndarray,f1new: float,fa: np.ndarray,fb: np.ndarray) -> float:
    """
    Calculate the partial derivative of Free energy difference function (f_next_point) with respect 
    to the new coexistence point (f2)


    Parameters
    ----------

    f1new : float
        The next integration point f1 in the coexistence line
    
    f2new : float
        The next integration point f2 in the coexistence line
    
    free_energyb: np.ndarray [iphase]
        Free energy of phases I and II in the current point 
    
    zb: np.ndarray [2,iphase]
        Conjugate varibales (z1,z2) of currrent point (f1,f2) for both I and II phases
    
    za: np.ndarray [2,iphase]
        Conjugate varibales (z1,z2) of previous point (f1,f2) for both I and II phases
    
    f_b : np.ndarray [2,iphase]
        Current integration points f1 and f2
    
    f_b : np.ndarray [2,iphase]
        Previous integration points f1 and f2

     Returns
    -------
    fun[1] - fun[2] : float
        Partial derivative of Free energy difference function  
    """
    fun = np.zeros(2)
    dfab = delta_f(fb[0,:],fb[1,:],fa)
    omega2 = (za[1,:]-zb[1,:]-covb[1,:]*dfab[1]
              -covb[2,:]*dfab[0])/(3*dfab[1]**2)
    df = delta_f(f1new,f2new,fb)
    fun   = ( zb[1,:] - covb[1,:]*df[1] -
             covb[2,:]*df[1] + 3*omega2*df[1]**2)
    return fun[0] - fun[1]

def f_next_point_zeta(f2new: float,free_energyb: np.ndarray,za: np.ndarray,zb: np.ndarray,cova: np.ndarray,covb: np.ndarray,f1new: float,fa: np.ndarray,fb: np.ndarray) -> float:

    """
    Free energy difference function of the next point to be optimized if the Number of points is greater than 1
    Specialized for the case that when f1 can not be linearly uncoupled from the free energy 

    The polynomial form of the free energy is different from that of f_next_point with
    n = 3 in that an asymmetry in the f 1 and f 2 terms is introduced to forsake the need to evaluate second (or higher) order
    derivatives of free energy with respect to f 1. 
    REF : Escobedo, F. J. Chem. Phys. 147, 214501 (2017)

    Parameters
    ----------

    f1new : float
        The next integration point f1 in the coexistence line
    
    f2new : float
        The next integration point f2 in the coexistence line
    
    free_energyb: np.ndarray [iphase]
        Free energy of phases I and II in the current point 
    
    zb: np.ndarray [2,iphase]
        Conjugate varibales (z1,z2) of currrent point (f1,f2) for both I and II phases
    
    za: np.ndarray [2,iphase]
        Conjugate varibales (z1,z2) of previous point (f1,f2) for both I and II phases
    
    covb: np.ndarray [3,iphase]
        Covariances [cov(z1,Z1),cov(z2,Z2),cov(z1,Z2)] of current point for both I and II phases

    cova: np.ndarray [3,iphase]
        Covariances [cov(z1,Z1),cov(z2,Z2),cov(z1,Z2)] of previous point for both I and II phases
    
    f_b : np.ndarray [2,iphase]
        Current integration points f1 and f2
    
    f_b : np.ndarray [2,iphase]
        Previous integration points f1 and f2

     Returns
    -------
    fun[1] - fun[2] : float
        Free energy difference in the next point   
    """
    fun = np.zeros(2)
    dfab = delta_f(fb[0,:],fb[1,:],fa)
    omega1 = (-covb[1,:] + cova[1,:])/(6*dfab[1])
    omega2 = (zb[1,:]-za[1,:] + 0.5*(cova[1,:]+covb[1,:])*dfab[1]
              )/dfab[0]
    omega3 = (zb[0,:]-za[0,:])/(2*dfab[0]) - dfab[1]*omega2/(2*dfab[0])

    df = delta_f(f1new,f2new,fb)
    coef=poly_coefficients(df,zb,covb)
    fun   = (free_energyb+coef[0,:] +
             coef[1,:]-0.5*coef[3,:] + omega3*df[0]**2
             + omega1*df[1]**3+omega2*df[0]*df[1])
    return fun[0] - fun[1]

def df_next_point_zeta(f2new: float,free_energyb: np.ndarray,za: np.ndarray,zb: np.ndarray,cova: np.ndarray,covb: np.ndarray,f1new: float,fa: np.ndarray,fb: np.ndarray) -> float:
    """
    Calculate the partial derivative of Free energy difference function (f_next_point_zeta) with respect 
    to the new coexistence point (f2)


    Parameters
    ----------

    f1new : float
        The next integration point f1 in the coexistence line
    
    f2new : float
        The next integration point f2 in the coexistence line
    
    free_energyb: np.ndarray [iphase]
        Free energy of phases I and II in the current point 
    
    zb: np.ndarray [2,iphase]
        Conjugate varibales (z1,z2) of currrent point (f1,f2) for both I and II phases
    
    za: np.ndarray [2,iphase]
        Conjugate varibales (z1,z2) of previous point (f1,f2) for both I and II phases
    
    f_b : np.ndarray [2,iphase]
        Current integration points f1 and f2
    
    f_b : np.ndarray [2,iphase]
        Previous integration points f1 and f2

     Returns
    -------
    fun[1] - fun[2] : float
        Partial derivative of Free energy difference function  
    """
    fun = np.zeros(2)
    dfab = delta_f(fb[0,:],fb[1,:],fa)
    omega1 = (-covb[0,:] + covb[1,:])/(6*dfab[1])
    omega2 = (zb[1,:]-za[1,:] + (cova[1,:]+covb[1,:])*dfab[1]
              )/(2*dfab[0])
    df = delta_f(f1new,f2new,fb)
    fun   = (zb[1,:]-covb[1,:]*df[1] +
             + 3*omega1*df[1]**2+omega2*df[0])
    return fun[0] - fun[1]

def estimate_coexistence(int_type ,int_1) -> 'tuple[np.ndarray,float]' :
    """
    Calculate the  next point (f2) in the integration and its free energies in the Nf1f2 ensemble 
    Parameters
    ----------
    int_type:character
        Integration type: defines the function used to calculate the next point
        Options are coupled or decoupled
    Attributes of class int_1:
    f1new : float
        The next integration point f1 in the coexistence line
    
    f: np.ndarray [2,npoints,iphase]
        Integration points f1 and f2
    
    free_energy: np.ndarray [npoints,iphase]
        Free energy of the initial integartion point  [:,0]
        the rest of the array is zeros
    
    z: np.ndarray [2,npoints,iphase]
        Conjugate varibales (z1,z2) of previous point (f1,f2) for both I and II phases
    
    cov: np.ndarray [3,npoints,iphase]
        Covariances [cov(z1,Z1),cov(z2,Z2),cov(z1,Z2)] of previous point for both I and II phases


    Returns
    -------
    f2new : float
        The next property 2 point to carry out a simulation in the Nf1f2 ensemble 

    free_energy: np.ndarry [Npoints,iphase]
        Free energies for all the points in the intagetion
    """
    
    
    Npoints = np.shape(int_1.f)[1]
    
    if int(Npoints)==1:
        f2new=calculate_first_point(int_1.f1new,int_1.f,int_1.free_energy, int_1.z, int_1.cov)
    else:
        for ipt in range(Npoints-1):
            bpt = ipt +1
            int_1.free_energy[bpt,:] = (cal_free_energy(int_1.f[:,ipt,:],int_1.f[:,bpt,:],
            int_1.z[:,ipt,:],int_1.z[:,bpt,:],int_1.cov[:,ipt,:],int_1.cov[:,bpt,:],int_1.free_energy[ipt,:]))
        apt = Npoints - 2
        bpt = Npoints - 1
        f2new0 = first_guess_newton(int_1.free_energy[bpt,:],int_1.z[:,bpt,:],int_1.f1new,int_1.f[:,bpt,:])
    
        if (int_type=='coupled'):
            f2new = optimize.newton(f_next_point_zeta,f2new0,fprime=df_next_point_zeta,
                              args=(int_1.free_energy[bpt,:],int_1.z[:,apt,:],int_1.z[:,bpt,:],
                              int_1.cov[:,apt,:],int_1.cov[:,bpt,:],int_1.f1new,int_1.f[:,apt,:],int_1.f[:,bpt,:],),
                              maxiter=500)
        elif (int_type=='decoupled'):
            f2new = optimize.newton(f_next_point,f2new0,fprime=df_next_point,
                              args=(int_1.free_energy[bpt,:],int_1.z[:,apt,:],int_1.z[:,bpt,:],
                              int_1.cov[:,bpt,:],int_1.f1new,int_1.f[:,apt,:],int_1.f[:,bpt,:],),
                              maxiter=500)
        else:
            return print('Integration types are coupled or decoupled')
    return int_1.free_energy,f2new



