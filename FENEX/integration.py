
from queue import Empty
import numpy as np
from FENEX import cal_free_energy,calculate_first_point,first_guess_newton,f_next_point_zeta,df_next_point_zeta,f_next_point,df_next_point
from FENEX import f_first_point,df_first_point,delta_f,poly_coefficients
import sys
from scipy import optimize
class Integrate:
    def __init__(self, Npoints, f1new, f, free_energy, z,cov,stats,int_type):
      self.free_energy = free_energy
      self.z = z
      self.Npoints = Npoints
      self.f1new = f1new
      self.f = f
      self.cov = cov
      self.stats =stats
      self.int_type = int_type
    def estimate_coexistence(self):
      """

         Calculate the  next point (f2) in the integration and its free energies in the Nf1f2 ensemble 
      Parameters
      ----------
      int_type:character
          Integration type: defines the function used to calculate the next point
          Options are coupled or decoupled
      Attributes of class self:
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


      Npoints = np.shape(self.f)[1]

      if int(Npoints)==1:
          self.f2new=calculate_first_point(self.f1new,self.f,self.free_energy, self.z, self.cov)
      else:
          for ipt in range(Npoints-1):
              bpt = ipt +1
              self.free_energy[bpt,:] = (cal_free_energy(self.f[:,ipt,:],self.f[:,bpt,:],
              self.z[:,ipt,:],self.z[:,bpt,:],self.cov[:,ipt,:],self.cov[:,bpt,:],self.free_energy[ipt,:]))
          apt = Npoints - 2
          bpt = Npoints - 1
          f2new0 = first_guess_newton(self.free_energy[bpt,:],self.z[:,bpt,:],self.f1new,self.f[:,bpt,:])

          if (self.int_type=='coupled'):
              self.f2new = optimize.newton(f_next_point_zeta,f2new0,fprime=df_next_point_zeta,
                                args=(self.free_energy[bpt,:],self.z[:,apt,:],self.z[:,bpt,:],
                                self.cov[:,apt,:],self.cov[:,bpt,:],self.f1new,self.f[:,apt,:],self.f[:,bpt,:],),
                                maxiter=500)
          elif (self.int_type=='decoupled'):
              self.f2new = optimize.newton(f_next_point,f2new0,fprime=df_next_point,
                                args=(self.free_energy[bpt,:],self.z[:,apt,:],self.z[:,bpt,:],
                                self.cov[:,bpt,:],self.f1new,self.f[:,apt,:],self.f[:,bpt,:],),
                                maxiter=500)
          else:
              print('Integration types are coupled or decoupled')
              sys.exit(1)
      
    def refine_coexistence(self):
      """
      Refines the near coexistence data

      Parameters
      ----------
      Attributes of class self:
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
      free_energy_sat : np.ndarray
          Refined Free energies for all the points in the intagetion
      """
      ene = self.stats[1,:,:]
      Npoints = np.shape(self.f)[1]
      self.enesat = np.zeros((Npoints,2))
      self.zsat = np.zeros((2,Npoints,2))
      self.f2sat = np.zeros(Npoints)
      self.free_energy_sat = np.zeros((Npoints,2))
      for i in range(Npoints):
          f2new0 = self.f[1,i,0] 
          f1sat  = self.f[0,i,0]
          self.f2sat[i] = optimize.newton(f_first_point, f2new0,fprime=df_first_point,
                                args=(self.free_energy[i,:],self.z[:,i,:],self.cov[:,i,:],f1sat,self.f[:,i,:]),
                                maxiter=500)

          df = delta_f(f1sat,self.f2sat[i],self.f[:,i,:])
          self.zsat[0,i,:] = self.z[0,i,:] -self.cov[0,i,:]*df[0,:] - self.cov[2,i,:]*df[1,:]
          self.zsat[1,i,:] = self.z[1,i,:] -self.cov[1,i,:]*df[1,:] - self.cov[2,i,:]*df[0,:]
          coef=poly_coefficients(df,self.zsat[:,i,:],self.cov[:,i,:])
          #redundancy 

          self.free_energy_sat[i,:]   = (self.free_energy[i,:]+coef[0,:] +
                   coef[1,:]-0.5*(coef[2,:] +
                   coef[3,:]) - coef[4,:])

          self.enesat[i,:] = ene[i,:] - self.cov[2,i,:]*df[1,:]