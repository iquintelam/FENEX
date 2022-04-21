
import numpy as np
class SimulationData:
    def __init__(self, Npoints, f1new, f, free_energy, z,cov,stats):
      self.free_energy = free_energy
      self.z = z
      self.Npoints = Npoints
      self.f1new = f1new
      self.f = f
      self.cov = cov
      self.stats =stats
