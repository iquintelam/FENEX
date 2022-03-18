"""
Unit and regression test for the FENEX package.
"""

# Import package, test suite, and other packages as needed
import sys
import os
from numpy.testing import assert_almost_equal

import pytest


from FENEX import *
# sys.path.append('../')
# from FENEX.Data_Proc import *
# from FENEX.functions import *


import matplotlib.pyplot as plt

def test_FENEX_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "FENEX" in sys.modules

#f,z,cov,stats = read_covdata_MC(r"data/iso.dat")
def test_first_point():
    dir_path = os.path.dirname(os.path.realpath(__file__))
    data1 = os.path.join(dir_path, "..", "data", "simulation_data_one.dat")
    Npoints,f1new,f,free_energy, z, cov,stats = read_input_MC(data1)
    f2new = calculate_first_point(f1new,f,free_energy, z, cov) 
    f2check = 8.27367550491069
    delta = 0.0001
    # assert function() to check if values are almost equal
    assert_almost_equal(f2new, f2check,decimal=delta)

def test_fenex():

    # Change file path
    dir_path = os.path.dirname(os.path.realpath(__file__))
    data1 = os.path.join(dir_path, "..", "data", "simulation_data.dat")
    #data1 = r"data\simulation_data.dat"
    Npoints,f1new,f,free_energy, z, cov,stats = read_input_MC(data1)

    free_energy,f2new = calculate_next_point(1 ,f1new,f,free_energy, z, cov) 
    f2check=8.37242894118646
    zsat,free_energy_zsat,enesat,f2sat = calc_zsat(Npoints,free_energy,z,cov,f,stats[1,:,:])
    delta =6e-4
    f2sat_check = [8.1685 ,8.2736]
    zsat_check =[[[-.2095241E+01 ,-.1595087E+01], [ -.1967657E+01 ,-.1498863E+01]],
 [[ 0.2225709E+01 ,0.2463211E+01], [ 0.2225398E+01 ,0.2459797E+01]]]
    # assert function() to check if values are almost equal
    assert_almost_equal(f2new, f2check,decimal=delta)
    assert_almost_equal(f2sat, f2sat_check,decimal=delta)
    assert_almost_equal(zsat,zsat_check,decimal=delta)
    
    
