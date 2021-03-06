"""
Unit and regression test for the FENEX package.
"""

# Import package, test suite, and other packages as needed
import sys
import os
from numpy.testing import assert_almost_equal

import pytest


from FENEX import *



def test_FENEX_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "FENEX" in sys.modules
def test_delta_f():
    f1sat = 1.10000000000000
    f2sat = 8.13763782732879
    f = np.array([[[1.15     ,  1.15      ]],[[8.02277426, 8.02277426]]])
    df = delta_f(f1sat,f2sat,f)
    dfcheck = np.array([[-4.999999999999982E-002,-4.999999999999982E-002],[ 0.114863567636570 , 0.114863567636570]])
    delta = 0.0001
    assert_almost_equal(df, dfcheck,decimal=delta)
def test_poly_coefficients():
    z = np.array([[[-2.08372931917015, -1.57569717110201]],[[ 2.23293244807807 , 2.47249643756598]]])
    cov = np.array([[[0.96084424, 0.99139901]],[[0.04958482, 0.06373951]],[[0.07901728 ,0.13309644]]])
    f = np.array([[[1.15     ,  1.15      ]],[[8.02277426, 8.02277426]]])
    df = np.array([[-4.999999999999982E-002,-4.999999999999982E-002],[ 0.114863567636570 , 0.114863567636570]])
    coef=poly_coefficients(df,z[:,0,:],cov[:,0,:])
    coefcheck =np.array([[  0.104186465958507     ,  7.878485855510022E-002  ],
                     [  0.256482587277708     ,  0.283999761787539  ],
                     [  2.402110607439918E-003,  2.478497537156032E-003  ],
                     [  6.542042738191447E-004,  8.409560465083100E-004  ],
                     [ -4.538103315873886E-004, -7.643965721313706E-004  ],
                     [ -0.00395086            , -0.00665482 ]])
    delta = 1e-6
    assert_almost_equal(coef, coefcheck,decimal=delta)
def test_first_guess_newton():
    f1new = 1.1
    free_energy = np.array([[0.0349, 0.    ]])
    f = np.array([[[1.15     ,  1.15      ]],[[8.02277426, 8.02277426]]])
    z = np.array([[[-2.08372931917015, -1.57569717110201]],[[ 2.23293244807807 , 2.47249643756598]]])
    f2new= first_guess_newton(free_energy[0,:],z[:,0,:],f1new,f[:,0,:])
    delta = 1e-6
    f2check = 8.27448824864485
    assert_almost_equal(f2new, f2check,decimal=delta)
def test_f_first_point():
    f1new = 1.1
    f2new = 8.274488248952627
    free_energy = np.array([[0.0349, 0.    ]])
    f = np.array([[[1.15     ,  1.15      ]],[[8.02277426, 8.02277426]]])
    z = np.array([[[-2.08372931917015, -1.57569717110201]],[[ 2.23293244807807 , 2.47249643756598]]])
    cov = np.array([[[0.96084424, 0.99139901]],[[0.04958482, 0.06373951]],[[0.07901728 ,0.13309644]]])
    fun = f_first_point(f2new,free_energy[0,:],z[:,0,:],cov[:,0,:],f1new,f[:,0,:])
    funcheck = -1.940106807892805E-004
    delta = 1e-9
    assert_almost_equal(fun, funcheck,decimal=delta)
def test_df_first_point():
    f1new = 1.1
    f2new = 8.274488248952627
    free_energy = np.array([[0.0349, 0.    ]])
    f = np.array([[[1.15     ,  1.15      ]],[[8.02277426, 8.02277426]]])
    z = np.array([[[-2.08372931917015, -1.57569717110201]],[[ 2.23293244807807 , 2.47249643756598]]])
    cov = np.array([[[0.96084424, 0.99139901]],[[0.04958482, 0.06373951]],[[0.07901728 ,0.13309644]]])
    fun = df_first_point(f2new,free_energy[0,:],z[:,0,:],cov[:,0,:],f1new,f[:,0,:])
    funcheck = -0.238705015716258
    delta = 1e-6
    assert_almost_equal(fun, funcheck,decimal=delta)
def test_calculate_first_point():
    f1new = 1.1
    free_energy = np.array([[0.0349, 0.    ]])
    f = np.array([[[1.15     ,  1.15      ]],[[8.02277426, 8.02277426]]])
    z = np.array([[[-2.08372931917015, -1.57569717110201]],[[ 2.23293244807807 , 2.47249643756598]]])
    cov = np.array([[[0.96084424, 0.99139901]],[[0.04958482, 0.06373951]],[[0.07901728 ,0.13309644]]])
    stats = np.array([[[ 0.37352467 , 0.224034  ]],[[-2.08372932 ,-1.57569717]]])
    f2new = calculate_first_point(f1new,f,free_energy, z, cov) 
    f2check = 8.27367550491069
    delta = 0.0001
    # assert function() to check if values are almost equal
    assert_almost_equal(f2new, f2check,decimal=delta)

def test_estimate_coexistence():
    Npoints = 2
    f1new = 1.05
    free_energy = np.array([[0.0349, 0.    ],[0.  ,   0.    ]])
    f = np.array([[[1.15,1.15],[1.1,1.1]],[[8.02277426,8.02277426],[8.13763783,8.13763783]]])
    z = np.array([[[-2.08372932,-1.57569717],[-1.95799492,-1.48165106]],[[2.23293245,2.47249644],[2.23128398,2.46899051]]])
    cov = np.array([[[0.96084424,0.99139901],[0.89934423,0.93119316]],[[0.04958482,0.06373951],[0.04330027,0.06763215]],[[0.07901728,0.13309644],[0.07107514,0.12661185]]])
    stats = np.array([[[0.37352467,0.224034],[0.366731,0.22013267]],[[-2.08372932,-1.57569717],[-1.95799492,-1.48165106]]])
    int_1 = Integrate(f1new,f,free_energy,z,cov,'decoupled',stats)
    Integrate.estimate_coexistence(int_1) 
    free_energy_check=[[0.0349  ,   0.        ],[0.3923189 , 0.36023006]]
    f2check=8.37242894118646
    delta =6e-4
    # assert function() to check if values are almost equal
    assert_almost_equal(int_1.f2new, f2check,decimal=delta)
    assert_almost_equal(int_1.free_energy, free_energy_check,decimal=delta)
    
    
def test_refine_coex():
    Npoints = 2
    f1new = 1.05
    free_energy=np.array([[0.0349  ,   0.        ],[0.3923189 , 0.36023006]])
    f = np.array([[[1.15,1.15],[1.1,1.1]],[[8.02277426,8.02277426],[8.13763783,8.13763783]]])
    z = np.array([[[-2.08372932,-1.57569717],[-1.95799492,-1.48165106]],[[2.23293245,2.47249644],[2.23128398,2.46899051]]])
    cov = np.array([[[0.96084424,0.99139901],[0.89934423,0.93119316]],[[0.04958482,0.06373951],[0.04330027,0.06763215]],[[0.07901728,0.13309644],[0.07107514,0.12661185]]])
    stats = np.array([[[0.37352467,0.224034],[0.366731,0.22013267]],[[-2.08372932,-1.57569717],[-1.95799492,-1.48165106]]])
    int_1 = Integrate(f1new,f,free_energy,z,cov,'decoupled',stats)
    Integrate.refine_coexistence(int_1)
    delta =6e-4
    f2sat_check = [8.1685 ,8.27357712677913]
    zsat_check =[[[-.2095241E+01 ,-.1595087E+01], [ -.1967657E+01 ,-.1498863E+01]],
    [[ 0.2225709E+01 ,0.2463211E+01], [ 0.2225398E+01 ,0.2459797E+01]]]
    free_energysat_check = [[0.36001651, 0.35971349],[0.69443783, 0.69398819]]
    assert_almost_equal(int_1.f2sat, f2sat_check,decimal=delta)
    assert_almost_equal(int_1.zsat,zsat_check,decimal=delta)
    assert_almost_equal(int_1.free_energy_sat,free_energysat_check,decimal=delta)