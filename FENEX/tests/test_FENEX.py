"""
Unit and regression test for the FENEX package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

#import FENEX
sys.path.append('../')
from FENEX.Data_Proc import *
from FENEX.functions import *


import matplotlib.pyplot as plt

def test_FENEX_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "FENEX" in sys.modules

#f,z,cov,stats = read_covdata_MC(r"data/iso.dat")
data1 = sys.argv[1]
Npoints,f1new,f,free_energy, z, cov,stats = read_input_MC(data1)

free_energy,f2new = calculate_next_point(1 ,f1new,f,free_energy, z, cov) 


zsat,free_energy_zsat,enesat,f2sat = calc_zsat(Npoints,free_energy,z,cov,f,stats[1,:,:])

np.savetxt(r'data/data_array.csv', (np.transpose([f[0,:,0],f[1,:,0],f2sat,free_energy_zsat[:,0],zsat[0,:,0],zsat[0,:,1],zsat[1,:,0],zsat[1,:,1],enesat[:,0],enesat[:,1]])), delimiter=' ', header='f1 f2 f2sat FEsat z1satI z1satII z2satI z2satII enesatI enesatII', comments="")
with open(r"data/data_array.csv", "a") as file:
    file.write(format(f1new, '.8f') + " ")
    file.write(format(f2new, '.8f'))
fig1, ax1 = plt.subplots()
ax1=plt.gca()
ax1.plot(f[0,:,0], f[1,:,0], 'o',color='#1690B8',label='Near coex')
ax1.plot(f[0,:,0], f2sat, 'o', color="black",label='Coex')
ax1.legend(frameon=False,ncol=1,fontsize=20)
ax1.tick_params(labelsize=20)
ax1.set_ylabel(r"$f_2$", fontsize=22)
ax1.set_xlabel(r"$f_1$", fontsize=22)

fig2, ax2 = plt.subplots()
ax2=plt.gca()
ax2.plot(f[0,:,0], free_energy_zsat[:,0], 'o',color='#1690B8',label='I')
ax2.plot(f[0,:,0], free_energy_zsat[:,1], 'o', color="black",label='II')
ax2.legend(frameon=False,ncol=1,fontsize=20)
ax2.tick_params(labelsize=20)
ax2.set_ylabel(r"$\phi$", fontsize=22)
ax2.set_xlabel(r"$f_1$", fontsize=22)



fig3, ax3 = plt.subplots()
ax3=plt.gca()
ax3.plot(f[0,:,0], zsat[0,:,0], 'o',color='#1690B8',label='I')
ax3.plot(f[0,:,0], zsat[0,:,1], 'o', color="black",label='II')
ax3.legend(frameon=False,ncol=1,fontsize=20)
ax3.tick_params(labelsize=20)
ax3.set_ylabel(r"$z_1$", fontsize=22)
ax3.set_xlabel(r"$f_1$", fontsize=22)

fig4, ax4 = plt.subplots()
ax4=plt.gca()
ax4.plot(f[0,:,0], zsat[1,:,0], 'o',color='#1690B8',label='I')
ax4.plot(f[0,:,0], zsat[1,:,1], 'o', color="black",label='II')
ax4.legend(frameon=False,ncol=1,fontsize=20)
ax4.tick_params(labelsize=20)
ax4.set_ylabel(r"$z_2$", fontsize=22)
ax4.set_xlabel(r"$f_1$", fontsize=22)

fig5, ax5 = plt.subplots()
ax5=plt.gca()
ax5.plot(f[0,:,0], enesat[:,0], 'o',color='#1690B8',label='I')
ax5.plot(f[0,:,0], enesat[:,1], 'o', color="black",label='II')
ax5.legend(frameon=False,ncol=1,fontsize=20)
ax5.tick_params(labelsize=20)
ax5.set_ylabel(r"$<u>$", fontsize=22)
ax5.set_xlabel(r"$f_1$", fontsize=22)


plt.show()