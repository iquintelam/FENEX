import sys
import os
from FENEX import read_test_system 
import FENEX as FENEX




integ_type =1
Npoints,f1new,f,free_energy,z,cov,stats = FENEX.read_test_system() 
int_1 = FENEX.Integrate(f1new,f,free_energy,z,cov,'decoupled')
FENEX.Integrate.estimate_coexistence(int_1)
print(int_1.free_energy,int_1.f2new)
results_sat = FENEX.Integrate.refine_coexistence(int_1)
print(int_1.free_energy_sat)

# int_1 = FENEX.SimulationData(Npoints,f1new,f,free_energy,z,cov,stats)
# print(f"{int_1.f[0,:]}")

# integ_type='decoupled'
# #free_energy,f2new = FENEX.estimate_coexistence(integ_type ,int_1.f1new,int_1.f,int_1.free_energy, int_1.z, int_1.cov) 
# int_1.free_energy,f2new = FENEX.estimate_coexistence(integ_type ,int_1)


# print(int_1.free_energy,f2new)
# results_sat = FENEX.refine_coex(int_1)
# print(results_sat['free_energy_sat'])
# dir_path = os.path.dirname(os.path.realpath(__file__))
# data1 = os.path.join(dir_path, "..", "data", "simulation_data_one.dat")
# Npoints,f1new,f,free_energy, z, cov,stats = read_input_MC(data1)
# # print(Npoints)
# # print(f1new)
# # print(f)
# # print(free_energy)
# # print(z)
# # print(cov)
# # print(stats)
# # print(np.shape(f),np.shape(free_energy),np.shape(z),np.shape(cov))
# #print(f[:,0,:])

# f2new = calculate_first_point(f1new,f,free_energy, z, cov) 
# #print(f2new)


# dir_path = os.path.dirname(os.path.realpath(__file__))
# data1 = os.path.join(dir_path, "..", "data", "simulation_data.dat")
# #data1 = r"data\simulation_data.dat"
# Npoints,f1new,f,free_energy, z, cov,stats = read_input_MC(data1)
# print(Npoints)
# print(f1new)
# print(f)
# print(free_energy)
# print(z)
# print(cov)
# print(stats)

# free_energy,f2new = calculate_next_point(1 ,f1new,f,free_energy, z, cov) 
# print(free_energy,f2new)
# zsat,free_energy_zsat,enesat,f2sat = calc_zsat(Npoints,free_energy,z,cov,f,stats[1,:,:])
# print(zsat)
# zsat_check =([[-.2095241E+01 ,-.1595087E+01],
#        [ -.1967657E+01 ,-.1498863E+01],
#        [ 0.2225709E+01 ,0.2463211E+01],
#        [ 0.2225398E+01 ,0.2459797E+01]])

# zsat_check =[[[-.2095241E+01 ,-.1595087E+01], [ -.1967657E+01 ,-.1498863E+01]],
#  [[ 0.2225709E+01 ,0.2463211E+01], [ 0.2225398E+01 ,0.2459797E+01]]]

# f2sat_check = [8.1685 ,8.2736]
# print((zsat_check)-(zsat))

# #df = delta_f(f1new,f2new,f)
# f1sat = 1.10000000000000
# f2sat = 8.13763782732879
# z = np.array([[[-2.08372931917015, -1.57569717110201]],[[ 2.23293244807807 , 2.47249643756598]]])
# cov = np.array([[[0.96084424, 0.99139901]],[[0.04958482, 0.06373951]],[[0.07901728 ,0.13309644]]])
# f = np.array([[[1.15     ,  1.15      ]],[[8.02277426, 8.02277426]]])
# df = np.array([[-4.999999999999982E-002,-4.999999999999982E-002],[ 0.114863567636570 , 0.114863567636570]])
# free_energy = np.array([[0.0349, 0.    ]])
# coef=poly_coefficients(df,z[:,0,:],cov[:,0,:])

# coefcheck =np.array([[  0.104186465958507     ,  7.878485855510022E-002  ],
#                      [  0.256482587277708     ,  0.283999761787539  ],
#                      [  2.402110607439918E-003,  2.478497537156032E-003  ],
#                      [  6.542042738191447E-004,  8.409560465083100E-004  ],
#                      [ -4.538103315873886E-004, -7.643965721313706E-004  ],
#                      [ -0.00395086            , -0.00665482 ]])
# delta = 1e-6
# print(coef-coefcheck)

# f1new = 1.10000000000000
# f2new = first_guess_newton(free_energy[0,:],z[:,0,:],f1new,f[:,0,:])

# fun = f_first_point(f2new,free_energy[0,:],z[:,0,:],cov[:,0,:],f1new,f[:,0,:])
# funcheck = -1.940106807892805E-004
# print( df_first_point(f2new,free_energy[0,:],z[:,0,:],cov[:,0,:],f1new,f[:,0,:]))
# f2new = optimize.newton(f_first_point, f2new0,fprime=df_first_point,
#                           args=(free_energy[0,:],z[:,0,:],cov[:,0,:],f1new,f[:,0,:]),
#                           maxiter=500)

# np.savetxt(r'data/data_array.csv', (np.transpose([f[0,:,0],f[1,:,0],f2sat,free_energy_zsat[:,0],zsat[0,:,0],zsat[0,:,1],zsat[1,:,0],zsat[1,:,1],enesat[:,0],enesat[:,1]])), delimiter=' ', header='f1 f2 f2sat FEsat z1satI z1satII z2satI z2satII enesatI enesatII', comments="")
# with open(r"data/data_array.csv", "a") as file:
#         file.write(format(f1new, '.8f') + " ")
#         file.write(format(f2new, '.8f'))

# fig1, ax1 = plt.subplots()
# ax1=plt.gca()
# ax1.plot(f[0,:,0], f[1,:,0], 'o',color='#1690B8',label='Near coex')
# ax1.plot(f[0,:,0], f2sat, 'o', color="black",label='Coex')
# ax1.legend(frameon=False,ncol=1,fontsize=20)
# ax1.tick_params(labelsize=20)
# ax1.set_ylabel(r"$f_2$", fontsize=22)
# ax1.set_xlabel(r"$f_1$", fontsize=22)

# fig2, ax2 = plt.subplots()
# ax2=plt.gca()
# ax2.plot(f[0,:,0], free_energy_zsat[:,0], 'o',color='#1690B8',label='I')
# ax2.plot(f[0,:,0], free_energy_zsat[:,1], 'o', color="black",label='II')
# ax2.legend(frameon=False,ncol=1,fontsize=20)
# ax2.tick_params(labelsize=20)
# ax2.set_ylabel(r"$\phi$", fontsize=22)
# ax2.set_xlabel(r"$f_1$", fontsize=22)



# fig3, ax3 = plt.subplots()
# ax3=plt.gca()
# ax3.plot(f[0,:,0], zsat[0,:,0], 'o',color='#1690B8',label='I')
# ax3.plot(f[0,:,0], zsat[0,:,1], 'o', color="black",label='II')
# ax3.legend(frameon=False,ncol=1,fontsize=20)
# ax3.tick_params(labelsize=20)
# ax3.set_ylabel(r"$z_1$", fontsize=22)
# ax3.set_xlabel(r"$f_1$", fontsize=22)

# fig4, ax4 = plt.subplots()
# ax4=plt.gca()
# ax4.plot(f[0,:,0], zsat[1,:,0], 'o',color='#1690B8',label='I')
# ax4.plot(f[0,:,0], zsat[1,:,1], 'o', color="black",label='II')
# ax4.legend(frameon=False,ncol=1,fontsize=20)
# ax4.tick_params(labelsize=20)
# ax4.set_ylabel(r"$z_2$", fontsize=22)
# ax4.set_xlabel(r"$f_1$", fontsize=22)

# fig5, ax5 = plt.subplots()
# ax5=plt.gca()
# ax5.plot(f[0,:,0], enesat[:,0], 'o',color='#1690B8',label='I')
# ax5.plot(f[0,:,0], enesat[:,1], 'o', color="black",label='II')
# ax5.legend(frameon=False,ncol=1,fontsize=20)
# ax5.tick_params(labelsize=20)
# ax5.set_ylabel(r"$<u>$", fontsize=22)
# ax5.set_xlabel(r"$f_1$", fontsize=22)


# plt.show()