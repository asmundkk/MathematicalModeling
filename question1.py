import numpy as np
import sympy as sp
from numpy import array
from sympy import symbols, solve, Matrix, nonlinsolve, diff, sqrt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# defining symbols used in sympy
Te, Ta = symbols('Te, Ta')

# constants from figure 1
p_sun = 341                   # W m**-2
heat_flux_to_surface = 161                 # W m**-2, radiation coming through atm, and back radiation from atm should not be used
back_radiation_a = 333                     # W m**-2, radiation coming through atm, and back radiation from atm should not be used
surface_radiation = 356                    # W m**-2, lw radiation from earth                                   should not be used
p_thermal = 17                             # W m**-2, thermal heat
p_evap = 80                                # W m**-2, evaporation heat

# coefficients from Table 1
cloud_cower = 0.66                         # -
sw_molecular_scattering_coeff = 0.1065     # -
sw_cloud_scattering_coeff = 0.22           # -
sw_earth_reflectivity = 0.17               # -
sw_absorbtivity_ozone = 0.08               # -
sw_absorptitity_clouds = 0.1239            # -
sw_absorbtivity_gases = 0.1451             # -

lw_cloud_scattering_coeff = 0.195          # -
lw_earth_reflectivity = 0                  # -
lw_cloud_absorbtivity = 0.62               # -
lw_absorbtivity_gases = 0.8258             # -

earth_emmisivity = 1                       # - , BB assumption
atm_emmisivity = 0.875                     # - , no BB assumption for the atmosphere
asymmerty_factor = 0.618                   # - , does this number come from the figure

sensible_heat_flux = 3                     # W m**-2 K**-1
latent_heat_flux = 4                       # W m**-2 K**-1

sigma = 5.67*10**-8                        # W m**-2 K**-4, Boltzmanns constant

# basic equations
p_e = earth_emmisivity * sigma * Te**4
p_a = atm_emmisivity * sigma * Ta**4

print("p_e: ", p_e)
print("p_a: ", p_a)

p_sensible = sensible_heat_flux*(Te -Ta)   # Newtons law of cooling, making these positive
p_latent = latent_heat_flux*(Te - Ta)      # Newtons law of cooling

print("p_sensible:", p_sensible)
print("p_latent:", p_latent)

# reflection
sw_reflection = cloud_cower * sw_cloud_scattering_coeff + sw_molecular_scattering_coeff
print("sw_reflection:", sw_reflection)

# atmosphere
# power in and absorption
sw_p_in = p_sun * (1 - sw_reflection)
sw_p_absorbed = sw_p_in * (cloud_cower * sw_absorptitity_clouds + sw_absorbtivity_ozone + sw_absorbtivity_gases)
print("sw_p_in:", sw_p_in)
print("sw_p_absorbed:", sw_p_absorbed)

# defining the equations for atmosphere
lw_p_in = p_e + p_latent + p_sensible # + p_reflected_e, no reflection because of BB assumption
lw_p_absorped = lw_p_in * (cloud_cower * lw_cloud_absorbtivity + lw_absorbtivity_gases)
print("lw_p_in:", lw_p_in)
print("lw_p_absorped:", lw_p_absorped)

p_absorbed = sw_p_absorbed + lw_p_absorped
p_in = lw_p_in + sw_p_in
print("p_absorbed:", p_absorbed)

p_out_a = heat_flux_to_surface + p_a

fa = p_in - p_out_a - p_absorbed
# power out
#  p_in = sw_p_in + lw_p_in + sources
#  p_out = p_in - p_out_up - p_absorbed                      # p_stored = p_in - p_out + sources
#  p_out_up = p_out * asymmerty_factor
#  p_out_down = p_in - p_out_up                              # how much of the radiated heat that is going down
#  fa = p_out_down - heat_flux_to_surface

# earth surface
p_out_e = p_e + p_latent + p_sensible + p_evap + p_thermal
fe = -p_out_e

#  fe = p_e + p_latent + p_sensible + sources - surface_radiation

print("fa: ", fa)
print("fe: ", fe)
#  fe = atm_emmisivity * Ta**4 - 100

#  print("diff" ,diff(fa, Te))

J = Matrix([[diff(fa, Ta), diff(fa, Te)],
            [diff(fe, Ta), diff(fe, Te)]])

F = Matrix([fa, fe])
Ta0 = 273 - 20              # -20 degC
Te0 = 273                   #  0  degC
T0 = Matrix([Ta0, Te0])
e = abs((T0[0] - T0[1]) / (T0[0] * 100))          # relative error
print("e = ", e, "\n")

def newton(J, F, T0, err = 0.001):

    J_inv = J.inv()
    i = 0
    T = T0
    error = 100             # some high number

    while error > err and i <= 100:

        Jn = J_inv.subs([(Ta, T[0]), (Te, T[1])])
        Fn = F.subs([(Ta, T[0]), (Te, T[1])])
        print(Jn)
        print(Fn)
        print(T.shape)
        Tp1 = T - Jn * Fn

        print("i = ", i)
        print("Ta = ", T[0])
        print("Te = ", T[1])
        print("relative error = ", str(error / err), "\n")

        error = sqrt(((T[0]-Tp1[0])/T[0])**2 + ((T[1]-Tp1[1])/T[1])**2)
        T = Tp1
        i += 1

    return T

#  T_vec = np.linspace(200, 400, num=100)
#  T_vec = np.linspace(0, 500, num=100), # uncomment if plotting the fa, and fe

def plot_functions(F, T):

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    Te_mat = np.zeros([len(T), len(T)])
    Ta_mat = np.zeros_like(Te_mat)
    Z1 = np.zeros_like(Te_mat)
    Z2 = np.zeros_like(Te_mat)

    # homemade meshgrid
    for i in range(len(T)):
        Te_mat[:, i] = T[i]
        Ta_mat[i, :] = T[i]

    #  Xm, Ym = np.meshgrid(T, T)

    for i in range(len(T)):
        for j in range(len(T)):

            Z1[i, j] = fa.subs([(Te,Te_mat[i, j]), (Ta, Ta_mat[i, j])])    # this is very slow
            Z2[i, j] = fe.subs([(Te,Te_mat[i, j]), (Ta, Ta_mat[i, j])])

    ax = Axes3D(plt.gcf())
    ax.plot_surface(Te_mat, Ta_mat, Z1) # , label='fa')
    ax.plot_surface(Te_mat, Ta_mat, Z2) #, label='fe')
    ax.set_xlabel('Temperature [K]')
    ax.set_ylabel('Temperature [K]')
    ax.set_zlabel('f')
    #  ax.legend()
    plt.show()


if __name__ == '__main__':

    T = newton(J, F, T0, err=e) #TODO: test the newton implementation

    #  plot_functions(F, T_vec), so far yield very little result

    #  solve([fa, fe], [Ta, Te])
    #  solve([fa, fe], [Ta, Te], implisit=True)
    #  nonlinsolve([fa, fe], [Ta, Te])
