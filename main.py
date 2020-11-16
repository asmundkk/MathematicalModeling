from sympy import symbols, Matrix
from functions import newton_temperature, plot_sensitivity

# defining symbols used in sympy
Te, Ta = symbols('Te, Ta')                                                                                  # task 1
lw_absorbtivity_gases, sw_absorbtivity_gases = symbols('lw_absorbtivity_gases, sw_absorbtivity_gases')      # task 2

# constants from figure 1
p_sun = 341                                # W m**-2
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
sw_absorbtivity_gases_num = 0.1451         # - , numerical value

lw_cloud_scattering_coeff = 0.195          # -
lw_earth_reflectivity = 0                  # -
lw_cloud_absorbtivity = 0.62               # - , numerical value
lw_absorbtivity_gases_num = 0.8258         # -

earth_emmisivity = 1                       # - , BB assumption
atm_emmisivity = 0.875                     # - , no BB assumption for the atmosphere
asymmerty_factor = 0.618                   # - , does this number come from the figure

sensible_heat_flux = 3                     # W m**-2 K**-1
latent_heat_flux = 4                       # W m**-2 K**-1

sigma = 5.67*10**-8                        # W m**-2 K**-4, Boltzmann's constant

# basic equations
p_e = earth_emmisivity * sigma * Te**4
p_a = atm_emmisivity * sigma * Ta**4

p_sensible = sensible_heat_flux*(Te -Ta)   # Newtons law of cooling, making these positive
p_latent = latent_heat_flux*(Te - Ta)      # Newtons law of cooling

# reflection
sw_reflection = cloud_cower * sw_cloud_scattering_coeff + sw_molecular_scattering_coeff

# atmosphere
# power in and absorption
sw_p_in = p_sun * (1 - sw_reflection)
sw_p_absorbed = sw_p_in * (cloud_cower * sw_absorptitity_clouds + sw_absorbtivity_ozone + sw_absorbtivity_gases)

# defining the equations for atmosphere
lw_p_in = p_e + p_latent + p_sensible
lw_p_absorped = lw_p_in * (cloud_cower * lw_cloud_absorbtivity + lw_absorbtivity_gases)
p_absorbed = sw_p_absorbed + lw_p_absorped
p_in = lw_p_in + sw_p_in
p_out_a = heat_flux_to_surface + p_a

fa = p_in - p_out_a - p_absorbed

# earth surface
p_out_e = p_e + p_latent + p_sensible + p_evap + p_thermal

fe = -p_out_e

F = Matrix([fa, fe])
T = Matrix([Ta, Te])

#  initial conditions
Ta0 = 273 - 20              # -20 degC
Te0 = 273                   #  0  degC
T0 = Matrix([Ta0, Te0])

# relative error
e = abs((T0[0] - T0[1]) / (T0[0] * 100))

A = Matrix([lw_absorbtivity_gases, sw_absorbtivity_gases])

if __name__ == '__main__':

    #  question 1
    T_num = newton_temperature(F, T0, T, lw_absorbtivity_gases, sw_absorbtivity_gases, lw_absorbtivity_gases_num, sw_absorbtivity_gases_num, err=e)

    #  question 2
    plot_sensitivity(F, A, T, T0, lw_absorbtivity_gases, sw_absorbtivity_gases, lw_absorbtivity_gases_num, sw_absorbtivity_gases_num)
