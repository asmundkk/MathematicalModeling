from sympy import Matrix, diff, sqrt
from numpy import linspace
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

def jacobian(F, x):
    """
    The jacobian of a 2x2 matrix

    F: list of functions, Matrix(2, 1)
    x: list of variables, Matrix(2, 1)
    """

    J = Matrix([[diff(F[0], x[0]), diff(F[0], x[1])],
               [diff(F[1], x[0]), diff(F[1], x[1])]])

    print("jacobian:\n",J)
    return J

def newton_temperature(F, T0, T, lw_absorbtivity_gases, sw_absorbtivity_gases, lw_absorbtivity_gases_num, sw_absorbtivity_gases_num, err=0.001):
    """
    Finds the roots of a system of equations using the newton method

    F: list of function, Matrix(2, 1)
    T0: list of initial temperatures, Matrix(2, 1)
    T: list of temperatures, Matrix(2, 1)
    lw_absorbtivity_gases: coefficient of absorbation for long waves
    sw_absorbtivity_gases: coefficient of absorbation for short waves
    lw_absorbtivity_gases_num: coefficient of absorbation for long waves numerical value
    sw_absorbtivity_gases_num: coefficient of absorbation for short waves numerical value
    err: error to check if the convergence criterion is reached
    """

    Ta, Te  = T[0] ,T[1]

    J = jacobian(F, T)
    J_inv = J.inv()
    T = T0

    i = 0
    error = 100             # some high number

    while error > err and i <= 100:

        Jn = J_inv.subs([(Ta, T[0]), (Te, T[1]), (lw_absorbtivity_gases, lw_absorbtivity_gases_num), (sw_absorbtivity_gases, sw_absorbtivity_gases_num)])
        Fn = F.subs([(Ta, T[0]), (Te, T[1]), (lw_absorbtivity_gases, lw_absorbtivity_gases_num), (sw_absorbtivity_gases, sw_absorbtivity_gases_num)])
        Tp1 = T - Jn * Fn

        print("i = ", i)
        print("Ta = ", T[0])
        print("Te = ", T[1])
        print("relative error = ", str(error / err), "\n")

        error = sqrt(((T[0]-Tp1[0])/T[0])**2 + ((T[1]-Tp1[1])/T[1])**2)
        T = Tp1
        i += 1

    return T

def plot_sensitivity(F, A, T, T0, lw_absorbtivity_gases, sw_absorbtivity_gases, lw_absorbtivity_gases_num, sw_absorbtivity_gases_num, e=0.1, n=20):
    """
    Plot the sensitivity of temperature model made in task 1

    F: list of functions, Matrix(2, 1)
    A: list of absorptivities, Matrix(2, 1)
    T: list of temperatures, Matrix(2, 1)
    T0: list of initial temperatures, Matrix(2, 1)
    lw_absorbtivity_gases: absorotivity of long wave radiation
    sw_absorbtivity_gases: absorotivity of long short radiation
    lw_absorbtivity_gases_num: absorotivity of long wave radiation numrical value
    sw_absorbtivity_gases_num: absorotivity of short wave radiation numrical value
    e: how much to vary the absoptivity coefficient
    n: number of nodes of the x-axis
    """

    a_lw = linspace(lw_absorbtivity_gases_num - e, lw_absorbtivity_gases_num + e, num=n)
    a_sw = linspace(sw_absorbtivity_gases_num - e, sw_absorbtivity_gases_num + e, num=n)

    Ta_lw = [0] * n
    Te_lw = [0] * n
    Ta_sw = [0] * n
    Te_sw = [0] * n

    for i in range(n):
        Alw = Matrix([a_lw[i], sw_absorbtivity_gases_num])
        Asw = Matrix([lw_absorbtivity_gases_num, a_sw[i]])
        #  Ta_lw[i], Te_lw[i] = newton_radiation(F ,T, Alw, T0)
        #  Ta_sw[i], Te_sw[i] = newton_radiation(F ,T, Asw, T0)
        Ta_lw[i], Te_lw[i] = newton_temperature(F ,T0, T, lw_absorbtivity_gases, sw_absorbtivity_gases, a_lw[i], sw_absorbtivity_gases_num)
        Ta_sw[i], Te_sw[i] = newton_temperature(F ,T0, T, lw_absorbtivity_gases, sw_absorbtivity_gases, lw_absorbtivity_gases_num, a_sw[i])

    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.set_title(r'Sensitivity long and short wave radiation $H_2O$, $CO_2$ and $CH_4$')
    ax1.plot(a_lw, Ta_lw, label=r'$T_a$')
    ax1.plot(a_lw, Te_lw, label=r'$T_e$')
    ax1.set_ylabel(r'Temperature [K]')
    ax1.set_xlabel(r'Absorbtivity long wave radiation')
    ax1.legend()

    ax2.plot(a_sw, Ta_sw, label=r'$T_a$')
    ax2.plot(a_sw, Te_sw, label=r'$T_e$')
    ax2.set_ylabel(r'Temperature [K]')
    ax2.set_xlabel(r'Absorbtivity short wave radiation')
    ax2.legend()
    fig.savefig('/home/asmund/Uni/MattematiskModellering/Project/MathematicalModeling/sensitivity.png')
    plt.show()
