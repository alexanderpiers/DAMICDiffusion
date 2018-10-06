# fitting the hole drift velocity (Jacoboni 1977) to an exponential function

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def f(x, a, b):
    """ Defines the function to fit to. (1 - exp(x/a)) """
    return b * (1 - np.exp(-x/a))

def fJacob(x, vm, Ec, beta):
    """ Function to fit from Jacoboni paper """
    return vm * (x/Ec) / (1 + (x/Ec)**beta) ** (1./beta)

# Defining data (from Jacoboni, A Review of Some Charge Tranport Properties of SIlicon, 1977)
x = np.array([400, 700, 900, 2000, 3000, 4000, 5000, 6000, 7000, 10000, 20000])
y = np.array([0.9, 1.3, 1.6, 2.8, 3.5, 4.1, 4.6, 5.05, 5.5, 6, 7.2]) * 1e6

# find fit parameters
initialparam = [5e3, 1e7] # guess using paper data
popt,pcov = curve_fit(f, x, y, p0=initialparam)
print(popt)
xaxis = np.logspace(1, 5, 200)

# coefficients listed in the Jacoboni paper
T = 130
vmEst = 1.62e8 * T**(-0.52)
EcEst = 1.24 * T**(1.68)
betaEst = 0.46 * T**(0.17)

poptJac, pcovJac = curve_fit(fJacob, x, y, p0=[vmEst, EcEst, betaEst])

fig, ax = plt.subplots()
ax.loglog(x, y, 'X', markersize=10)
ax.loglog(xaxis, f(xaxis, popt[0], popt[1]), linewidth=3, color='k')
ax.loglog(xaxis, fJacob(xaxis, poptJac[0], poptJac[1], poptJac[2]), linewidth=3, color='r')
ax.loglog(xaxis, fJacob(xaxis, vmEst, EcEst, betaEst), linewidth=3, color='b')

# ax.set_xlim(1, 1e5)
# ax.set_ylim( 4e4, 2e7)
ax.legend(['data', r'exponential: $\beta(1-e^{-\frac{E}{\alpha}})$', r'Jacoboni: $v_m\frac{E/E_c}{(1+(E/E_c)^{\beta})^{1/\beta}}$', 'Jacoboni Paper Fit Parameters'], fontsize=16)
ax.set_xlabel('E Field (V/cm)', fontsize=16)
ax.set_ylabel('Hole Drift Velocity (cm/s)', fontsize=16)
ax.set_title('Hole Drift Velocity for High E Fields with Different Function Fits. Data from Jacoboni Paper.', fontsize=18)

# calculating low efield slopes
x1, x2 = 1, 2
slopeJacFit = (fJacob(x2, poptJac[0], poptJac[1], poptJac[2]) - fJacob(x1, poptJac[0], poptJac[1], poptJac[2]))/(x2-x1)
slopeJacEst = (fJacob(x2, vmEst, EcEst, betaEst) - fJacob(x1, vmEst, EcEst, betaEst))/(x2-x1)
expSlope = ( f(x2, popt[0], popt[1]) -  f(x1, popt[0], popt[1]))/(x2-x1)
print('Jacoboni Fit slope = %.2f' % slopeJacFit)
print('Jacoboni Estimate slope = %.2f' % slopeJacEst)
print('Exponential Fit slope = %.2f' % expSlope)
plt.show()
