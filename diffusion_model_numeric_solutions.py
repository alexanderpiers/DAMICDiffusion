# Numerical Methods for Solving Diffusion Equation
from scipy.integrate import odeint, cumtrapz, trapz
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import brewer2mpl

def param():
    """
    Set the parameters necessary for solving the PDE numerically
    """

    mu = 0.200 # m^2/Vs
    q = 1.6e-19 # C
    T = 130 # K
    kb = 1.38e-23 # m^2 kg / K s^2
    epsSi = 11.68 * 8.85e-12 # C/Vm
    Dp = mu * kb * T / q
    qInit = 200*q

    # mu = 1
    # q = 1
    # T = 1
    # kb = 1
    # epsSi = 1
    # Dp = 1

    return mu, q, T, kb, epsSi, Dp, qInit


def GroomDiffusionEquation(length=1e-6):
    """
    Solves the PDE outlined in the Groom papers in spherical coordinates (for simple case of no E field)
    """
    rmax = 60.e-6 # 60 um
    # rmax = 10
    nr = 1000
    dr = rmax / nr
    t = np.linspace(0, 0.1e-6, 1000) # 10 us max time

    mu, q, T, kb, eps, Dp, qInit = param()

    # Define initial conditions
    p0 = np.zeros(nr)
    p0[0:int(length/dr)] = qInit/length**3
    # p0[1] = p0[0]

    chargeDensity = odeint(GroomDiffusionFunc, p0, t, args=(dr, Dp))

    # Normalize the charge density distribution
    chargeDensitySum = np.sum(chargeDensity, axis=1)*dr
    chargeDensity /= np.reshape(np.repeat(chargeDensitySum, nr), chargeDensity.shape)

    return chargeDensity, dr, t

def GroomDiffusionFunc(p, t, dr, Dp):
    """
    Generates the array of ODEs to be used by ODEint
    """
    N = len(p) - 1
    r = np.arange(0, dr*(N+1), dr)
    func = np.zeros(len(p))

    func[1:N] = (Dp / dr**2) * (p[2:N+1] - 2*p[1:N] + p[0:N-1]) + Dp / (r[1:N] * dr) * (p[2:N+1] - p[0:N-1])
    func[0] = 2 * Dp / dr**2 * (p[1] - p[0])
    func[N] = 0 #2 * Dp / dr**2 * (p[N-1] - p[N])

    return func

def CoulombDiffusionFunc(p, t, dr, Dp, eps, mu):
    """
    Generates the array of ODEs for use with the coulomb repulsion between carriers case
    """

    # Defining values and parameters
    N = len(p) - 1
    r = np.arange(0, dr * (N+1), dr)
    odes = np.zeros(len(p))

    # Defining ODE equations
    # ODE broken into the d^2p/dr^2, dp/dr, and p lines
    # odes[1:N] = (Dp / dr**2) * (p[2:N+1] - 2*p[1:N] + p[0:N-1]) + (Dp / (r[1:N] * dr) - (mu * p[1:N] * r[1:N]) / (3* eps * dr)) * (p[2:N+1] - p[0:N-1]) - mu * p[1:N]**2 / eps # for constant rho gauss' law

    odes[1:N] = \
            (Dp / dr**2) * (p[2:N+1] - 2*p[1:N] + p[0:N-1]) + \
            ((Dp / (r[1:N] * dr) - (mu * cumtrapz(r[0:N]**2 * p[0:N], r[0:N], dr)) / (eps * r[1:N]**2 * dr))) * (p[2:N+1] - p[0:N-1]) - \
             mu * p[1:N]**2 / eps
    odes[0] = 2 * Dp / dr**2 * (p[1] - p[0]) - mu * p[0]**2 / eps
    odes[N] = 0

    return odes

def integrateRho(rho, r, dr):
    """
    Performs the integral integral[rho*r'^2*dr''] from r'=0 to r'=r for all values of r passed
    """
    integral = np.zeros(r.size)
    for i in range(r.size):
        integral[i] = np.sum(rho[0:i+1] * r[0:i+1]**2 * dr)

    return integral

def CoulombDiffusionEquation(qInit=125*1.6e-19, length=2.e-6):
    """
    Solves the PDE outlined in the Groom papers in spherical coordinates for case of mutual repulsion between generated charge
    """

    # define constants
    rmax = 100.e-6 # 100 um
    nr = 1000
    dr = rmax / nr
    t = np.linspace(0, 0.1e-6, 2000) # 10 us max time


    mu, q, T, kb, eps, Dp, _ = param()

    # Define initial conditions
    p0 = np.zeros(nr)
    p0[0:int(length/dr)] = qInit/length**3

    p = odeint(CoulombDiffusionFunc, p0, t, args=(dr, Dp, eps, mu))

    # Normalize the charge density distribution
    pSum = np.sum(p, axis=1)*dr
    p /= np.reshape(np.repeat(pSum, nr), p.shape)


    return p, dr, t

def PiecewiseDiffusionEquation(teff=10e-9, qinit=1.6e-17, length=1e-6, rmax=100.e-6, nr=1000, tend=0.1e-6):
    """
    Solving the diffusion equation for two regimes 1) when electron and hole charge clouds overlap and therefore no efield present 2) charge clouds seperate and mutual coulomb repulsion matters. Regimes seperates by some t_eff
    """

    dr = rmax / nr
    t = np.linspace(0, tend, 1000)
    dt = t[1] - t[0]
    r = np.arange(0, rmax, dr)

    mu, q, T, kb, eps, Dp, _ = param()

    # define the two time regimes
    t1regime = t[0:int(teff/dt)]
    t2regime = t[int(teff/dt):]

    # Define initial conditions
    p0 = np.zeros(nr)
    p0[0:int(length/dr)] = 3 * qinit / (4 * np.pi * length**3)


    # solving equations in first regime
    if t1regime.size:
        p1regime = odeint(GroomDiffusionFunc, p0, t1regime, args=(dr, Dp))

        # initial conditions for next regime, including normalizing rho (making sure charge is conserved)
        p1 = p1regime[-1,:]
    else:
        p1 = p0
        p1regime = np.array([])

    # solve the second regime with
    p2regime =  odeint(CoulombDiffusionFunc, p1, t2regime   , args=(dr, Dp, eps, mu))

    if p1regime.size:
        p = np.concatenate((p1regime, p2regime), axis=0)
    else:
        p = p2regime

    # normalize p
    alpha = qinit / trapz(4 * np.pi * p * r**2, r, dr)
    print(alpha)
    p *= np.repeat(np.reshape(alpha, (alpha.size, 1)), nr, axis=1)

    return p, r, t

def testPiecewiseSolution():

    fig, ax = plt.subplots()

    teff = np.linspace(0, 15e-9, 10)
    pall = []
    tsample = 30e-9
    E = 5000
    w = 3.8
    ncharge = E/w
    legstring = []
    halfmaxdist = []

    for te in teff:
        p, rr, t = PiecewiseDiffusionEquation(teff=te, qinit=ncharge*1.6e-19)
        dt = np.diff(t)[0]
        ptime = p[int(tsample/dt),:]
        pall.append(ptime)
        halfmaxdist = rr[np.nonzero(ptime < ptime[0]/2)[0][0]]
        legstring.append('teff=%.2f ns. FWHM=%.2f $\mu m$' % (te*1e9, 2*halfmaxdist*1e6))

    pall = np.array(pall)
    ax.plot(rr*1e6, pall.T, linewidth=2)
    ax.legend(legstring, fontsize=16)
    ax.set_ylabel('Charge Density', fontsize=16)
    ax.set_xlabel('Radial Distance ($\mu m$)', fontsize=16)
    ax.set_title('Charge Distribution from Using Piecewise Diffusion model. E=1keV, Ri=1$\mu m$', fontsize=18)

    plt.show()
    return None

def ComparisonOfEandNoEField():

    fig, axs = plt.subplots(2,2)
    axs = axs.flatten()

    pNoE, dr, t = GroomDiffusionEquation()
    pE, dr, t = CoulombDiffusionEquation()

    rAxis = np.arange(0, p.shape[1]*dr*1e6, dr*1e6)
    dt = t[1] - t[0]

    tstartIndex = 1
    tstepIndex = 30

    for i, ax in enumerate(axs):
        ax.set_ylabel('Charge Density (au)', fontsize=14)
        if i > 1:
            ax.set_xlabel('Radial Distance ($\mu m$)', fontsize=14)
        ax.set_title('Time=%.2f $\mu s$' % ((tstartIndex + i*tstepIndex)*dt*1e6), fontsize=16)

        ax.plot(rAxis, pE[(tstartIndex + i*tstepIndex),:]/1e6, color='k', linewidth=3, label='Coulomb')
        ax.plot(rAxis, pNoE[(tstartIndex + i*tstepIndex),:]/1e6, color='r', linewidth=3, label='No Coulomb')

        ax.legend()

    fig.suptitle('Radial Probability Distribution of Charge Density of Diffusion for 1keV deposited Energy', fontsize=18)
    plt.show()

    return None

def EnergyDependenceOfCoulomb():

    fig, ax = plt.subplots()
    wSi = 4.26
    energy = [100, 200,  500, 1000, 2000, 5000,]# 10000]
    ncharge = (np.array(energy)/wSi).astype(int)*1.6e-19
    length=5.e-6
    timeslice = 0.05e-6 # 50 ns

    bmap = brewer2mpl.get_map('Set2', 'Qualitative', 8).mpl_colors


    for i, q in enumerate(ncharge):
        p, dr, t = CoulombDiffusionEquation(q, length=length)
        rAxis = np.arange(0, p.shape[1]*dr* 1e6, dr* 1e6)
        dt = t[1] - t[0]
        tIndex = int(timeslice/dt)

        ax.plot(rAxis, p[tIndex,:]/1e6, linewidth=3, color=bmap[i])

    ax.legend([str(x) for x in energy])
    ax.set_ylabel('Charge Density', fontsize=16)
    ax.set_xlabel('Radial Distance ($\mu m$)', fontsize=16)
    ax.set_title('Charge Distribution from Diffusion at 50 ns for varying Energy (legend in eV)', fontsize=18)

    plt.show()
    return None

def InitialRadiusDependence():

    fig, ax = plt.subplots()
    length = [0.25, 0.5, 1, 2, 3, 5]
    timeslice = 2e-9 # 50 ns


    bmap = brewer2mpl.get_map('Set2', 'Qualitative', 8).mpl_colors

    for i, l in enumerate(length):
        p, dr, t = GroomDiffusionEquation(length=l*1e-6)
        rAxis = np.arange(0, p.shape[1]*dr* 1e6, dr* 1e6)
        dt = t[1] - t[0]
        tIndex = int(timeslice/dt)

        ax.plot(rAxis, p[tIndex,:]/1e6, linewidth=3, color=bmap[i])

    ax.legend([str(x) + ' $\mu m$'  for x in length])
    ax.set_ylabel('Charge Density', fontsize=16)
    ax.set_xlabel('Radial Distance ($\mu m$)', fontsize=16)
    ax.set_title('Charge Distribution from Diffusion at 50 ns for varying Initial Radius Charge', fontsize=18)

    plt.show()
    return None

def energyDependencePiecewise():

    fig, ax = plt.subplots()
    teff = 2e-9 # nanoseconds
    rinit = 0.5e-6
    teffArray = np.linspace(1e-9, 3e-9, 5)
    tsample = 20e-9
    E = np.linspace(100, 15000, 10)
    w = 3.8
    ncharge = E/w
    pall = []
    legstring = []
    fwhmMatrix =  []

    # for teff in teffArray:
    fwhmArray = []
    for n in ncharge:
        p, rr, t = PiecewiseDiffusionEquation(teff=teff, qinit=n*1.6e-19, length=rinit)
        dt = np.diff(t)[0]
        print(dt)
        pfinal = p[int(tsample/dt),:]
        fwhm = rr[np.nonzero(pfinal < pfinal[0]/2)[0][0]]
        fwhmArray.append(fwhm)
        legstring.append('E=%.2f keV, FWHM=%.1f'%(n*w/1000, fwhm*2*1e6))
        pall.append(pfinal/n)
    fwhmMatrix.append(fwhmArray)

    # fwhmMatrix = np.array(fwhmMatrix)
    # ax.plot(E/1000, fwhmMatrix.T*1e6, linewidth=3)
    # ax.legend(['teff=%0.1f ns'%(x*1e9) for x in teffArray], fontsize=16)
    # ax.set_ylabel('FWHM ($\mu m$)', fontsize=16)
    # ax.set_xlabel('Energy (keV)', fontsize=16)
    # ax.set_title('Energy vs FWHM for different $t_{eff}$', fontsize=18)


    pall = np.array(pall)
    ax.plot(rr*1e6, pall.T, linewidth=2)
    ax.legend(legstring, fontsize=16)
    ax.set_ylabel('Charge Density Distribution', fontsize=16)
    ax.set_xlabel('Radial Distance ($\mu m$)', fontsize=16)
    ax.set_title('Charge Distribution from Using Piecewise Diffusion model After 25 ns. teff=%.1f ns, Ri=1$\mu m$'%(teff*1e9), fontsize=18)
    plt.show()


def projectSphericalonXY(f, r, dphi, dtheta, dx, dy):
    """
    Takes the radial distribution function and projects it onto the xy plane.
    Returns a slice on the x-axis (but for spherically symmetric distributions it shouldn't matter)

    Convention of physcicist. Phi is the azimuthal angle and theta is the polar angle
    """

    # define the infintesimal changes
    phi = np.arange(0, 2*np.pi, dphi)
    # theta = np.arange(0, np.pi, dtheta)
    dr = np.diff(r)[0]

    xval = []
    yval = []
    weight = []

    for i in r:
        for j in phi:
            for k in theta:
                xval.append(i * np.sin(j) * np.cos(k))
                yval.append(i * np.sin(j) * np.sin(k))
                weight.append(i**2 * dr * dphi)


    f, edges = np.histogram(xval, bins=int(max(xval)/dx), weights=weight)
    x = edges[1:] - np.diff(edges)/2

    return f, x

def testProjectSphere():

    r = np.linspace(0, 10, 2000)
    dphi = 2*np.pi/2000
    f = np.ones(r.size)
    dtheta = 1
    dx = 0.1

    f, x = projectSphericalonX(f, r, dphi, dtheta, dx)

    fig, ax = plt.subplots()
    ax.plot(x, f, linewidth=3)
    plt.show()


if __name__ == '__main__':





    # ComparisonOfEandNoEField()
    # EnergyDependenceOfCoulomb()
    # InitialRadiusDependence()
    # PiecewiseDiffusionEquation(qinit=1)
    # testPiecewiseSolution()
    # testProjectSphere()
    energyDependencePiecewise()
    plt.show()




