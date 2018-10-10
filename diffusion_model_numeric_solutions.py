# Numerical Methods for Solving Diffusion Equation
from scipy.integrate import odeint, cumtrapz, trapz
import numpy as np
import matplotlib.pyplot as plt

import brewer2mpl


# Solver functions
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

    return mu, q, T, kb, epsSi, Dp


def GroomDiffusionEquation(qinit=1.6e-17, length=1e-6, rmax=100.e-6, nr=1000, tend=0.1e-6, nt=1000):
    """
    Solves the PDE outlined in the Groom papers in spherical coordinates (for simple case of no E field)
    """

    # define constants
    dr = rmax / nr
    r = np.arange(0, rmax, dr)
    t = np.linspace(0, tend, nt) # 10 us max time
    mu, q, T, kb, eps, Dp = param()

    # Define initial conditions
    p0 = np.zeros(nr)
    p0[0:int(length/dr)] = qinit/length**3
    # p0[1] = p0[0]

    p = odeint(GroomDiffusionFunc, p0, t, args=(dr, Dp))

    # Normalize the charge density distribution to qinit
    alpha = qinit / trapz(4 * np.pi * p * r**2, r, dr)
    p *= np.repeat(np.reshape(alpha, (alpha.size, 1)), nr, axis=1)


    return p, r, t

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

def CoulombDiffusionEquation(qinit=1.6e-17, length=1e-6, rmax=100.e-6, nr=1000, tend=0.1e-6, nt=1000):
    """
    Solves the PDE outlined in the Groom papers in spherical coordinates for case of mutual repulsion between generated charge
    """

    # define constants
    dr = rmax / nr
    r = np.arange(0, rmax, dr)
    t = np.linspace(0, tend, nt) # 10 us max time
    mu, q, T, kb, eps, Dp = param()

    # Define initial conditions
    p0 = np.zeros(nr)
    p0[0:int(length/dr)] = qinit/length**3

    p = odeint(CoulombDiffusionFunc, p0, t, args=(dr, Dp, eps, mu))

    # Normalize the charge density distribution
    alpha = qinit / trapz(4 * np.pi * p * r**2, r, dr)
    p *= np.repeat(np.reshape(alpha, (alpha.size, 1)), nr, axis=1)


    return p, r, t

def PiecewiseDiffusionEquation(teff=10e-9, qinit=1.6e-17, length=1e-6, rmax=100.e-6, nr=1000, tend=0.1e-6, nt=1000):
    """
    Solving the diffusion equation for two regimes 1) when electron and hole charge clouds overlap and therefore no efield present 2) charge clouds seperate and mutual coulomb repulsion matters. Regimes seperates by some t_eff

    Inputs:
        teff - time at which the model switches from diffusion only to diffusion and coulomb interactions
        qinit - initial amount of charge distibuted in a sphere (C)
        length - radius of the initial distribution of charge (m)
        rmax - maximumt size of the simulation (m)
        nr - number of distance points. rmax/nr is the position step
        tend - final time of the simulation (s)
        nt - number of time points. tend/nt is the time step

    Outputs:
        p - radial charge density distribution. NxM matrix (timexdistance)
        r - radial axis Mx1 array
        t - time axis Nx1 array
    """

    dr = rmax / nr
    t = np.linspace(0, tend, nt)
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
    p *= np.repeat(np.reshape(alpha, (alpha.size, 1)), nr, axis=1)

    return p, r, t


# Analysis functions that can be called in main

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

    pNoE, rr, t = GroomDiffusionEquation()
    pE, rr, t = CoulombDiffusionEquation()

    dt = t[1] - t[0]

    tstartIndex = 1
    tstepIndex = 30

    for i, ax in enumerate(axs):
        ax.set_ylabel('Charge Density (au)', fontsize=14)
        if i > 1:
            ax.set_xlabel('Radial Distance ($\mu m$)', fontsize=14)
        ax.set_title('Time=%.2f $\mu s$' % ((tstartIndex + i*tstepIndex)*dt*1e6), fontsize=16)

        ax.plot(rr*1e6, pE[(tstartIndex + i*tstepIndex),:]/1e6, color='k', linewidth=3, label='Coulomb')
        ax.plot(rr*1e6, pNoE[(tstartIndex + i*tstepIndex),:]/1e6, color='r', linewidth=3, label='No Coulomb')

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
        p, rr, t = CoulombDiffusionEquation(qinit=q, length=length)
        dt = t[1] - t[0]
        tIndex = int(timeslice/dt)

        ax.plot(rr*1e6, p[tIndex,:]/1e6, linewidth=3, color=bmap[i])

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
        p, rr, t = GroomDiffusionEquation(length=l*1e-6)
        dt = t[1] - t[0]
        tIndex = int(timeslice/dt)

        ax.plot(rr*1e6, p[tIndex,:]/1e6, linewidth=3, color=bmap[i])

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
    E = np.linspace(0, 15000, 10)
    w = 3.8
    drho = 0.5e-6
    dz = 0.5e-6
    ncharge = E/w
    pall = []
    legstring = []
    fwhmMatrix =  []

    for teff in teffArray:
        fwhmArray = []
        for n in ncharge:

            if n == 0:
                p, rr, t = GroomDiffusionEquation(qinit=100*1.6e-19, length=rinit)
                n = 100
            else:
                p, rr, t = PiecewiseDiffusionEquation(teff=teff, qinit=n*1.6e-19, length=rinit)
            dt = np.diff(t)[0]
            print(dt)
            pfinal = p[int(tsample/dt),:]

            # Convert from radial to x
            theta = np.linspace(0, np.pi, 4000)
            pfinalMatrix = np.reshape(np.repeat(pfinal, theta.size), (pfinal.size, theta.size))
            rrmesh, thetamesh = np.meshgrid(rr, theta, indexing='ij')

            px, x = ProjectSphere2Rho(pfinalMatrix, rrmesh, thetamesh, drho, dz)

            fwhm = x[np.nonzero(px < px[0]/2)[0][0]]
            fwhm = np.interp(px[0]/2, np.flipud(px), np.flipud(x))
            fwhmArray.append(fwhm)
            legstring.append('E=%.2f keV, FWHM=%.1f'%(n*w/1000, fwhm*2*1e6))
            pall.append(px/n)
        fwhmMatrix.append(fwhmArray)

    fwhmMatrix = np.array(fwhmMatrix)
    ax.plot(E/1000, fwhmMatrix.T*1e6, linewidth=3)
    ax.legend(['teff=%0.1f ns'%(x*1e9) for x in teffArray], fontsize=16)
    ax.set_ylabel('FWHM ($\mu m$)', fontsize=16)
    ax.set_ylim([0,1.05*fwhmMatrix.max()*1e6])
    ax.set_xlabel('Energy (keV)', fontsize=16)
    ax.set_title('Energy vs FWHM for different $t_{eff}$', fontsize=18)


    # pall = np.array(pall)
    # ax.plot(x[2:]*1e6, pall[:,2:].T, linewidth=2)
    # ax.legend(legstring, fontsize=16)
    # ax.set_ylabel('Charge Density Distribution', fontsize=16)
    # ax.set_xlabel('X ($\mu m$)', fontsize=16)
    # ax.set_title('Charge Distribution from Using Piecewise Diffusion model After 20 ns. teff=%.1f ns, Ri=1$\mu m$'%(teff*1e9), fontsize=18)
    plt.show()

def Sphere2Cyl(f, rr, theta, drho, dz):
    """
    Converts a function in spherical coordinates to cylindrical coordinates. Theta is the polar angle

    f - two dimensional matrix describing the function at each r and theta value
    r, theta - 2d meshgrid of coordinates.

    All mesh grids should be in the {ij} indexing scheme. Ie first dimension is r, second dimension is theta
    """

    # Find the cylindrical coordinates of all spherical pairs
    rho = rr * np.sin(theta)
    z = rr * np.cos(theta)
    dr = np.diff(rr, axis=0)[0][0]
    fWeight = f*rr*dr

    # Calculate the number of bins to histogram given the range and the wanted spacing of the axis
    nbinsRho = int((rho.max() - rho.min()) / drho)
    nbinsZ = int((z.max() - z.min()) / dz)

    # Find the value of the cylindrical function by histogramming the data we have
    fCyl, rhoEdges, zEdges = np.histogram2d(rho.flatten(), z.flatten(), bins=(nbinsRho, nbinsZ), weights=fWeight.flatten())

    rhoAx = rhoEdges[1:] - np.diff(rhoEdges)/2
    zAx = zEdges[1:] - np.diff(zEdges)/2

    return fCyl, rhoAx, zAx

def ProjectSphere2Rho(f, rr, theta, drho, dz):
    """
    Takes the function in spherical coordinates and projects it onto the xy plane (ie rho). Assumes azimuthal symmetry is maintained
    """

    fCyl, rhoAx, _ = Sphere2Cyl(f, rr, theta, drho, dz)

    # sum over the z axis. Effect to integrating out the variable and projecting it on the rho plane
    fRho = np.sum(fCyl, axis=1)

    return fRho, rhoAx

def testCylConversion():

    r = np.linspace(0, 1, 2000)
    theta = np.arange(0, np.pi, np.pi/4000)

    rmesh, thetamesh = np.meshgrid(r, theta, indexing='ij')

    f = np.reshape(np.repeat(np.linspace(0,1,r.size)**0, theta.size), (rmesh.shape))

    drho = 0.01
    dz = 0.01

    fCyl, rhoAx, zAx = Sphere2Cyl(f, rmesh, thetamesh, drho, dz)

    fig, ax = plt.subplots()

    y, x = ProjectSphere2Rho(f, rmesh, thetamesh, drho, dz)
    ax.plot(x,y)
    plt.show()

if __name__ == '__main__':

    # ComparisonOfEandNoEField()
    # EnergyDependenceOfCoulomb()
    # InitialRadiusDependence()
    # PiecewiseDiffusionEquation(qinit=1)
    # testPiecewiseSolution()
    # testCylConversion()
    energyDependencePiecewise()
    plt.show()




