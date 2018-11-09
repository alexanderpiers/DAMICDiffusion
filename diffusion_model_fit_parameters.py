#Diffusion Model Fit Parameters
import sys
sys.path.append(r'C:\Users\alexp\Documents\UW\Research\DAMIC\diffusion_model')
import diffusion_model_numeric_solutions as dm
import numpy as np
import matplotlib.pyplot as plt

def generateData(x, pixelwidth=15.e-6):
    y = (0.0035*x/1.e3+0.5)*pixelwidth
    return y

def fitParameters():
    tdriftGuess = 11.0e-9
    teffGuess = 7.8e-9

    # define parameters
    rinit = 0.5e-6

    tend = 0.03e-6 # 50 ns
    nt = 300

    w = 3.8

    # energy range
    energy = np.linspace(1e3,20e3, 8)
    qinit =  energy*1.6e-19/w

    error = 1e9
    errorflag = True
    dtdrift = 1.e-9
    dteff = 0.4e-9

    theta = np.linspace(0, np.pi, 4000)
    drho = 0.5e-6
    dz = 0.5e-6

    gamma = 50

    tDriftArray = []
    tEndArray = []

    while errorflag:

        sigmax=[]
        sigmaxP = []
        sigmaxN = []
        sigmaEffP = []
        sigmaEffN = []
        for q in qinit:
            p, rr, t = dm.PiecewiseDiffusionEquation(teff=teffGuess, qinit=q, length=rinit, rmax=30.e-6, nr=500, tend=tend, nt=nt)
            dt = np.diff(t)[0]
            ptime = p[int(tdriftGuess/dt),:]
            ptimeP = p[int((tdriftGuess+dtdrift)/dt),:]
            ptimeN = p[int((tdriftGuess-dtdrift)/dt),:]
            ptimeMat = np.reshape(np.repeat(ptime, theta.size), (ptime.size, theta.size))
            ptimeMatP = np.reshape(np.repeat(ptimeP, theta.size), (ptimeP.size, theta.size))
            ptimeMatN = np.reshape(np.repeat(ptimeN, theta.size), (ptimeN.size, theta.size))

            peffP, _, _ = dm.PiecewiseDiffusionEquation(teff=(teffGuess+dteff), qinit=q, length=rinit, rmax=30.e-6, nr=500, tend=tend, nt=nt)
            peffN, _, _ = dm.PiecewiseDiffusionEquation(teff=(teffGuess-dteff), qinit=q, length=rinit, rmax=30.e-6, nr=500, tend=tend, nt=nt)
            peffP = peffP[int(tdriftGuess/dt),:]
            peffN = peffN[int(tdriftGuess/dt),:]

            peffMatP = np.reshape(np.repeat(peffP, theta.size), (peffP.size, theta.size))
            peffMatN = np.reshape(np.repeat(peffN, theta.size), (peffN.size, theta.size))

            rrmesh, thetamesh = np.meshgrid(rr, theta, indexing='ij')
            sigmax.append(halfMax(dm.ProjectSphere2Rho(ptimeMat, rrmesh, thetamesh, drho, dz)))
            sigmaxP.append(halfMax(dm.ProjectSphere2Rho(ptimeMatP, rrmesh, thetamesh, drho, dz)))
            sigmaxN.append(halfMax(dm.ProjectSphere2Rho(ptimeMatN, rrmesh, thetamesh, drho, dz)))
            sigmaEffP.append(halfMax(dm.ProjectSphere2Rho(peffMatP, rrmesh, thetamesh, drho, dz)))
            sigmaEffN.append(halfMax(dm.ProjectSphere2Rho(peffMatN, rrmesh, thetamesh, drho, dz)))

        realspread = generateData(energy)
        e = np.sum((realspread-sigmax)**2)
        eP = np.sum((realspread-sigmaxP)**2)
        eN = np.sum((realspread-sigmaxN)**2)
        eEffP = np.sum((realspread-sigmaEffP)**2)
        eEffN = np.sum((realspread-sigmaEffN)**2)


        if e > error:
            errorflag = False
            pass

        error = e

        timeslope = (eP - eN) / (2 * dtdrift)
        effslope = (eEffP - eEffN) / (2 * dteff)

        teffGuess -= np.sign(effslope)*np.ceil(np.abs(gamma*effslope*dteff/dt))*dt
        tdriftGuess -= np.sign(timeslope)*np.ceil(np.abs(gamma*timeslope*dtdrift/dt))*dt

        tDriftArray.append(tdriftGuess)
        tEndArray.append(teffGuess)

        print(error)
        print(teffGuess)
        print(tdriftGuess)

    return teffGuess, tdriftGuess, [tDriftArray, tEndArray]

def testDiffusionOnly(tdrift):
    rinit = 0.5e-6
    tend = 0.02e-6
    nt = 400

    w = 3.8
    qinit = 100*1.6e-19

    theta = np.linspace(0,np.pi,4000)
    drho = 0.1e-6
    dz = 0.1e-6

    p, rr, t = dm.GroomDiffusionEquation(qinit=qinit, length=rinit, rmax=50.e-6, nr=1000, tend=tend, nt=nt)
    dt = np.diff[0]
    ptime = p[int(tdrift/dt),:]

    rrmesh, thetamesh = np.meshgrid(rr, theta, indexing='ij')

def fwhmvstime():

    energy = 1000 # 5 kev

    # define parameters
    rinit = 0.5e-6
    tend = 0.03e-6
    nt = 600

    tsample = np.linspace(0.1, 12, 25)*1e-9

    w = 3.8
    qinit = energy/w*1.6e-19
    teff = np.linspace(1,8,8)*1e-9
    teff = np.concatenate((teff, [20e-9]))

    theta = np.linspace(0, np.pi, 4000)
    drho = 0.5e-6
    dz = 0.5e-6

    allt = []
    for te in teff:
        p, rr, t = dm.PiecewiseDiffusionEquation(teff=te, qinit=qinit, length=rinit, rmax=25.e-6, nr=500, tend=tend, nt=nt)
        dt = np.diff(t)[0]
        ptime = p[(tsample/dt).astype('int'),:]
        rrmesh, thetamesh = np.meshgrid(rr, theta, indexing='ij')
        sigma = []
        for i in range(tsample.size):
            pMat = np.reshape(np.repeat(ptime[i,:], theta.size), (rr.size, theta.size))
            sigma.append(halfMax(dm.ProjectSphere2Rho(pMat, rrmesh, thetamesh, drho, dz)))

        allt.append(sigma)

    fig, ax = plt.subplots()
    ax.plot(tsample*1e9, np.array(allt).T*1e6, linewidth=2)
    legstr = ['$t_{eff}$=%i'%(x*1e9) for x in teff]
    legstr.pop()
    legstr.append('$t_{eff}=\infty$')
    ax.legend(legstr, fontsize=16)
    ax.set_ylabel('$\sigma_x$ ($\mu m$)', fontsize=16)
    ax.set_xlabel('Drift Time (ns)', fontsize=16)
    ax.set_title('$\sigma_x$ vs time for a %i keV deposition'%(energy/1e3), fontsize=18)
    plt.show()

def fwhmvsteff():


    # define parameters
    rinit = 0.5e-6
    tend = 0.03e-6
    nt = 800

    tsample = 12.e-9

    w = 3.8
    theta = np.linspace(0, np.pi, 4000)
    drho = 0.5e-6
    dz = 0.5e-6

    energy = np.linspace(0,15, 5)*1e3
    qinit = energy/w*1.6e-19
    teff = np.linspace(1, 10, 10)*1e-9

    sigmaAll=[]
    for q in qinit:
        sigma = []
        for te in teff:
            if q == 0:
                p, rr, t = dm.GroomDiffusionEquation(qinit=100*1.6e-19, length=rinit, rmax=50.e-6, nr=500, tend=tend, nt=nt)
            else:
                p, rr, t = dm.PiecewiseDiffusionEquation(teff=te, qinit=q, length=rinit, rmax=50.e-6, nr=500, tend=tend, nt=nt)

            dt = np.diff(t)[0]
            ptime = p[int(tsample/dt),:]
            rrmesh, thetamesh = np.meshgrid(rr, theta, indexing='ij')
            pMat = np.reshape(np.repeat(ptime, theta.size), (rr.size, theta.size))
            sigma.append(halfMax(dm.ProjectSphere2Rho(pMat, rrmesh, thetamesh, drho, dz)))

        sigmaAll.append(sigma)

    fig, ax = plt.subplots()
    ax.plot(teff*1e9, np.array(sigmaAll).T*1e6, linewidth=2)
    ax.legend(['Energy = %i keV'%(x/1e3) for x in energy])
    ax.set_xlabel('$t_{eff}$ (ns)', fontsize=16)
    ax.set_ylabel('$\sigma_x$ ($\mu m$)', fontsize=16)
    ax.set_title('$\sigma_x$ vs $t_{eff}$ at sampled time %i ns'%(tsample*1e9), fontsize=18)
    plt.show()








def comparePlot(tdriftGuess=11.1e-9, teffGuess=7.9e-9):

    # define parameters
    rinit = 0.5e-6

    tend = 0.03e-6 # 50 ns
    nt = 300

    w = 3.8

    # energy range
    energy = np.linspace(1e3,100e3, 20)
    energyKnown = np.linspace(1e3, 20e3, 10)
    qinit =  energy*1.6e-19/w

    theta = np.linspace(0, np.pi, 4000)
    drho = 0.5e-6
    dz = 0.5e-6

    sigma = []

    for q in qinit:
        p, rr, t = dm.PiecewiseDiffusionEquation(teff=teffGuess, qinit=q, length=rinit, rmax=50.e-6, nr=500, tend=tend, nt=nt)
        dt = np.diff(t)[0]
        ptime = p[int(tdriftGuess/dt),:]
        ptimeMat = np.reshape(np.repeat(ptime, theta.size), (ptime.size, theta.size))
        rrmesh, thetamesh = np.meshgrid(rr, theta, indexing='ij')
        sigma.append(halfMax(dm.ProjectSphere2Rho(ptimeMat, rrmesh, thetamesh, drho, dz)))

    fig, ax = plt.subplots()
    ax.plot(energyKnown/1e3, generateData(energyKnown)*1e6, '+', markersize=12, linewidth=2, color='black', label='Guillermo Data')
    ax.plot(energy/1e3, np.array(sigma)*1e6, linewidth=2, color='blue', label='Diffusion Model')
    ax.plot(energy[5:]/1e3, generateData(energy[5:])*1e6, '+', markersize=12, linewidth=2, color='red', label='Guillermo Extended')
    # ax.set_ylim([0,12])
    ax.legend(fontsize=16)
    ax.set_xlabel('Energy (keV)', fontsize=16)
    ax.set_ylabel('$\sigma_x$ ($\mu m$)', fontsize=16)
    ax.set_title('Diffusion Drift Parameter Fitting. $t_{eff}$=%.1f ns, $t_{drift}$=%.1f ns'%(teffGuess*1e9, tdriftGuess*1e9), fontsize=18)
    plt.show()



def halfMax(f):
    return 2*np.interp(f[0][0]/2, np.flipud(f[0]), np.flipud(f[1]))/2.355




if __name__ == '__main__':
    # teff, tdrift, array = fitParameters()

    # np.save('./coefficients.npy', np.array(array))
    # print(teff)
    # print(tdrift)
    # comparePlot(tdrift, teff)
    # comparePlot()
    # fwhmvstime()
    fwhmvsteff()
