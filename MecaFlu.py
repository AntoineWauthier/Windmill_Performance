import numpy as np
import matplotlib.pyplot as plt
import scipy

global R
global H
global c
R = 0.15
H = 0.4
c = 0.25

essaiRapport3 = [[0,3,6,9,9.5,10,10.5,11,13,15,20,25,30,90],
                 [0.040,0.040,0.061,0.099,0.111,0.122,0.133,0.144,0.194,0.244,0.327,0.417,0.506,1.315],
                 [0,0.158,0.327,0.485,0.521,0.539,0.557,0.575,0.647,0.675,0.618,0.639,0.661,0]]

essaiRapport10 = [[0,3,6,9,11,11.5,12,12.5,13,15,20,25,30,90],
                  [0.021,0.031,0.050,0.079,0.108,0.126,0.129,0.151,0.169,0.219,0.298,0.384,0.471,1.331],
                  [0,0.216,0.359,0.485,0.553,0.564,0.575,0.582,0.589,0.611,0.578,0.575,0.596,0]]

essaiRapport20=[[0,    3,    6,    9,    12,   12.5, 13,   15,   20,   25,   30,   90],
                [0.023,0.032,0.050,0.079,0.119,0.129,0.147,0.201,0.291,0.366,0.453,1.350],
                [0,0.251,0.377,0.503,0.575,0.589,0.603,0.629,0.585,0.571,0.575,0]]

Eiffel = [[0,3,5,7,10,15,20,25,30,60,90],
          [0.006,0.008,0.009,0.011,0.013,0.017,0.024,0.033,0.038,0.065,0.081],
          [0.022,0.048,0.064,0.078,0.085,0.09,0.094,0.084,0.07,0.039,0]]

def valeurTheorique(valeurExperimentale,deg,affichage=True):

    angles = np.array(valeurExperimentale[0])
    expCx = np.array(valeurExperimentale[1])
    expCz = np.array(valeurExperimentale[2])
    
    
    coeffX = np.polyfit(angles,expCx,deg)
    coeffZ = np.polyfit(angles,expCz,deg)

    def formuleCx(angle):
        res = 0
        for k in range(deg+1):
            res+=coeffX[k]*angle**(deg-k)
        return res
    
    def formuleCz(angle):
        res = 0
        for k in range(deg+1):
            res+=coeffZ[k]*angle**(deg-k)
        if res<0:
            return 0
        return res

    if affichage == True:
        
        continuumAngle = np.linspace(0,90,10001)
        continuumCx = np.zeros(10001)
        continuumCz = np.zeros(10001)
        
        for k in range(10001): #actualisation de chaque élément des listes en calculant
            continuumCx[k] = formuleCx(continuumAngle[k]) 
            continuumCz[k] = formuleCz(continuumAngle[k])
        
        plt.plot(continuumAngle,continuumCx)
        plt.plot(continuumAngle,continuumCz)
        plt.show()
    return(continuumAngle,continuumCx,continuumCz)
    
def alpha0(theta,vSpec,a): #calcul de l'angle alpha
    return np.arctan(np.sin(theta)/((vSpec/(1-a))+np.cos(theta)))

def alpha(theta,vSpec,a):
    return alpha0(theta,vSpec,a)-theta/2

def calculCpPaleUnique(valeurExperimentale,a,deg):
    res = valeurTheorique(valeurExperimentale,deg)
    Cx = res[1]
    Cz = res[2]
    def Couple(theta,vSpec,a):
        if np.pi/2<theta%(2*np.pi)<3*np.pi/2:
            C=((1-a)**2-2*vSpec*(1-a)+vSpec**2)*Cx[-1]
        else:
             incidence = alpha(theta, vSpec, a)
             indiceCoeffs = int(((incidence%(np.pi/2)*180/np.pi))//10001)
             C=((1-a)**2+2*vSpec*(1-a)*np.cos(theta)+vSpec**2)*(-Cx[indiceCoeffs]*np.cos(incidence)+Cz[indiceCoeffs]*np.sin(incidence))
        return C
    
    Cp = np.zeros(1001)
    listeVspec = np.linspace(0,0.7,1001)
    
    for k in range(1001):
        vSpec = listeVspec[k]
        integrale = scipy.integrate.quad(Couple,0,2*np.pi/2,args=(vSpec,a))
        Cp[k] = c*vSpec*integrale[0]/(R*4*np.pi)
    plt.plot(listeVspec,Cp)
    #plt.plot(listeVspec,16/27 * np.ones(1001))
    #plt.show()
    return(detectionMaximum(listeVspec,Cp,0.6))

def detectionMaximum(listeV,listeCp,Vlim):
    maxi = listeCp[0]
    V = listeV[0]
    k = 0
    while listeV[k]<Vlim:
        k+=1
        if listeCp[k]>maxi:
            maxi = listeCp[k]
            V = listeV[k]
    return V,maxi
            
    
    
    
