import numpy as np
import matplotlib.pyplot as plt
import scipy

global R
global H
global c
R = 0.15
H = 0.35
c = 0.25
rho = 1.2


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

ventMaxVecteur = np.array([0,8])
ventMax = 8

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
        
        plt.scatter(angles,expCx, color = "orange")
        plt.scatter(angles,expCz, color = "b")
        plt.plot(continuumAngle,continuumCx, label = "Cx", color = "orange")
        plt.xlabel("Angle (degrés)")
        plt.plot(continuumAngle,continuumCz,label = "Cz", color = "b")
        plt.legend()
        plt.show()
    return(continuumAngle,continuumCx,continuumCz)
