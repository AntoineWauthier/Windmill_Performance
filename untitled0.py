import numpy as np
import matplotlib.pyplot as plt
import math
import scipy

global R
global H
global c
global rho
global V0
global Pdispo    
R = 0.15
H = 0.4
c = 0.25
rho = 1.2
V0 = 8
Pdispo = (rho*R*H*V0**3)/2

listeAngle = [0,3,5,7,10,15,20,25,30,60,90]
listeKx = [0.006,0.008,0.009,0.011,0.013,0.017,0.024,0.033,0.038,0.065,0.081]
listeKy = [0.022,0.048,0.064,0.078,0.085,0.09,0.094,0.084,0.07,0.039,0] #valeur Kx et Ky d'Eiffel

rho = 1.2 #masse volumique de l'air

listeCx = []
listeCz = []

for x in listeKx:
    listeCx.append(2*x/rho)
for x in listeKy:
    listeCz.append(2*x/rho)

angles = np.array(listeAngle)
Cx = np.array(listeCx)
Cz = np.array(listeCz)

deg = 3

coeffX = np.polyfit(angles,Cx,deg)
coeffZ = np.polyfit(angles,Cz,deg)
    
def formuleCx(angle,deg):
    res = 0
    for k in range(deg+1):
        res+=coeffX[k]*angle**(deg-k)
    return res

def formuleCz(angle,deg):
    res = 0
    for k in range(deg+1):
        res+=coeffZ[k]*angle**(deg-k)
    return res

continuumAngle = np.linspace(0,90,10001)
continuumCx = np.zeros(10001)
continuumCz = np.zeros(10001)

for k in range(10001): #actualisation de chaque élément des listes en calculant
    continuumCx[k] = formuleCx(continuumAngle[k],deg) 
    continuumCz[k] = formuleCz(continuumAngle[k],deg)
    
def alpha(theta,vSpec,a): #calcul de l'angle alpha
    return np.arctan(np.sin(theta)/((vSpec/(1-a))+np.cos(theta)))

