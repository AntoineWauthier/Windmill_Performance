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
R = 0.4
H = 0.4
c = 0.2
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

for k in range(10001):
    continuumCx[k] = formuleCx(continuumAngle[k],deg)
    continuumCz[k] = formuleCz(continuumAngle[k],deg)

plt.plot(continuumAngle,continuumCx)
plt.plot(continuumAngle,continuumCz)
#plt.plot(continuumCx,continuumCz)
#plt.scatter(Cx,Cz)
plt.show()




def coeffAlpha(theta,a,vSpec):
    return np.arctan(np.sin(theta)/((vSpec/(1-a))+np.cos(theta)))

def puissance(theta,a,vSpec):
    alpha = coeffAlpha(theta,a,vSpec)
    v = (1-a)**2 + 2*vSpec*(1-a)*np.cos(theta)+vSpec**2
    trainee = formuleCx((180*theta/np.pi)%180,deg)*math.cos(alpha+theta/2)
    portance = formuleCz((180*theta/np.pi)%180,deg)*math.sin(alpha+theta/2)
    v*= portance - trainee
    return v

def CpTournante(vSpec,a):
    return c*vSpec*(1/(2*np.pi*R))*(scipy.integrate.quad(puissance,np.pi,3*np.pi/2,args=(a,vSpec))[0])

def CpNormale(vSpec,a):
    v = (1-a)**2 + 2*vSpec*(1-a)+vSpec**2
    alpha = coeffAlpha(np.pi/2,a,vSpec)
    trainee = formuleCx(90,deg)*math.cos(alpha+np.pi)
    return -0.5*c*vSpec*v*trainee/R
    
    
listvSpec= np.linspace(0,10,1001)
listCp = np.zeros(1001)

for a in [1/10,1/8,1/6,1/5,1/4,1/3]:
    a = 1/3
    for k in range(1001):
       listCp[k]=(CpTournante(listvSpec[k],a))/3
    plt.plot(listvSpec,listCp)
    plt.plot(listvSpec,16/27 * np.ones(1001))
plt.show()

def testFonction(f,a,b,coeff,vspec):
    listx = np.linspace(a,b,10000)
    listy = np.zeros(10000)
    for k in range(10000):
        listy[k]=f(listx[k],coeff,vspec)
    plt.plot(listx,listy)
    plt.show()

testFonction(puissance, 0, 5*np.pi,1/3,2)


    

        
        
        
        

