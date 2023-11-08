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

def norme(vecteur):
    return math.sqrt(vecteur[0]**2+vecteur[1]**2)

def vitesseRelative(theta,omega):
    return ventMax+R*omega*np.array([-np.sin(theta),np.cos(theta)])

def forceTrainee(theta,omega):
    return -vitesseRelative(theta,omega)*norme(vitesseRelative(theta,omega))*rho*H*c/2

def forcePortance(theta,omega):
    return 

Pdispo = (rho*0.4*0.4*ventMax**3)/2

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
    
def alpha0(theta,vSpec,a): #calcul de l'angle alpha
    return np.arctan(np.sin(theta)/((vSpec/(1-a))+np.cos(theta)))

def alpha(theta,vSpec,a):
    return alpha0(theta,vSpec,a)-theta/2

J=1

def CoupleResistif(omega):
    return 0.02*np.sqrt(100*30*omega/(2*np.pi))

def evolution(valeurExperimentale,a,deg,nom):
    
    res = valeurTheorique(valeurExperimentale,deg)
    Cx = res[1]
    Cz = res[2]
    
    def vSpec(omega):
        return R*omega/ventMax
    
    def vitesseCarree(theta,omega):
        if np.pi/2<theta%(2*np.pi)<3*np.pi/2:
            v = (1-a)*ventMax - np.cos(theta)*R*omega
            return (v*abs(v))
        else:
            return (ventMax**2)*((1-a)**2 + 2*vSpec(omega)*(1-a)*np.cos(theta) + vSpec(omega)**2)
    
    def Puissance(omega):
        return 3*omega/(2*np.pi)*scipy.integrate.quad(Couple,0,2*np.pi,args = (omega))[0]
        
    def Couple(theta,omega):
        if np.pi/2<theta%(2*np.pi)<3*np.pi/2:
            C=vitesseCarree(theta,omega)*rho*0.5*0.25*0.35*Cx[-1]*(-np.cos(theta))
        else:
            incidence = alpha(theta, R*omega/ventMax, a)
            indiceCoeffs = int(((incidence%(np.pi/2)*180/np.pi))//10001)
            
            C=rho*0.5*0.25*0.35*vitesseCarree(theta,omega)*(-Cx[indiceCoeffs]*np.cos(incidence%np.pi)+Cz[indiceCoeffs]*np.sin(incidence%np.pi))
        return C
    
    # N = 10001
    # omega = 0
    # t = 0
    # dt = 160/(N-1)
    # listeT = np.linspace(0,160,N)
    # listeOmega = np.zeros(N)
    # listeCouple = np.zeros(N)
    # listePuissance = np.zeros(N)
    # listeThetas = np.zeros((3,N))
    
    # listeCouples = np.zeros((3,N))
    
    # TrioCouple = [0,0,0]
    # theta = np.array([0,np.pi/3,2*np.pi/3])
    
    # listeThetas[0][0] = 0
    # listeThetas[1][0] = np.pi/3
    # listeThetas[2][0] = 2*np.pi/3
    
    # for k in range(1,N):
    #     for i in range(3):
    #         TrioCouple[i] = Couple(theta[i],omega,t)
    #         listeCouples[i][k] = TrioCouple[i]

    #     C = sum(TrioCouple)
    #     listeCouple[k] = C
    
    #     theta+=omega*dt*np.ones(3)
        
    #     for i in range(3):
    #         listeThetas[i][k] = theta[i]
    #         listeThetas[i][k] = theta[i]
    #         listeThetas[i][k] = theta[i]
    #     t+=dt
    #     omega += dt*C/J
    #     listeOmega[k]=omega
    #     listePuissance[k] = omega*C
        
    # plt.plot(listeT,60*listeOmega/(2*np.pi))
    # plt.xlabel("Temps (s)")
    # plt.ylabel("Vitesse de rotation (rpm)")
    # plt.show()
    
    # for i in range(1):  
    #     plt.plot(listeThetas[0],listeCouples[i])
    # for i in range(20):
    #     plt.plot([(2*i+1)*np.pi/2,(2*i+1)*np.pi/2],[-2,3])
     #plt.show()
    facteurPerte = 0.5
    #plt.plot((R/ventMax)*listeOmega[0:-N//1000],(1/Pdispo)*listePuissance[0:-N//1000]-facteurPerte*((R/ventMax)*listeOmega[0:-N//1000])**2, label=nom)
    #plt.plot(listeT,listePuissance)
    OMEGA = np.linspace(0,40,4001)
    CPfluide = np.zeros(4001)
    CPreel = np.zeros(4001)
    for k in range(4001):
        CPfluide[k] = Puissance(OMEGA[k])/Pdispo
        CPreel[k] = CPfluide[k] - facteurPerte*vSpec(OMEGA[k])**2
    plt.plot(vSpec(OMEGA),CPfluide, label ="CpFluide")
    plt.plot(vSpec(OMEGA),CPreel, label ="CpReel")


    

def moyenneur(liste,tranche,longueur):
    
    listemoyenne = np.zeros(longueur)
    
    for i in range(longueur):
        moyenne = 0
        if i<(tranche//2+1):
            for j in range(i):
                moyenne+=liste[j]
            listemoyenne[i]=moyenne/(i+1)
        elif i>longueur-(tranche//2+1):
            for j in range(longueur-i):
                moyenne+=liste[longueur-j-1]
            listemoyenne[i]=moyenne/(longueur-i)
        else:
            for j in range(tranche):
                moyenne+=liste[j+i-tranche//2]
            listemoyenne[i]=moyenne/tranche
    return listemoyenne
    
# evolution(Eiffel,1/3,3,"Eiffel")
# evolution(essaiRapport3,1/3,3,"Rapport 3")
# evolution(essaiRapport10,1/3,3,"Rapport 10")
evolution(essaiRapport20,1/3,3,"Rapport 20")
plt.title("Cp de la pâle réalisée")
plt.xlabel("Vitesse spécifique")
plt.ylabel("Cp")
plt.legend()

#plt.plot(listeT,listeCouple)

plt.show()