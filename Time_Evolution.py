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
