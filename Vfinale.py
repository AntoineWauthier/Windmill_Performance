import numpy as np
import matplotlib.pyplot as plt
import scipy
import math 

rho = 1.2

def fitAerodynamicCoefficients(expValues,iteration,deg,display=True):

    angles = np.array(expValues[0])
    expCx = np.array(expValues[1])
    expCz = np.array(expValues[2])
    
    
    coeffX = np.polyfit(angles,expCx,deg)
    coeffZ = np.polyfit(angles,expCz,deg)

    def computeCx(angle):
        res = 0
        for k in range(deg+1):
            res+=coeffX[k]*angle**(deg-k)
        return res
    
    def computeCz(angle):
        res = 0
        for k in range(deg+1):
            res+=coeffZ[k]*angle**(deg-k)
        if res<0:
            return 0
        return res

    if display == True:      
        
        continuumAngle = np.linspace(0,90,iteration)
        continuumCx = np.zeros(iteration)
        continuumCz = np.zeros(iteration)
            
        for k in range(iteration):
            continuumCx[k] = computeCx(continuumAngle[k]) 
            continuumCz[k] = computeCz(continuumAngle[k])
    
        plt.scatter(angles,expCx, color = "orange")
        plt.scatter(angles,expCz, color = "b")
        plt.plot(continuumAngle,continuumCx, label = "Cx", color = "orange")
        plt.xlabel("Angle (degrés)")
        plt.plot(continuumAngle,continuumCz,label = "Cz", color = "b")
        plt.legend()
        plt.show()
        
    return(computeCx,computeCz)

def computeTipSpeedRatio(angularSpeed,radius,windSpeed):
    return radius*angularSpeed/windSpeed

def convertTSRtoRPM(TSR,radius,windSpeed):
    return 30*TSR/(radius*np.pi)

def convertRadtoDegree(angle):
    return 180*angle/np.pi

def computeDarrieusIncidenceAngle(angle,TRS,inductionFactor):
    return np.arctan(np.sin(angle)/((TRS/(1-inductionFactor))+np.cos(angle)))

def computeIncidenceAngle(angle,TRS,inductionFactor):
    return computeDarrieusIncidenceAngle(angle, TRS, inductionFactor) - angle/2

def computeSquaredSpeed(angle,angularSpeed,radius,windSpeed,inductionFactor):
    
    if np.pi/2 < angle%(2*np.pi)< 3*np.pi/2 : 
        speed = (1 - inductionFactor) * windSpeed + np.cos(angle)*radius*angularSpeed
        return speed**2
    
    else:
        return windSpeed**2 * ( (1-inductionFactor)**2 + 2*computeTipSpeedRatio(angularSpeed,radius,windSpeed)*(1-inductionFactor)*np.cos(angle) + computeTipSpeedRatio(angularSpeed,radius,windSpeed) )

def computeTorque(angle,angularSpeed,radius,height,width,windSpeed,inductionFactor,computeCx,computeCz):
    
    if np.pi/2<angle%(2*np.pi)<3*np.pi/2:
        torque = 0.5*rho*width*height*computeSquaredSpeed(angle,angularSpeed,radius,windSpeed,inductionFactor)*computeCx(90)*(-np.cos(angle))
    
    else:
        incidenceAngle = computeIncidenceAngle(angle,computeTipSpeedRatio(angularSpeed,radius,windSpeed),inductionFactor)
        torque = 0.5*rho*width*height*computeSquaredSpeed(angle,angularSpeed,radius,windSpeed,inductionFactor)*(-computeCx(convertRadtoDegree(incidenceAngle%(np.pi/2)))*np.cos(incidenceAngle%np.pi)+computeCz(convertRadtoDegree(incidenceAngle%(np.pi/2)))*np.sin(incidenceAngle%np.pi))
    return torque

def computePower(angle,angularSpeed,radius,height,width,windSpeed,inductionFactor,computeCx,computeCz):
    return angularSpeed * computeTorque(angle, angularSpeed, radius, height, width, windSpeed, inductionFactor, computeCx, computeCz)

def computePowerOnFullRotation(angularSpeed,radius,height,width,windSpeed,inductionFactor,numberOfBlades,computeCx,computeCz):
    return numberOfBlades*angularSpeed/(2*np.pi)*scipy.integrate.quad(computeTorque,0,2*np.pi,args = (angularSpeed,radius,height,width,inductionFactor,computeCx,computeCz))[0]

def computeTimeEvolution(radius,height,width,momentOfInertia,windSpeed,inductionFactor,numberOfBlades,Tmax,dt,computeCx,computeCz,display = True):
    
    Nit = int(Tmax//dt)
    bladeAngleList = np.linspace(0,2*np.pi,numberOfBlades+1)[:-1]
    bladeAngleStorage = []
    angularSpeed = 0
    angularSpeedList = []
    torqueList = []
    timeList = np.linspace(0,Tmax,Nit)

    for k in range(Nit):
        k=0
        bladeAngleStorage.append(bladeAngleList) 
        angularSpeedList.append(angularSpeed)
        torque = 0
        listT = []
        for angle in bladeAngleList:
            t = computeTorque(angle,angularSpeed,radius,height,width,windSpeed,inductionFactor,computeCx,computeCz)
            listT.append(t)
            torque+=t
        torqueList.append(listT)
        bladeAngleList += dt*angularSpeed*np.ones(numberOfBlades)
        angularSpeed += dt*torque/momentOfInertia
        
    if display == True:
        
        plt.plot(timeList,torqueList)
        plt.show()
        
    return timeList,angularSpeedList

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

res = fitAerodynamicCoefficients(essaiRapport3, 1000, 3,display = False)
computeTimeEvolution(0.15,0.4,0.25,1,8,0.3,3,3,0.001,res[0],res[1])

omega = np.linspace(0,200,10)
theta = np.linspace(0,2*np.pi,1000)

for o in omega:
    speed  = []
    for t in theta:
        speed.append(computeSquaredSpeed(t, o, 0.15, 8, 0.3))
    plt.plot(theta,speed)
plt.show()
