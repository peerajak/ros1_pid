#!/usr/bin/env python3
import numpy as np
from matplotlib import pyplot as plt 

def D_robot(m1,m2,l1,l2,lc1,lc2,Izz1,Izz2,th):
  D=np.zeros((2,2))
  t1=th[0,0]
  t2=th[1,0]
  D[0,0]=Izz1+Izz2+lc1**2*m1+l1**2*m2+2*np.cos(t2)*l1*lc2*m2+lc2**2*m2
  D[0,1]=Izz2+np.cos(t2)*l1*lc2*m2+lc2**2*m2
  D[1,0]=Izz2+np.cos(t2)*l1*lc2*m2+lc2**2*m2
  D[1,1]=Izz2+lc2**2*m2
  return D

def C_robot(m1,m2,l1,l2,lc1,lc2,Izz1,Izz2,th,dth):
  C=np.zeros((2,2))
  t1=th[0,0]
  t2=th[1,0]

  dt1=dth[0,0]
  dt2=dth[1,0]
  
  C[0,0]=-2 * np.sin(t2) * l1 * lc2 * m2 * dt2
  C[0,1]=-np.sin(t2) * l1 * lc2 * m2 * dt2
  C[1,0]= np.sin(t2) * l1 * lc2 * m2 * dt1;
  C[1,1]= 0.0

  return C

def g_robot(m1,m2,l1,l2,lc1,lc2,Izz1,Izz2,th,grav):
  gravity=np.zeros((2,1))
  
  t1=th[0,0]
  t2=th[1,0]

  gravity[0,0]=0.0
  gravity[1,0]=0.0

  return gravity

  dt=0.01 #step time for the simulation

m1=1.0 #mass link 1
m2=1.0 #mass link 2

l1=1.0 #length link 1
l2=1.0 #length link 2

lc1=0.5 #distance to the center of mass link 1
lc2=0.5 #distance to the center of mass link 2

Izz1=(1.0/12.0)*(0.1*0.01+l1**2) #inertia link 1
Izz2=(1.0/12.0)*(0.1*0.01+l2**2) #inertia link 2

grav=9.81 #gravity

Jeff1=Izz1+m1*lc1**2 #effective inertia of link 1 in the motor 
Jeff2=Izz2+m2*lc2**2 #effective inertia of link 2 in the motor

K1=0.1 #voltage vs torque, motor 1
K2=0.1 #voltage vs torque, motor 2

Beff1=0.001 #angular velocity, friction gain, motor 1
Beff2=0.001 #angular velocity, friction gain, motor 2

red1=1 #reduction relation, link-motor 1 (we assume motor is connected directly to motor)
red2=1 #reduction relation, link-motor 2 (we assume motor is connected directly to motor)

  # inertia matrix
  D=D_robot(m1,m2,l1,l2,lc1,lc2,Izz1,Izz2,q)
  # coriolis term
  C=C_robot(m1,m2,l1,l2,lc1,lc2,Izz1,Izz2,q,dq)
  # gravity vector
  g=g_robot(m1,m2,l1,l2,lc1,lc2,Izz1,Izz2,q,grav)
  # matrix of inertias
  JM=np.diagflat(np.array([(Jeff1)/(red1**2),(Jeff2)/(red2**2)]))
  # matrix of friction
  BM=np.diagflat(np.array([(Beff1)/(red1**2),(Beff2)/(red2**2)]))

  # dynamic simulation, equation (1)
  ddq=np.linalg.inv(D+JM)@(input-C@dq-BM@q-g)

  # update angular velocities values using the derivative approximation
  dt = 0.01
  dq=dq+dt*ddq 
  # update theta angles using the derivative approximation
  q=q+dt*dq