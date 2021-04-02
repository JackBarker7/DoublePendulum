import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
from vpython import *
import time
#define constants
m1 = 2
m2 = 1
L1 = 1
L2 = 2

g=9.81
dt = 0.001
tmax = 30

theta1_0 = np.pi/1
theta2_0 = 0.1
omega1_0 = 0.0
omega2_0 = 0.0

initial = (theta1_0, omega1_0, theta2_0, omega2_0)

#set up scene
mainScene = canvas(title = "Double Pendulum", height = 700, width = 1000, background = color.green, align = "left")
mainScene.camera.pos = vector(0, -(L1+L2)/2, 0)

def randomise():
    '''generates a pendulum with randomised starting conditions'''
    global m1, m2, L1, L2
    m1 = np.random.uniform(low = 0.1, high = 10)
    m2 = np.random.uniform(low = 0.1, high = 10)
    L1 = np.random.uniform(low = 0.1, high = 10)
    L2 = np.random.uniform(low = 0.1, high = 10)


def sin(x):
    return(np.sin(x))

def cos(x):
    return(np.cos(x))

t=np.linspace(0, tmax, int(tmax/dt)) #array with all times

def int_dpendulum_sim(initial, t):

    theta1, omega1, theta2, omega2 = initial

    c = cos(theta1-theta2)
    s = sin(theta1-theta2)

    theta1_dot = omega1
    omega1_dot = ( (m2*g*sin(theta2)*c) - (m2*s*(L1*omega1**2*c + L2*omega2**2)) - ((m1+m2)*g*sin(theta1))) / (L1* (m1 + m2*s**2))

    theta2_dot = omega2
    omega2_dot = ( (m1+m2)*(L1*omega1**2*s - g*sin(theta2) + g*sin(theta1)*c) + m2*L2*omega2**2*s*c) / ( L2*(m1 + m2*s**2) )

    return theta1_dot, omega1_dot, theta2_dot, omega2_dot

#perform numerical integration
y = integrate.odeint(int_dpendulum_sim, initial, t)

#Get omegas and thetas as functions of time
theta1 = np.zeros(int(tmax/dt))
omega1 = np.zeros(int(tmax/dt))
theta2 = np.zeros(int(tmax/dt))
omega2 = np.zeros(int(tmax/dt))
for i in range(int(tmax/dt)):
    theta1[i] = y[i][0]
    omega1[i] = y[i][1]
    theta2[i] = y[i][2]
    omega2[i] = y[i][3]

#convert to cartesian coords
x1 = L1*sin(theta1)
y1 = -L1*cos(theta1)
x2 = x1+L2*sin(theta2)
y2 = y1-L2*cos(theta2)
print(x1)

def graphing():
    plt.plot(t, x1)
    plt.show()

def animate():
    gd = graph(title = "energy graph", 
              ytitle = "Energy (J)", 
              xtitle = "Time (s)",
              align = "right")
    gpe_graph = gcurve(color = color.red, label = "GPE")
    ke_graph = gcurve(color = color.blue, label = "KE")
    total_graph = gcurve(color = color.green, label = "Total Energy")

    rod1 = curve(vector( x1[0], y1[0], 0 ), vector(0,0,0), color = color.red)
    rod2 = curve(vector( x2[0], y2[0], 0 ), vector( x1[0], y1[0], 0 ), color = color.blue)
    ball1 = sphere(pos = vector( x1[0], y1[0], 0 ),
                  radius = 0.1, color = color.red, make_trail = True, retain = 50)

    ball2 = sphere(pos = vector( x2[0], y2[0], 0 ),
                  radius = 0.1, color = color.blue, make_trail = True, retain = 50)
    
    for i in range(int(tmax/dt)):

        rod1.modify(0, pos = vector( x1[i], y1[i], 0 ))
        rod2.modify(0, pos = vector( x2[i], y2[i], 0 ))
        rod2.modify(1, pos = vector( x1[i], y1[i], 0 ))
        ball1.pos = vector( x1[i], y1[i], 0 )
        ball2.pos = vector( x2[i], y2[i], 0 )

        #plot graphs:
        GPE = m1*g*(L1+L2+y1[i]) + m2*g*(L1+L2+y2[i])
        gpe_graph.plot(t[i], GPE)
        KE = 0.5*( m1*L1**2*omega1[i]**2 + m2*( L1**2*omega1[i]**2 + L2**2*omega2[i]**2 + 2*L1*L2*omega1[i]*omega2[i]*cos(theta1[i]-theta2[i])) )
        ke_graph.plot(t[i], KE)
        total_graph.plot(t[i], KE+GPE)
        time.sleep(dt/100)

        



animate()