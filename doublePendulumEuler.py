import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#define constants
m1 = 10
m2 = 1
L1 = 1
L2 = 2

g=9.81
dt = 0.01
tmax = 30
t=0
times=np.linspace(0, tmax, int(tmax/dt))

theta1 = np.pi/1
theta2 = 0.1
y1 = 0.0 #omega1
y2 = 0.0

def sin(x):
    return(np.sin(x))

def cos(x):
    return(np.cos(x))

def euler(theta1, theta2, y1, y2):
    theta1 = theta1 + dt*y1
    theta2 = theta2 + dt*y2

    c=cos(theta1-theta2)
    s=sin(theta1-theta2)

    y1 = y1 + dt * (( (m2*g*sin(theta2)*c) - (m2*s*(L1*y1**2*c + L2*y2**2)) - ((m1+m2)*g*sin(theta1))) / (L1* (m1 + m2*s**2)))
    y2 = y2 + dt * (( (m1+m2)*(L1*y1**2*s - g*sin(theta2) + g*sin(theta1)*c) + m2*L2*y2**2*s*c) / ( L2*(m1 + m2*s**2) ))
    return (theta1, theta2, y1, y2)

theta1_vals = np.zeros(int(tmax/dt))
omega1_vals = np.zeros(int(tmax/dt))
theta2_vals = np.zeros(int(tmax/dt))
omega2_vals = np.zeros(int(tmax/dt))

for i in range(int(tmax/dt)):
    (theta1, theta2, y1, y2) = euler(theta1, theta2, y1, y2)
    theta1_vals[i] = theta1
    theta2_vals[i] = theta2
    omega1_vals[i] = y1
    omega2_vals[i] = y2
    t += dt


#convert to cartesian coords
x1 = L1*sin(theta1_vals)
y1 = -L1*cos(theta1_vals)
x2 = x1+L2*sin(theta2_vals)
y2 = y1-L2*cos(theta2_vals)



fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-1.2*(L1+L2), 1.2*(L1+L2)), ylim=(-1.2*(L1+L2), 1.2*(L1+L2)))
ax.set_aspect('equal')
ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)


def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text


def animate(i):
    thisx = [0, x1[i], x2[i]]
    thisy = [0, y1[i], y2[i]]

    line.set_data(thisx, thisy)
    time_text.set_text(time_template % (i*dt))
    return line, time_text


ani = animation.FuncAnimation(fig, animate, range(1, len(x1)),
                              interval=dt*1, blit=True, init_func=init)
plt.show()
