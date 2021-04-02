# DoublePendulum
Uses the coupled second order differential equaitons obtained from the Euler-Lagrange equation to simulate the motion of a double pendulum.

`doublePendulum.py` uses `scipy.integrate.odeint()` to numerically solve the equations, and generates the simulation with `vpython`, and offers graphing from `matplotlib`

`doublePendulumEuler` uses Euler's method for solving differential equations numerically, and generates the simulation using `matplotlib`'s `animation.FuncAnimation()` function
