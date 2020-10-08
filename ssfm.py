# %%
%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack

#CHAOTIC SOLITONS
D = 0.5
gamma1 = 11
gamma2 = 9
Na =gamma1/(6*D)
M = (3*gamma2)/(4*gamma1)
epsilon = 0.14
omega =3
R0 = 0.1
N0 = 8.8
C=Na*(0.85*M)**3
M0 = (C/N0)**(1/3)
# M0 = 0.85*M
G = 8.1
B  = 0.53
F = (C/G)**(1/3)


# Number of samplepoints
N = 2**13
dt = 0.002 # Time Step
finalTime = 200

numberSteps = finalTime/dt

L = 500
# sample spacing (space step)
h = 0.0625

n = np.linspace(-N/2, N/2, N)
x = n*h
k = 2*n*np.pi/L

#u = np.sin(50.0 * 2.0*np.pi*x) + 0.5*np.sin(80.0 * 2.0*np.pi*x)
u = u = 1/(M0+N0*x**2)*np.exp(1j*R0*x**2)
for m in range(finalTime):

    u = np.exp(dt*1j*(-gamma1*(1+epsilon*np.sin(omega*m*dt))*np.absolute(u)+gamma2*(np.absolute(u)*np.absolute(u))))*u # Solve non-linear part of NLSE
    c = np.fft.fftshift(np.fft.fft(u))
    c = np.exp(dt*1j*(-D*k*k))*c
    u = np.fft.ifft(np.fft.fftshift(c))

uAbs = np.absolute(u)
fig, ax = plt.subplots()
#plt.xlim(-N*h/2,N*h/2)
plt.plot(x,uAbs)

plt.show()
# %%
