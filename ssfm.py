# %%
#%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack
from mpl_toolkits.mplot3d import Axes3D


epsilon = 0.1
p = 3
q = 6
alpha = 0.1
omega = 1

# Number of samplepoints
N = 2**12
dt = 0.02 # Time Step
finalTime = 10

numberSteps = finalTime/dt

L = 500
# sample spacing (space step)
h = 0.0625

n = np.linspace(-80, 80, N)
t = np.linspace(0,finalTime,finalTime/dt)
x = n*h

#X,T = np.meshgrid(x,t)
k = 2*n*np.pi/L
u = np.zeros(shape=(N,int(finalTime/dt)))
u[:,0] = np.sin(50.0 * 2.0*np.pi*x) + 0.5*np.sin(80.0 * 2.0*np.pi*x)
#u[:,0] = (omega/(1+np.sqrt(1-alpha*omega)*np.cosh(2*p*np.sqrt(omega*x))))**(1/(2*p))+epsilon
#fig, ax = plt.subplots()
#plt.xlim(-N*h/2,N*h/2)
#plt.plot(x,np.absolute(u))

#plt.show()

for m in range(finalTime):

    u[:,m+1] = np.exp(dt*1j*(-2*(p+1)*np.absolute(u[:,m])*np.absolute(u[:,m])**p +alpha*(q+1)*np.absolute(u[:,m])*np.absolute(u[:,m])**q))*u[:,m] # Solve non-linear part of NLSE
    c = np.fft.fftshift(np.fft.fft(u[:,m+1]))
    c = np.exp(dt*1j*(-k*k))*c
    u[:,m+1] = np.fft.ifft(np.fft.fftshift(c))

#uAbs = np.absolute(u[:,0])
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


for timeCounter in range(int(finalTime/dt)):
    
    ax.plot(x, abs(u[:,timeCounter]), zs=timeCounter, zdir='y', alpha=0.8)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# On the y axis let's only label the discrete values that we have data for.
#ax.set_yticks(yticks)

plt.show()

#fig, ax = plt.subplots()
#plt.xlim(-N*h/2,N*h/2)
#plt.plot(x,uAbs)

#lt.show()
# %%
