# pylint: disable=invalid-name
# pylint: disable=too-many-instance-attributes
# pylint: disable=too-few-public-methods
# pylint: disable=too-many-arguments
from math import pi
from math import sqrt

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from matplotlib import cm


from mpl_toolkits.mplot3d import Axes3D


class Param:
    """Container for holding all simulation parameters."""
    def __init__(self,
                 xmax: float,
                 res: int,
                 dt: float,
                 timesteps: int,
                 im_time: bool) -> None:

        self.xmax = xmax
        self.res = res
        self.dt = dt
        self.timesteps = timesteps
        self.im_time = im_time

        self.dx = 2 * xmax / res
        self.x = np.arange(-xmax + xmax / res, xmax, self.dx)
        self.dk = pi / xmax
        self.k = np.concatenate((np.arange(0, res / 2),
                                 np.arange(-res / 2, 0))) * self.dk


class Operators:
    """Container for holding operators and wavefunction coefficients."""
    def __init__(self, res: int) -> None:

        self.V = np.empty(res, dtype=complex)
        self.R = np.empty(res, dtype=complex)
        self.K = np.empty(res, dtype=complex)
        self.wfc = np.empty(res, dtype=complex)


def init(par: Param, voffset: float, wfcoffset: float) -> Operators:
    """Initialize the wavefunction coefficients and the potential."""
    opr = Operators(len(par.x))
    opr.V = 0.5 * (par.x - voffset) ** 2
    opr.wfc = np.exp(-((par.x - wfcoffset) ** 2) / 2, dtype=complex)

    if par.im_time:
        opr.K = np.exp(-0.5 * (par.k ** 2) * par.dt)
        opr.R = np.exp(-0.5 * opr.V * par.dt)
    else:
        opr.K = np.exp(-0.5 * (par.k ** 2) * par.dt * 1j)
        opr.R = np.exp(-0.5 * opr.V * par.dt * 1j)
    return opr


def split_op(par: Param, opr: Operators) -> None:

    u = np.zeros(shape=(par.timesteps,par.res))
    u[0,:] = opr.wfc
    for i in range(par.timesteps):

        # Half-step in real space
        opr.wfc *= opr.R

        # FFT to momentum space
        opr.wfc = np.fft.fft(opr.wfc)

        # Full step in momentum space
        opr.wfc *= opr.K

        # iFFT back
        opr.wfc = np.fft.ifft(opr.wfc)

        # Final half-step in real space
        opr.wfc *= opr.R

        u[i,:] = opr.wfc
        # Density for plotting and potential
        density = np.abs(opr.wfc) ** 2

        # Renormalizing for imaginary time
        if par.im_time:
            renorm_factor = sum(density) * par.dx
            opr.wfc /= sqrt(renorm_factor)

        # Outputting data to file. Plotting can also be done in a
        # similar way. This is set to output exactly 100 files, no
        # matter how many timesteps were specified.
        #if i % (par.timesteps // 100) == 0:
        #    filename = "output{}.dat".format(str(i).rjust(5, str(0)))
        #    with open(filename, "w") as outfile:
        #        # Outputting for gnuplot. Any plotter will do.
        #        for j in range(len(density)):
        #            template = "{}\t{}\t{}\n".format
        #            line = template(par.x[j], density[j].real, opr.V[j].real)
        #            outfile.write(line)
        print("Outputting step: ", i + 1)
    return u


def calculate_energy(par: Param, opr: Operators) -> float:
    """Calculate the energy <Psi|H|Psi>."""
    # Creating real, momentum, and conjugate wavefunctions.
    wfc_r = opr.wfc
    wfc_k = np.fft.fft(wfc_r)
    wfc_c = np.conj(wfc_r)

    # Finding the momentum and real-space energy terms
    energy_k = 0.5 * wfc_c * np.fft.ifft((par.k ** 2) * wfc_k)
    energy_r = wfc_c * opr.V * wfc_r

    # Integrating over all space
    energy_final = sum(energy_k + energy_r).real

    return energy_final * par.dx


def main() -> None:
    N = 2**12
    XMAX = 5.0
    DT = 0.05
    TIMESTEPS = 100
    par = Param(XMAX, N, DT, TIMESTEPS, False)

    # Starting wavefunction slightly offset so we can see it change
    opr = init(par, 0.0, -1.00)
    u = split_op(par, opr)
    #u = np.reshape(u,(int(TIMESTEPS/DT),N))
    fig = plt.figure(figsize=plt.figaspect(0.5))
    fig.suptitle('$i\\frac{\\partial \\psi}{\\partial t} = -\\frac{\\partial^2 \psi}{\\partial x^2}+0.5(x-x_0)^2$')

    ax = fig.add_subplot(1, 2, 1, projection='3d')
    X,T = np.meshgrid(par.x,np.arange(0,TIMESTEPS*DT,DT))


    #for timeCounter in range(TIMESTEPS):
        
        #ax.plot(par.x, abs(u[:,timeCounter]), zs=timeCounter, zdir='y', alpha=0.8,color=colors[])
    surf = ax.plot_surface(X, T, abs(u)**2, cmap=cm.jet, \
                       linewidth=1, antialiased=False,alpha=0.6)
    fig.colorbar(surf, shrink=0.5, aspect=5)

    ax.set_xlabel('X')
    ax.set_ylabel('T')
    ax.set_zlabel('$|u|^2$')


    ax = fig.add_subplot(1, 2, 2, projection='3d')

    surf = ax.plot_surface(X, T, abs(u)**2, cmap=cm.jet, \
                       linewidth=1, antialiased=False,alpha=0.6)
    fig.colorbar(surf, shrink=0.5, aspect=5)

    ax.set_xlabel('X')
    ax.set_ylabel('T')
    ax.set_zlabel('$|u|^2$')

    #ax = fig.add_subplot(1, 3, 3)
    #ax.plot(X, abs(u[-1,:])**2)



# On the y axis let's only label the discrete values that we have data for.
#ax.set_yticks(yticks)

    plt.show()

    #energy = calculate_energy(par, opr)
    #print("Energy is: ", energy)


if __name__ == "__main__":
    main()
