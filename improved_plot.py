import ase
from ase.calculators.espresso import Espresso
from ase.visualize import view

PbTe_atoms = ase.io.read("PbTe.cif")

import cellconstructor as CC, cellconstructor.Phonons
import cellconstructor.ForceTensor
# Import the numerical libraries and those for plotting
import numpy as np
import matplotlib.pyplot as plt
import sys, os


def plot_bubble(filename):
        d = np.loadtxt(filename,skiprows=3)
        plt.figure(dpi = 150)
        ax = plt.gca()
        nband = int((len(d[0])-1)/2)
        for i in range(nband):
                if i==0:
                        ax.plot(d[:,0],d[:,i+1]*cm_mev,ls=':',color='k',label='v2')
                        ax.plot(d[:,0],d[:,i+1+nband]*cm_mev,color='r',label='v2+static_bubble')
                else:
                        ax.plot(d[:,0],d[:,i+1]*cm_mev,ls=':',color='k')
                        ax.plot(d[:,0],d[:,i+1+nband]*cm_mev,color='r')
        xticks = [d[0,0],d[-1,0]]
        ax.set_xlim(xticks[0], xticks[1])
        ax.set_xticks(xticks)
        ax.set_xticklabels(["$\\Gamma$", "L"])
        ax.legend()
        ax.set_ylabel("Phonons [meV]");ax.set_ylim(-1,18)
        plt.show()
# Let us define the PATH in the brilluin zone and the total number of points
def plot(iplot,PATH,sscha_dyn):
        global ax
        if nplot==1:
                axs=[ax]
        else:
                axs=ax
        qpath, data = CC.Methods.get_bandpath(sscha_dyn[0].structure.unit_cell, PATH, SPECIAL_POINTS,N_POINTS)
        xaxis, xticks, xlabels = data # Info to plot correclty the x axis
        print('xaxis',xaxis)
        print('xticks',xticks)
        print('xlabels',xlabels)
        sscha_dispersion = [CC.ForceTensor.get_phonons_in_qpath(sscha_dyn[i], qpath)*cm_mev for i in range(len(sscha_dyn))]
        nmodes = sscha_dyn[0].structure.N_atoms * 3

        if iplot>0:
                axs[iplot].yaxis.set_visible(False)
        for i in range(nmodes):
                        LS=['dashed','dashed','solid']
                        if i == 0:
                                for n in range(len(sscha_dyn)):

                                        axs[iplot].plot(xaxis, sscha_dispersion[n][:,i], color = COLOR[n],ls=LS[n], label = SSCHA_DYN[n])
                        else:
                                for n in range(len(sscha_dyn)):
                                        axs[iplot].plot(xaxis, sscha_dispersion[n][:,i],ls=LS[n], color = COLOR[n])

        for x in xticks:
                        # ax.axvline(x, 0, 1, color = "k", lw = 0.4)
                        axs[iplot].axhline(0, 0, 1, color = 'k', ls = ':', lw = 0.4)
        axs[iplot].set_xlim(xticks[0], xticks[1])
        axs[iplot].set_xticks(xticks)
        axs[iplot].set_xticklabels(xlabels)
        axs[iplot].set_ylabel("Phonons [meV]");axs[iplot].set_ylim(-1,18)


cm_mev=1/8.066
PATH=['GL']
nplot=len(PATH)

N_POINTS = 200
# Here we define the position of the special points
SPECIAL_POINTS = {"G": [0,0,0],
"X": [0.0,0, 0.5],
"L": [.5, .5, .5],
"W": [.5, .25, .75],
"K": [3/8., 3/8., 3/4],
"M": [0,0.5,0.5]}

plot_bubble('v2_v2+d3static_freq.dat')
exit()


SSCHA_DYN = ['./target_666/harmonic_666_dynf2q_','target_121212/harmonic_121212_dynf2q_']
SSCHA_DYN = ['./444/dyn_pop3_']
nq = [8,72]
COLOR = ['k', 'r']
sscha_dyn = [CC.Phonons.Phonons(SSCHA_DYN[i], nq[i]) for i in range(len(SSCHA_DYN))]
SSCHA_DYN = ['444','121212']


if nplot==1:
        plt.figure(dpi = 150)
        ax = plt.gca()
        plot(0,PATH[0],sscha_dyn)
        ax.legend()
else:
        fig, ax = plt.subplots(1, nplot,dpi = 150)
        for i in range(nplot):
                plot(i,PATH[i],sscha_dyn)
        ax[0].legend()

plt.tight_layout()
plt.savefig("dispersion.png")
plt.show()
print('finish')