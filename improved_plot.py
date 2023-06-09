import ase
from ase.calculators.espresso import Espresso
from ase.visualize import view

PbTe_atoms = ase.io.read("PbTe.cif")

import cellconstructor as CC, cellconstructor.Phonons
import cellconstructor.ForceTensor
# Import the numerical libraries and those for plotting
import numpy as np
import matplotlib.pyplot as plt
import sys, os, time


def plot_gp(filename,color,LS):
        global ax
        d = np.loadtxt(filename)
        nband = int((len(d[0])-1))
        for i in range(nband):
                if i==0:
                        ax.plot(d[:,0],d[:,i+1]*cm_mev,color=color,ls=LS,label=f'{filename}')
                        
                else:
                        ax.plot(d[:,0],d[:,i+1]*cm_mev,ls=LS,color=color)
                      
        xticks = [d[0,0],d[200,0],d[-1,0]]
        # ax.set_xlim(xticks[0], xticks[1])
        ax.set_xticks(xticks)
        ax.set_xticklabels(["M", "$\\Gamma$", "L"])
        ax.legend()
        ax.set_ylabel("Phonons [meV]");ax.set_ylim(-1,18)
        # plt.show()

def plot_bubble(filename):
        global ax
        d = np.loadtxt(filename,skiprows=3)
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
        # plt.show()
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
        # axs[iplot].set_xlim(xticks[0], xticks[1])
        axs[iplot].set_xticks(xticks)
        axs[iplot].set_xticklabels(xlabels)
        axs[iplot].set_ylabel("Phonons [meV]");axs[iplot].set_ylim(-1,18)

start=time.time()
cm_mev=1/8.066
PATH=['GMKGALHAMLHK']
PATH=['GM']
nplot=len(PATH)

N_POINTS = 200
# Here we define the position of the special points
SPECIAL_POINTS = {'G': np.array([0., 0., 0.]),
 'M': np.array([ 0.5, -0.5,  0. ]),
 'K': np.array([ 0.66666667, -0.33333333,  0.        ]),
 'A': np.array([0. , 0. , 0.5]),
 'L': np.array([ 0.5, -0.5,  0.5]),
 'H': np.array([ 0.66666667, -0.33333333,  0.5       ])}


#plot_bubble('v2_v2+d3static_freq.dat')
#exit()


SSCHA_DYN = ['./target_666/harmonic_666_dynf2q_','target_121212/harmonic_121212_dynf2q_']
SSCHA_DYN = ['f2q_']
nq = [4]
COLOR = ['k', 'r']
sscha_dyn = [CC.Phonons.Phonons(SSCHA_DYN[i], nq[i]) for i in range(len(SSCHA_DYN))]
SSCHA_DYN = ['222_harmonic','121212']


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
end=time.time()
print(f'Time: {end-start} seconds')
plt.show()
print('finish')