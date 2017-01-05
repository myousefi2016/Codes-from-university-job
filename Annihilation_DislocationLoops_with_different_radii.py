import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pylab as pl
from mpl_toolkits.mplot3d import Axes3D
import collections
import pickle
import bisect
from scipy.optimize import curve_fit
from scipy import special as sp
import numpy.random
import math
import copy
from  plot_projections import *
from random import uniform
from scipy.optimize import fmin


class loop:
    """ dislocation loop class
    loop can expand freely or annihilate with another loop and become half 8 """

    def __init__(self, n, center=(0, 0, 0), r=1., v=1.):

        self.n = n
        self.center = center
        self.r0 = r
        self.v = v
        self.is_circle = True
        self.is_paired = False
        self.expand(0)
        self.qt = 2 * np.pi
        self.qt_c = 2 * np.pi

    def update(self):

        if self.is_paired:
            self.d = np.sqrt(
                (self.center[0] - self.paired_loop.center[0]) ** 2 + (self.center[1] - self.paired_loop.center[1]) ** 2)
            if (self.r + self.paired_loop.r > self.d):
                self.is_circle = False

        if (self.is_circle == True):
            self.p = 2 * self.r * np.pi
            self.area = np.pi * self.r ** 2


        elif (self.d < self.r):

            self.x = (self.r ** 2 - self.paired_loop.r ** 2 - self.d ** 2) / (self.d * 2)
            self.a = (self.d + self.x)
            self.h = np.sqrt(self.r ** 2 - self.a ** 2)

            self.phi = np.abs(np.arcsin(self.h / self.r))
            self.qt = 2 * (np.pi - self.phi)

            if (math.isnan(self.phi)):
                print self.n, self.r, self.paired_loop.r, self.d, self.h, self.a
                print self.r
                exit()
            self.p = 2 * self.r * (np.pi - self.phi)

            self.area = (np.pi - self.phi) * self.r ** 2 + self.h * self.a

        else:
            self.a = (self.d ** 2 + self.r ** 2 - self.paired_loop.r ** 2) / (self.d * 2)
            self.h = np.sqrt(self.r ** 2 - self.a ** 2)

            self.phi = np.abs(np.arcsin(self.h / self.r))
            self.qt = 2 * (np.pi - self.phi)

            if (math.isnan(self.phi)):
                print self.n, self.r, self.paired_loop.r, self.d, self.h, self.a
                print self.r
                exit()
            self.p = 2 * self.r * (np.pi - self.phi)

            self.area = (np.pi - self.phi) * self.r ** 2 + self.h * self.a

    def expand(self, t):
        self.r = self.r0 + self.v * t

    def pair(self, paired_loop):

        if (self.is_paired or paired_loop.is_paired):
            return False
        else:
            print "pairing: ", self.n, paired_loop.n

            self.paired_loop = paired_loop
            self.paired_loop.paired_loop = self
            self.is_paired = True
            self.paired_loop.is_paired = True
            self.update()
            self.paired_loop.update()
        return True


def check_intersection(new_loop, old_loops):
    if (len(old_loops) > 0):
        for i in range(len(old_loops)):
            d = np.sqrt(
                (old_loops[i].center[0] - new_loop.center[0]) ** 2 + (old_loops[i].center[1] - new_loop.center[1]) ** 2)
            if (d <= old_loops[i].r + new_loop.r):
                return True
    return False


def check_loop_inside_loop(new_loop, old_loops):
    if (len(old_loops) > 0):
        for i in range(len(old_loops)):
            d = np.sqrt(
                (old_loops[i].center[0] - new_loop.center[0]) ** 2 + (old_loops[i].center[1] - new_loop.center[1]) ** 2)
            if (old_loops[i].r >= d + new_loop.r):
                return True
    return False


def make_loops(nloops):
    loops = []
    total_r = 0.0
    i = 0
    while (i < nloops):
        r = r0
        center = np.random.random(3) * L

        loops.append(loop(i, center, r, v))
        total_r += r
        i += 1

    print ("avg r :", total_r / nloops)
    return loops


def make_loops_alternative(nloops):
    loops = []
    total_r = 0.0
    i = 0
    while (i < nloops):

        y_rand = uniform(0, 1)  # y_rand is a random value of the CDF
        r = -r0 * math.log(1 - y_rand)  # corresponding radius for a random value of the CDF
        center = np.random.random(3) * L

        while (check_intersection(loop(i, center, r, v), loops)):
            center = np.random.random(3) * L
            print "t"

        for i in range(nloops):
            loops.append(loop(i, center, r, v))
            total_r = total_r + r
            i += 1

    print ("avg r :", total_r / nloops)
    return loops


def pair_loops(loops, da):
    nloops = len(loops)

    for l1 in range(nloops):
        for l2 in range(l1 + 1, nloops):
            if (loops[l1].is_paired): break;
            dz = np.abs(loops[l1].center[2] - loops[l2].center[2])
            if (dz < da):
                if (loops[l1].pair(loops[l2])):
                    break;


def pair_loops_periodic(loops, da):
    nloops = len(loops)
    pair_candidates = []
    for l1 in range(nloops):
        for l2 in range(l1 + 1, nloops):

            dz = np.abs(loops[l1].center[2] - loops[l2].center[2])
            dx = np.abs(loops[l1].center[0] - loops[l2].center[0])
            dy = np.abs(loops[l1].center[1] - loops[l2].center[1])

            if (dz < da or L - dz < da):
                if (dx >= loops[l1].r - loops[l2].r and dy >= loops[l1].r - loops[l2].r):  # loops[l1].r>loops[l2].r and
                    # loop can be paired
                    # we check the minimum distance of l1 with l2 and its periodic images
                    d_array = []
                    for x0 in [-L, 0, L]:
                        for y0 in [-L, 0, L]:
                            d = np.sqrt((loops[l1].center[0] - loops[l2].center[0] - x0) ** 2 + (
                                loops[l1].center[1] - loops[l2].center[1] - y0) ** 2)
                            d_array.append((d, l1, l2, x0, y0))
                    pair_candidates.append(min(d_array))
    pair_candidates.sort()

    for p in pair_candidates:

        if (loops[p[1]].is_paired or loops[p[2]].is_paired):
            pass

        else:

            loops[p[2]].center = (loops[p[2]].center[0] + p[3], loops[p[2]].center[1] + p[4], loops[p[2]].center[2])
            loops[p[1]].pair(loops[p[2]])


def run():
    for i, t in enumerate(timesteps):
        for n, loop in enumerate(loops):
            loop.expand(t)
        for n, loop in enumerate(loops):
            loop.update()
            rhot[i] += loop.p
            qt[i] += loop.p / loop.r
            qt_c[i] += loop.qt_c

            gamma[i] += loop.area * bx
            ncircles[i] += loop.is_circle


class CDD:
    def __init__(self, n=0):
        self.rhot = np.zeros(n)
        self.qt = np.zeros(n)
        self.gamma = np.zeros(n)

    pass


do_DDD = True

nloops = 5  # n loops
nsimulations = 100  # number of DDD simulations

unit = 1.e-6  # meter
L = 3 * unit  # length
V = L ** 3  # Volume
r0 = .05 * unit  # radii
v = .001 * unit  # velocity
da = .02 * unit  # annihilation distance
bx = .256 * .001 * unit  # burgers vector
total_time = L / v * 0.03 * 5
nsteps = 200  # discritization time
timesteps = np.linspace(0, total_time, nsteps)
dt = total_time / nsteps

ax = []
hspace = .5
fig, ax = plt.subplots(3, 1, figsize=(8, 16), dpi=80)
fig.tight_layout()
pl.subplots_adjust(hspace=hspace)

alpha = 0.025

filename = 'nloops_{nloops:04d}_nsim_{nsimulations:04d}_r0_{r:.3e}_da_{da:.3e}.npz'.format(nloops=nloops,
                                                                                           nsimulations=nsimulations,
                                                                                           r=r0, da=da)

if do_DDD:

    rhot_avg = np.zeros(nsteps)
    qt_avg = np.zeros(nsteps)
    qt_c_avg = np.zeros(nsteps)

    gamma_avg = np.zeros(nsteps)
    ncircles_avg = np.zeros(nsteps)

    for i in range(nsimulations):
        print 'inside'
        rhot = np.zeros(nsteps)
        qt = np.zeros(nsteps)
        qt_c = np.zeros(nsteps)

        gamma = np.zeros(nsteps)
        ncircles = np.zeros(nsteps)

        loops = make_loops(nloops)  # r0 new argument

        pair_loops_periodic(loops, da)

        print"#######{0:03d}########".format(i)

        run()

        rhot = rhot / V
        qt = qt / V
        qt_c = qt_c / V

        gamma = gamma / V

        rhot_avg += rhot
        qt_avg += qt
        qt_c_avg += qt_c
        gamma_avg += gamma
        ncircles_avg += ncircles

        ncircles = np.zeros(nsteps)

        ax[0].plot(timesteps, rhot, color='r', ls='-', lw=1, alpha=alpha)
        ax[1].plot(timesteps, qt, color='r', ls='-', lw=1, alpha=alpha)
        ax[2].plot(timesteps, gamma, color='r', ls='-', lw=1, alpha=alpha)

    rhot_avg = rhot_avg / nsimulations
    qt_avg = qt_avg / nsimulations
    gamma_avg = gamma_avg / nsimulations
    ncircles_avg /= ncircles_avg / nsimulations

    print "average rhot  of DDD = ", rhot_avg

    simdata_nloops_nsimulations_r0_da_v_L_bx = np.array([nloops, nsimulations, r0, da, v, L, bx])
    np.savez_compressed(filename, rhot=rhot_avg, qt=qt_avg, gamma=gamma_avg, timesteps=timesteps,
                        simdata_nloops_nsimulations_r0_da_v_L_bx=simdata_nloops_nsimulations_r0_da_v_L_bx)

else:

    cdd = CDD(nsteps)

cdd_wo = CDD(nsteps)
cdd_wo.rhot = nloops * 2 * np.pi * (r0 + timesteps * v) / V
cdd_wo.qt = nloops * 2 * np.pi / V * (timesteps * 0 + 1)
cdd_wo.gamma = nloops * np.pi * (r0 + timesteps * v) ** 2 / V * bx

r = r0

C_qt = (4. + np.pi) / (2 * np.pi)
C_rhot = 1.

for i, t in enumerate(timesteps):
    r += v * dt


def qt_gamma(qt0, gamma, c):
    c2 = np.log(qt0) + (c * 2 * da / bx) * gamma[0]
    return (np.exp(-(c * 2 * da / bx) * gamma + c2))


ax[0].plot(timesteps, rhot_avg, color='.0', ls='--', lw=2, label=r"DDD: $\rho^t$")
ax[1].plot(timesteps, qt_avg, color='.0', ls='--', lw=2, label=r"DDD: $q^t$")
ax[2].plot(timesteps, gamma_avg, color='.0', ls='--', lw=2, label=r"DDD: $\gamma$")

ax[0].legend(loc="upper left")
ax[1].legend(loc="upper left")
ax[2].legend(loc="upper right")

# CDD simulation

rhot_cdd = np.zeros(nsteps)
qt_cdd = np.zeros(nsteps)
gamma_cdd = np.zeros(nsteps)

# intial values for CDD
rhot_cdd[0] = rhot_avg[0]
qt_cdd[0] = qt_avg[0]
gamma_cdd[0] = gamma_avg[0]

f = 2
time = len(timesteps) - 1

for i in range(time):
    rhot_cdd[i + 1] = rhot_cdd[i] + dt * (qt_cdd[i] * v - f * v * da * rhot_cdd[i] * rhot_cdd[i])
    qt_cdd[i + 1] = qt_cdd[i] - dt * (f * v * da * rhot_cdd[i] * qt_cdd[i])
    gamma_cdd[i + 1] = gamma_cdd[i] + dt * rhot_cdd[i] * v * bx

print "rhot of CDD = ", rhot_cdd

ax[0].plot(timesteps, rhot_cdd, color='g', ls='-', lw=1, label=r"CDD: $\rho^t$")
ax[1].plot(timesteps, qt_cdd, color='g', ls='-', lw=1, label=r"CDD: $q^t$")
ax[2].plot(timesteps, gamma_cdd, color='g', ls='-', lw=1, label=r"CDD: $\gamma$")

ax[0].legend(loc="upper left")
ax[1].legend(loc="lower left")
ax[2].legend(loc="upper left")
plt.rcParams['font.size'] = 18

# 3D DDD plot
unit = 10
da = .0001 * unit
fig3d = plt.figure(3, figsize=(15, 15), dpi=300)
ax3d = fig3d.add_subplot(111, projection='3d')

project_parallel()
plot_block3d(ax3d, [0, 0, 0], [i * unit for i in [1, 1, 1]], color='.0')
r0 = .05 * unit
L = 1 * unit
loops = make_loops(20)
pair_loops_periodic(loops, .1 * unit)

for loop in loops: loop.expand(0)
for loop in loops: loop.update()
for loop in loops: loop.update()
ncircles = 0
npaired = 0

for l in loops:

    if l.is_paired:
        npaired = npaired + 1
    if l.is_circle:
        ncircles = ncircles + 1

print ncircles, npaired
size = .8
lw = .001
fontsize = 6
for loop in loops:

    print
    if (not loop.is_paired):
        circle_3d(ax3d, loop.center, loop.r, 'z', ec='b', lw=.5, alpha=1, fill=False)
        ax3d.text(loop.center[0], loop.center[1], loop.center[2], "{}->{}".format(loop.n, '-'), fontsize=fontsize,
                  color='k')
        text3d(ax3d, loop.center, "{}->{}".format(loop.n, '-'), 'z', size=size, lw=lw)
    elif (not loop.is_circle):
        text3d(ax3d, loop.center, "{}->{}".format(loop.n, (loop.paired_loop.n)), 'z', size=size, lw=lw)
        ax3d.text(loop.center[0], loop.center[1], loop.center[2], "{}->{}".format(loop.n, loop.paired_loop.n),
                  fontsize=fontsize, color='k')

        circle_3d(ax3d, loop.center, loop.r, 'z', ec='r', lw=.5, alpha=1, fill=False)
    else:

        circle_3d(ax3d, loop.center, loop.r, 'z', ec='g', lw=.5, alpha=1, fill=False)
        text3d(ax3d, loop.center, "{}->{}".format(loop.n, (loop.paired_loop.n)), 'z', size=size, lw=lw)
        ax3d.text(loop.center[0], loop.center[1], loop.center[2], "{}->{}".format(loop.n, loop.paired_loop.n),
                  fontsize=fontsize, color='k')

        pass

print ncircles, npaired
ax3d.set_aspect('equal')

ax3d.axis('off')

plt.tight_layout()

fig.savefig("annihilation_density_max.png", dpi=300)

plt.show()

print timesteps[-1], total_time
