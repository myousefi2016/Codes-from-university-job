# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 16:45:43 2016

@author: dkalita
"""
from __future__ import division

import sys
import numpy as np
from numpy import sqrt, sin, cos, pi
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import matplotlib.image as mpimg
import pylab as pl
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import AxesGrid

import collections
import pickle
import bisect
from scipy import integrate
from scipy.integrate import quad, dblquad
from scipy.optimize import curve_fit
from scipy import special as sp
import scipy.integrate as integrate
import scipy.special as special
import numpy.random
import math
import copy
import datetime

from multiprocessing import Pool

from operator import attrgetter


# function [ stressArray_in ] = DislocStresses( L,nIm)
# Disloc_Rotated function evaluates
# the stress fields on discrete points created by a
# formation of  wall of edge dislocation located in the
# center of a grid with size (2L-1)^2. nIm is the number of images wall in
# (positive and negative) x


class Dislocation:
    def __init__(self, x, y, s, i, tau=0.):
        self.x = x
        self.y = y
        self.s = s
        self.i = i
        self.tau = tau
        self.tau_avg = tau
        self.flag = True  # whether the dislocation exists

    def __lt__(self, other):
        return np.abs(self.tau) < np.abs(self.tau)

    def update_stress(self, stress, L=200):
        self.tau = stress[self.x, self.y, 2]

        dx = np.sign(self.tau) * np.sign(self.s)
        x2 = np.mod(self.x + dx, L)

        self.tau_avg = .5 * (stress[x2, self.y, 2] + self.tau)

    def __str__(self):
        return "i:{}\tx:{}\ty:{}\ts:{}\ttau:{}\ttau_avg:{}\tflag:{}".format(self.i, self.x, self.y, self.s, self.tau,
                                                                            self.tau_avg, self.flag)


class DislocationsDynamics(Dislocation):
    def __init__(self, n, L, stress_core, m, nsteps):
        self.n = n
        self.L = L
        self.m = m
        self.stress_core = stress_core
        self.stress = np.zeros((L, L, 3))

        self.abs_tau_xy = np.zeros(n)
        self.tau_xy = np.zeros(n)
        self.nx = np.zeros(n, dtype=int)
        self.ny = np.zeros(n, dtype=int)
        self.sign = np.zeros(n)
        self.disloc_position = []
        self.disloc_tau_xy_array = []
        self.disloc_w_max_tau_xy = []
        self.displaced_disloc_array = ()

        self.dislocs = []
        self.unique_dislocs = []
        self.ee_density = np.zeros((L, L))

        self.total_energy = 0
        self.total_energy = np.empty(nsteps)
        self.lx = 10.e-6
        self.dx = self.lx / self.L
        self.min_total_energy = []
        self.annihilated_dislocs_list = []

        self.dis_history = []

        self.nsteps = nsteps

        self.generate_random_disloc()
        self.update()

    def generate_random_disloc(self):
        # generates n number of dislocs
        n_positive = (self.m * self.n + self.n) * 0.5
        seed = 1000
        numpy.random.seed(seed=999)
        for i in range(self.n):

            x = int(np.random.random() * self.L)
            y = int(np.random.random() * self.L)
            if i < n_positive:
                s = 1
            else:
                s = -1

            self.dislocs.append(Dislocation(x, y, s, i))

        npd = 0  # no. of positive dislocations
        nnd = 0  # no. of negative dislocations
        for d in self.dislocs:
            if (d.s == 1):
                npd += 1
            elif (d.s == -1):
                nnd += 1
            else:
                assert (False, 'you should not get here')
        print('# total positive and negative dislocs: ', len(self.dislocs), n_positive, self.m, npd, nnd)

    def generate_random_disloc_alternative(self):
        # generates n number of dislocs
        n_positive = (self.m * self.n + self.n) * 0.5
        seed = 1000
        numpy.random.seed(seed=999)
        for i in range(self.n):
            x = int(np.random.random() * self.L)
            y = int(np.random.random() * self.L)

            s = np.random.choice([1, -1])
            self.dislocs.append(Dislocation(x, y, s, i))

        npd = 0  # no. of positive dislocations
        nnd = 0  # no. of negative dislocations
        for d in self.dislocs:
            if (d.s == 1):
                npd += 1
            elif (d.s == -1):
                nnd += 1
            else:
                assert (False, 'you should not get here')
        print('# total positive and negative dislocs: ', len(self.dislocs), n_positive, self.m, npd, nnd)

    def annihilate_dislocations_at_generation(self):

        for d1 in self.dislocs:
            for d2 in self.dislocs:
                if (d1.x == d2.x and d1.y == d2.y and d1.s != d2.s and d1.flag and d2.flag):
                    print('annihilated dislocations:', d1, d2)

                    d1.flag = False
                    d2.flag = False

        for d in self.dislocs:
            if d.flag == False:
                self.annihilated_dislocs_list.append(d)

        self.dislocs[:] = [d for d in self.dislocs if d.flag]

        for i, d in enumerate(self.dislocs):
            d.i = i

    def annihilate_dislocations_during_relaxation(self):

        for d1 in self.dislocs:
            for d2 in self.dislocs:
                if (d1.x == d2.x and d1.y == d2.y and d1.s != d2.s and d1.flag and d2.flag):
                    print('annihilated dislocations:', d1, d2)

                    d1.flag = False
                    d2.flag = False

        annihilated_dislocs_list = []
        for d in self.dislocs:
            if d.flag == False:
                annihilated_dislocs_list.append(d)

        for d in self.annihilated_dislocs_list:
            self.stress -= self.stress_core[(self.L - d.x - 1):(2 * self.L - d.x - 1),
                           (self.L - d.y - 1):(2 * self.L - d.y - 1), :] * int(d.s)

        self.dislocs[:] = [d for d in self.dislocs if d.flag]
        for i, d in enumerate(self.dislocs):
            d.i = i

    def annihilate_dislocs_old_method(self):

        for d in self.dislocs:

            if (d.x == self.dislocs[self.disloc_w_max_tau_xy.i].x and d.y == self.dislocs[
                self.disloc_w_max_tau_xy.i].y and d.s != self.dislocs[self.disloc_w_max_tau_xy.i].s and d.flag and
                    self.dislocs[self.disloc_w_max_tau_xy.i].flag):
                d.flag = False
                self.dislocs[self.disloc_w_max_tau_xy.i].flag = False

        self.dislocs[:] = [d for d in self.dislocs if d.flag]  # and  self.dislocs[0].flag]

    def compute_stress_field(self):

        for d in self.dislocs:
            self.stress += self.stress_core[(self.L - d.x - 1):(2 * self.L - d.x - 1),
                           (self.L - d.y - 1):(2 * self.L - d.y - 1), :] * int(d.s)

        self.update_dislocations_stress()

    def update_dislocations_stress(self):
        for d in self.dislocs:
            d.update_stress(self.stress, self.L)

    def sort_dislocations(self):
        self.dislocs.sort(key=lambda x: abs(x.tau_avg), reverse=True)

    def weighted_choice(self):
        total = sum(abs(d.tau) for d in self.dislocs)
        r = np.random.uniform(0, total)
        upto = 0
        i = 0  # returns the index of a chosen dislocation in array
        for d in self.dislocs:
            w = abs(d.tau)
            if upto + w >= r:
                return d
            upto += w

        assert False, "Shouldn't get here"

    def maximum_choice(self):
        return (max(self.dislocs, key=lambda d: np.abs(d.tau_avg)))

    def disloc_moving_criteria(self):
        # we chose from the sorted list a dislocation to move

        # self.disloc_w_max_tau_xy = copy.copy(self.maximum_choice())
        self.disloc_w_max_tau_xy = copy.copy(self.weighted_choice())

    def max_tau_xy(self):
        # we chose from the sorted list a dislocation to move

        d = np.random.randint(0, 10)
        self.disloc_w_max_tau_xy = copy.copy(self.maximum_choice())

        print('disloc with max tau is :', self.disloc_w_max_tau_xy)
        self.disloc_w_max_tau_xy = copy.copy(self.dislocs[d])
        self.disloc_w_max_tau_xy.i = d

    def move_disloc(self):
        dx = np.sign(self.disloc_w_max_tau_xy.tau) * np.sign(self.disloc_w_max_tau_xy.s)
        # we print to keep a track of the dislocation being moved
        # print ("{}:  move {} -> {}".format(self.dislocs[self.disloc_w_max_tau_xy.i].i,self.dislocs[self.disloc_w_max_tau_xy.i].x ,dx))
        self.dislocs[self.disloc_w_max_tau_xy.i].x = np.mod(self.dislocs[self.disloc_w_max_tau_xy.i].x + dx, self.L)

    def move_disloc_alternate_func(self):
        dx = np.sign(self.disloc_w_max_tau_xy.tau) * np.sign(self.disloc_w_max_tau_xy.s)

        print ("{}:  move {} -> {}".format(self.dislocs[self.disloc_w_max_tau_xy.i].i,
                                           self.dislocs[self.disloc_w_max_tau_xy.i].x, dx))

        self.dislocs[0].x = np.mod(self.dislocs[0].x + dx, self.L)
        self.annihilate_dislocs()

    def update(self):
        self.annihilate_dislocations_at_generation()
        self.compute_stress_field()

    def update_stress(self):
        nx_old = self.disloc_w_max_tau_xy.x
        ny_old = self.disloc_w_max_tau_xy.y

        i = self.disloc_w_max_tau_xy.i

        stress_old = self.stress_core[(self.L - nx_old - 1):(2 * self.L - nx_old - 1),
                     (self.L - ny_old - 1):(2 * self.L - ny_old - 1), :] * int(self.disloc_w_max_tau_xy.s)

        stress_new = self.stress_core[(self.L - self.dislocs[i].x - 1):(2 * self.L - self.dislocs[i].x - 1),
                     (self.L - self.dislocs[i].y - 1):(2 * self.L - self.dislocs[i].y - 1), :] * int(self.dislocs[i].s)

        self.stress += stress_new

        self.stress -= stress_old

        self.update_dislocations_stress()

    def compute_elastic_energy_density(self):
        # ee_density = elastic energy density
        # G= shear modulus for aluminium in Giga pascals ; nu=poisson's ratio
        G = 27.e9
        nu = 0.33

        self.ee_density = (
                              ((self.stress[:, :, 0] + self.stress[:, :, 1]) * (
                                  self.stress[:, :, 0] + self.stress[:, :, 1])) * (
                                  1 - nu) + 2 * (
                                  (self.stress[:, :, 2]) ** 2 - self.stress[:, :, 0] * self.stress[:, :, 1])) / 4. / G

        return self.ee_density

    def compute_total_energy(self):

        energy = np.sum(self.ee_density)  # * self.dx * self.dx

        return energy

    def run(self, nsteps=1):

        self.total_energy = np.empty(self.nsteps)
        self.total_energy.fill(numpy.nan)

        for i in np.arange(self.nsteps):
            # print("##################step {}#################".format(i))

            self.disloc_moving_criteria()
            self.move_disloc()
            self.update_stress()
            # self.annihilate_dislocations_at_generation()
            self.annihilate_dislocations_during_relaxation()

            self.compute_elastic_energy_density()

            self.total_energy[i] = self.compute_total_energy()
            # self.dis_history.append(copy.deepcopy(self.dislocs))


            if (i > 150 and False):
                E0 = self.total_energy[0]
                E1 = np.mean(self.total_energy[i - 5:i])
                dE = np.abs(self.total_energy[i - 1] - self.total_energy[i]) / np.abs(E0 - E1)

                if (dE < 5.e-5):
                    print("Energy reached minimum in {} steps".format(i))
                    break

        print('total energy evolution over one simulation of nsteps', self.total_energy)

        index_min_energy = np.argmin(self.total_energy)
        print("Minimum total energy for m= ", self.m, "is :", self.total_energy[index_min_energy])
        self.min_total_energy.append((self.total_energy[index_min_energy]))

    def print_dislocs(self):
        for d in self.dislocs:
            print(d)


def disloc_stresses(L, nIm):
    stress_array_xy = np.zeros((2 * L - 1, 2 * L - 1))
    stress_array_xx = np.zeros((2 * L - 1, 2 * L - 1))
    stress_array_yy = np.zeros((2 * L - 1, 2 * L - 1))

    b = 1

    # Evaluation
    eps = 1.e-6
    # Sum over all images
    domain = 2 * nIm + 1
    for n in range(domain):

        # Go through all cells
        for i in range(0, 2 * L - 1):
            for j in range(0, 2 * L - 1):
                # FIRST PART
                # wall towards y direction
                # images to x direction
                dx = i - L + L * (nIm + 1 - n) + 1
                dy = j - L + 1

                denom1 = (np.cosh(2 * np.pi * dx / L) - np.cos(2 * np.pi * dy / L)) ** 2
                if (denom1 > eps):
                    addStress_xy = b * 2 * np.pi * dx / L * (
                        np.cosh(2 * np.pi * dx / L) * np.cos(2 * np.pi * dy / L) - 1) / denom1
                    addStress_xx = b * np.sin(2 * np.pi * dy / L) * (
                        np.cosh(2 * np.pi * dx / L) - np.cos(2 * np.pi * dy / L) + 2 * np.pi * dx / L * np.sinh(
                            2 * np.pi * dx / L)) / denom1
                    addStress_yy = b * np.sin(2 * np.pi * dy / L) * (
                        np.cosh(2 * np.pi * dx / L) - np.cos(2 * np.pi * dy / L) - 2 * np.pi * dx / L * np.sinh(
                            2 * np.pi * dx / L)) / denom1
                else:
                    addStress_xy = 0.
                    addStress_xx = 0.
                    addStress_yy = 0.

                stress_array_xy[i, j] = stress_array_xy[i, j] + np.sum(addStress_xy)
                stress_array_xx[i, j] = stress_array_xx[i, j] + np.sum(addStress_xx)
                stress_array_yy[i, j] = stress_array_yy[i, j] + np.sum(addStress_yy)

    n_components = 3
    stressArray_in = np.zeros((2 * L - 1, 2 * L - 1, n_components))

    stressArray_in[:, :, 0] = stress_array_xx
    stressArray_in[:, :, 1] = stress_array_yy
    stressArray_in[:, :, 2] = stress_array_xy

    return stressArray_in


def demo_grid_with_each_cbar(fig):
    """
    A grid of 2x2 images. Each image has its own colorbar.
    """

    grid = AxesGrid(fig, 143,  # similar to subplot(143)
                    nrows_ncols=(3, 3),
                    axes_pad=0.1,
                    label_mode="1",
                    share_all=True,
                    cbar_location="top",
                    cbar_mode="each",
                    cbar_size="7%",
                    cbar_pad="2%",
                    )
    Z, extent = get_demo_image()
    for i in range(4):
        im = grid[i].imshow(Z, extent=extent, interpolation="nearest")
        grid.cbar_axes[i].colorbar(im)

    for cax in grid.cbar_axes:
        cax.toggle_label(False)

    # This affects all axes because we set share_all = True.
    grid.axes_llc.set_xticks([-2, 0, 2])
    grid.axes_llc.set_yticks([-2, 0, 2])


def main(argv):
    plt.close('all')
    now = datetime.datetime.now()
    nsim = 300
    # calling disloc_stresses function
    n = 100
    L = 200
    nIm = 5
    nsteps = 6000
    nM = 21

    def run_DDD(n, L, stress_core, m, nsteps):
        dd = DislocationsDynamics(n, L, stress_core, m, nsteps)
        dd.run()

        alive_dilocs = len(dd.dislocs)
        print('No of dislocs alive: ', alive_dilocs)
        return dd.total_energy, alive_dilocs

    def get_dislocs(n, L, stress_core, m, nsteps):

        dd = DislocationsDynamics(n, L, stress_core, m, nsteps)
        return dd.dislocs

    # load if disloc core exist- if not generate it and save it !

    fname = "disloc_core_L_{num1:06d}_nIm_{num2:03d}.npy".format(num1=L, num2=nIm)

    try:
        stress_core = np.load(fname)
        print("load core from file: ", fname)
    except:
        print("generate core and save to file: ", fname)
        stress_core = disloc_stresses(L, nIm)
        np.save(fname, stress_core)

    # calling DislocationsDynamics class

    mlist = np.linspace(0.0, 1.0, nM)

    print(mlist)

    alive_dislocs_data = np.empty(shape=(mlist.size, nsim))
    alive_dislocs_data.fill(numpy.nan)

    average_min_energy = np.empty(mlist.size)
    ensemble_energy = np.empty(shape=(nM, nsim, nsteps))
    ensemble_energy.fill(numpy.nan)
    for i, m in enumerate(mlist):

        min_energies = np.empty(nsim)
        for j in np.arange(nsim):
            print('######### m  ###### :', m, ' ######## sim no. ########  :', j)

            e, d = run_DDD(n, L, stress_core, m, nsteps)

            ensemble_energy[i][j] = e
            min_energies[j] = e[-1]
            alive_dislocs_data[i][j] = d

        average_min_energy[i] = min_energies.mean()

    data = {}
    day = now.day
    month = now.month
    year = now.year
    hour = now.hour
    minute = now.minute
    np.savez(
        'data_n_{:06d}_nsim_{:06d}_nsteps_{:06d}_date_{:02d}_{:02d}_{:04d}_{:02d}_{:02d}.npz'.format(n, nsim, nsteps,
                                                                                                     day, month, year,
                                                                                                     hour, minute),
        ensemble_energy=ensemble_energy, average_min_energy=average_min_energy, mlist=mlist,
        alive_dislocs_data=alive_dislocs_data)

    return 0


def plot_ensemble(filename='data.npz'):
    now = datetime.datetime.now()
    hour = now.hour
    minute = now.minute

    lx = 10.e-6
    L = 200
    dx = lx / L

    data = np.load(filename)
    print(data.keys())
    ensemble_energy = data['ensemble_energy']  # *dx*dx

    average_min_energy = data['average_min_energy']  # *dx*dx
    mlist = np.array(data['mlist'])
    alive_dislocs_data = data['alive_dislocs_data']
    # average_min_energy[:] = [energy * 10.e8 for energy in average_min_energy]

    no_of_dislocs_per_m = (np.sum(alive_dislocs_data, axis=1) / 500).tolist()  # /500 means divided by nsim for each m
    average_min_energy_per_dislocation = average_min_energy / no_of_dislocs_per_m
    print(mlist)
    if (False):
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)

        ax1.plot(mlist, average_min_energy, ls='-')
        p1 = np.poly1d(np.polyfit(mlist, average_min_energy, 1))
        print(p1(mlist))
        ax1.plot(mlist, p1(mlist), ls='--', lw=1.5)
        ax1.set_ylim([5.5e-5, 1.e-4])
        ax1.set_ylabel('Average minimum energy for every $M$ (J)')
        ax1.set_xlabel('Values of $M$')
        plt.grid(True)
        plt.gca().get_yaxis().get_major_formatter().set_powerlimits((0, 0))
        # fig1.savefig('/home/dkalita/Pictures/Average_min_energy_for_every_M.png', format='png', dpi=1500)


        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)

        ax2.plot(mlist, average_min_energy_per_dislocation, ls='-')
        p2 = np.poly1d(np.polyfit(mlist, average_min_energy_per_dislocation, 1))
        print(p2(mlist))
        ax2.plot(mlist, p2(mlist), ls='--', lw=1.5)
        ax2.set_ylim([1.e-7, 16.e-7])
        ax2.set_ylabel('Average minimum energy per dislocation for every $M$ (J)')
        ax2.set_xlabel('Values of $M$')
        plt.grid(True)
        plt.gca().get_yaxis().get_major_formatter().set_powerlimits((0, 0))
        # fig2.savefig('/home/dkalita/Pictures/Average_min_energy_per_disloc_for_every_M.png', format='png', dpi=1500)


        # fig1.savefig('Ensemble_energy_time_{:02d}_{:02d}.pdf'.format(hour,minute))#saves figure according to the time of running

        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)

        ax3.plot(mlist, no_of_dislocs_per_m, ls='-')
        p3 = np.poly1d(np.polyfit(mlist, no_of_dislocs_per_m, 2))
        ax3.plot(mlist, p3(mlist), ls='--', lw=1.5)
        ax3.set_ylabel('Average no. of dislocations alive per $M$')
        ax3.set_xlabel('Values of $M$')
        plt.grid(True)
        # fig3.savefig('/home/dkalita/Pictures/dislocs_alive_per_M.png', format='png', dpi=1500)

    data_shape = ensemble_energy.shape
    print(ensemble_energy.shape)

    color = cm.rainbow(np.linspace(0, 1, 21))
    if (True):
        fig4 = plt.figure(figsize=(5, 5))
        ax4 = fig4.add_subplot(111)
        box = ax4.get_position()

        if (True):
            for i, m in enumerate(mlist[::10]):
                avg_over_time = np.mean(ensemble_energy[i], axis=0)
                print(avg_over_time.shape)
                ax4.plot(avg_over_time, alpha=1, label='M={:0.2}'.format(m), color=color[i])
                ax4.set_ylabel('Average energies for different $M$ (J)')
                ax4.set_xlabel('Time (dislocation motion steps)')
                plt.gca().get_yaxis().get_major_formatter().set_powerlimits((0, 0))

                # ax4.set_position([box.x0, box.y0 - box.height * 0.1, box.width, box.height * 0.95])
                ax4.set_position([box.x0, box.y0, box.width * 0.8, box.height])
                ax4.legend(bbox_to_anchor=(1.025, 0.5), loc='center left', borderpad=0.5, fontsize='small')
                # fig4.savefig('/home/dkalita/Pictures/average_energy_evolution_evry_M.png', format='png', dpi=1500)

                if (True):
                    fig5 = plt.figure(figsize=(5, 5))
                    ax5 = fig5.add_subplot(111)
                    ax5.plot(avg_over_time, alpha=1, color='k', label='M={:0.2}'.format(m))
                    plt.gca().get_yaxis().get_major_formatter().set_powerlimits((0, 0))

                    for j in np.arange(data_shape[1]):
                        ax5.plot(ensemble_energy[i][j], alpha=.05, color='k')

                    ax5.set_ylabel('Total energy of the system (J)')
                    ax5.set_xlabel('Time (dislocation motion steps)')
                plt.legend()

    plt.show()


if __name__ == '__main__':
    now = datetime.datetime.now()

    day = now.day
    month = now.month
    year = now.year
    hour = now.hour
    minute = now.minute

    plt.close('all')
    main(sys.argv)  # to run the code and save the results
    # fname is to read from previuosly saved results according to the date and time of running the code
    fname = 'data_n_{:06d}_nsim_{:06d}_nsteps_{:06d}_date_{:02d}_{:02d}_{:04d}_{:02d}_{:02d}.npz'.format(100, 500, 6000,
                                                                                                         28, 12, 2016,
                                                                                                         18, 19)

    # to plot the saved results
    plot_ensemble(fname)
    sys.exit()
