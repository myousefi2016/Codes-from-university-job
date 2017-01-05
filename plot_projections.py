# -*- coding: utf-8 -*-
"""
Provides a set of function for  projecting  axe3d plots in orthogonal parallel view 
or projecting 3d plot on xy,xz or yz  plane
"""
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d
from matplotlib.patches import *
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.text import TextPath
from matplotlib.transforms import Affine2D
import matplotlib.lines as mlines
from itertools import combinations, product


def orthogonal_proj(zfront, zback):
    a = (zfront + zback) / (zfront - zback)
    b = -2 * (zfront * zback) / (zfront - zback)

    return np.array([[1, 0, 0, 0],
                     [0, 1, 0, 0],
                     [0, 0, a, b],
                     [0, 0, -0.0001, zback]])


# Later in your plotting code ...

def project_parallel():
    proj3d.persp_transformation = orthogonal_proj


def project_xy(ax):
    proj3d.persp_transformation = orthogonal_proj
    ax.view_init(elev=-90, azim=-90)
    ax.set_zticks([])
    ax.w_zaxis.line.set_lw(0.)
    ax.set_zlabel('')


def project_xz(ax):
    proj3d.persp_transformation = orthogonal_proj
    ax.view_init(elev=0, azim=-90)
    ax.set_yticks([])
    ax.w_yaxis.line.set_lw(0.)
    ax.set_ylabel('')


def project_yz(ax):
    proj3d.persp_transformation = orthogonal_proj
    ax.view_init(elev=0, azim=0)
    ax.set_xticks([])
    ax.w_xaxis.line.set_lw(0.)
    ax.set_xlabel('')


theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
z = np.linspace(-2, 2, 100)
r = z ** 2 + 1
x = z + 1
y = z / 2


def text3d(ax, xyz, s, zdir="z", size=None, angle=0, usetex=False, **kwargs):
    x, y, z = xyz
    if zdir == "y":
        xy1, z1 = (x, z), y
    elif zdir == "y":
        xy1, z1 = (y, z), x
    else:
        xy1, z1 = (x, y), z

    text_path = TextPath((0, 0), s, size=size, usetex=usetex, **kwargs)
    trans = Affine2D().rotate(angle).translate(xy1[0], xy1[1])

    p1 = PathPatch(trans.transform_path(text_path))
    ax.add_patch(p1)
    art3d.pathpatch_2d_to_3d(p1, z=z1, zdir=zdir)


def circle_3d(ax, xyz, r, zdir=None, **kwargs):
    x, y, z = xyz
    p = Circle((x, y), r, **kwargs)
    ax.add_patch(p)
    art3d.pathpatch_2d_to_3d(p, z=z, zdir=zdir)


def rectangle_3d(ax, xyz, lx, ly, zdir=None, **kwargs):
    x, y, z = xyz
    p = Rectangle((x, y), lx, ly, **kwargs)
    ax.add_patch(p)
    art3d.pathpatch_2d_to_3d(p, z=z, zdir=zdir)


def plot_edge_dislocaiton_2d(ax, xyz, l=.1, zdir=None, AR=.2, color='.0', positive_edge=True, **kwargs):
    lx = l * (2 * positive_edge - 1)
    ly = lx * AR
    x, y, z = xyz
    print x, y, z

    xy = [[x - .5 * lx, y], [x + .5 * lx, y], [x + .5 * lx, y + ly], [x + .5 * ly, y + ly], [x + .5 * ly, y + lx],
          [x - .5 * ly, y + lx], [x - .5 * ly, y + ly], [x - .5 * lx, y + ly]]

    p = Polygon(xy, closed=True, color=color, **kwargs)
    ax.add_patch(p)
    art3d.pathpatch_2d_to_3d(p, z=z, zdir=zdir)


def plot_arrow(ax, xyz1, xyz2, zdir=None, **kwargs):
    x1, y1, z1 = xyz1
    x2, y2, z2 = xyz2

    p = FancyArrowPatch((x1, y1), (x2, y2), **kwargs)
    ax.add_patch(p)
    art3d.pathpatch_2d_to_3d(p, z=z1, zdir=zdir)


def plot_edge_dislocaiton_3d(ax, xyz1, xyz2, l=.1, zdir=None, color='.0', positive_edge=True, **kwargs):
    AR = .2
    x1, y1, z1 = xyz1
    x2, y2, z2 = xyz2

    plot_edge_dislocaiton_2d(ax, xyz1, l, zdir, AR, color, positive_edge, **kwargs)
    plot_edge_dislocaiton_2d(ax, xyz2, l, zdir, AR, color, positive_edge, **kwargs)

    d = l * .5 * AR * (2 * positive_edge - 1)
    if (zdir == 'x'):
        ax.plot([z1, z2], [x1, x2], [y1 + d, y2 + d], color)
    elif (zdir == 'y'):
        ax.plot([x1, x2], [z1, z2], [y1 + d, y2 + d], color)
    elif (zdir == 'z'):
        ax.plot([x1, x2], [y1 + d, y2 + d], [z1, z2], color)


def plot_block3d(ax, xyz1, xyz2, l=.1, **kwargs):
    x1, y1, z1 = xyz1
    x2, y2, z2 = xyz2

    P = []

    for x in ([x1, x2]):
        for y in ([y1, y2]):
            for z in ([z1, z2]):
                P.append([x, y, z])

    d = [0, 1]
    D = np.array(list(product(d, d, d)))
    C = np.array(list(combinations(np.arange(2 ** 3), 2)))
    for i, j in C:
        if (np.abs((D[i] - D[j])).sum() == 1):
            ax.plot3D([P[i][0], P[j][0]], [P[i][1], P[j][1]], [P[i][2], P[j][2]], **kwargs)


def test():
    fig = plt.figure(1, figsize=(4, 3), dpi=300)

    ax = fig.add_axes([0, 0, 1, 1], projection='3d')

    plot_edge_dislocaiton_2d(ax, x, y, z, 0)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    project_yz(ax)

    plt.show()
    print ax.get_proj()
