import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from copy import deepcopy

from lib.show import *
from lib.math import *
from lib.norcb import *




if __name__ == '__main__':

    # Generate random points
    PE = 16
    N  = 8192

    rcb_imbalance   = []
    norcb_imbalance = []

    __points     = uniform_circle(N, R=1)
    bias         = np.zeros_like(__points)

    bias[:, 0]   = np.random.rand(N) * 2.
    bias[:, 1]   = np.random.rand(N) - 0.5

    velocities   = expansion_velocity(__points)
    average_vel  = np.mean(velocities, axis=0)

    rcb_points   = deepcopy(__points)
    norcb_points = deepcopy(__points)
    f, ax        = plt.subplots(nrows=3, ncols=2)

    ax[0][0].set_xlim((-4, 6))
    ax[0][0].set_ylim((-4, 6))

    ax[0][1].set_xlim((-4, 6))
    ax[0][1].set_ylim((-4, 6))

    ax[1][0].set_xlim((-4, 6))
    ax[1][0].set_ylim((-4, 6))

    ax[1][1].set_xlim((-4, 6))
    ax[1][1].set_ylim((-4, 6))

    rcb   = partitioning_rcb(None, PE, rcb_points,    velocities, i=1, depth=0)
    norcb =     partitioning(None, PE, norcb_points,  velocities, i=1, depth=0)

    eval_partition(norcb, norcb_points)
    show_partitioned_points(ax[0][0], norcb, norcb_points)
    show_partitioned_points(ax[0][1], rcb,     rcb_points)
    dt = 0.01

    for i in range(100):
        rcb_points = rcb_points + dt * velocities
        rcb_imbalance.append(compute_imbalance(eval_partition(rcb, rcb_points)))

    for i in range(100):
        norcb_points = norcb_points + dt * velocities
        norcb_imbalance.append(compute_imbalance(eval_partition(norcb, norcb_points)))

    ax[2][1].plot(rcb_imbalance, label='RCB')
    ax[2][1].set_ylabel('imbalance')
    ax[2][1].set_xlabel('iteration')
    ax[2][1].legend()

    ax[2][0].plot(norcb_imbalance, label='NoRCB')
    ax[2][0].set_xlabel('iteration')
    ax[2][0].set_ylabel('imbalance')
    ax[2][0].legend()

    ax[0][0].set_title(f"{PE} {N} {average_vel}")
    ax[2][0].set_ylim(0, 2 + max( max(norcb_imbalance), max(rcb_imbalance) ))
    ax[2][1].set_ylim(0, 2 + max( max(norcb_imbalance), max(rcb_imbalance) ))

    show_partitioned_points(ax[1][1], rcb, rcb_points)
    show_partitioned_points(ax[1][0], norcb, norcb_points)

    plt.show()

