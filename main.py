from lib.norcb import *
from skgeom.draw import draw
from numpy import abs
from copy import deepcopy
from numba import jit


def create_domain(minx, miny, maxx, maxy):
    p1 = sg.Point2(minx, miny)
    p2 = sg.Point2(maxx, miny)
    p3 = sg.Point2(maxx, maxy)
    p4 = sg.Point2(minx, maxy)
    corners = [p1,p2,p3,p4]
    chull = sg.convex_hull.graham_andrew(corners)
    pseg = zip(corners, rotate(corners, 1))
    segments = []
    for b, e in pseg:
        segments.append(sg.Segment2(b, e))

    return sg.Polygon(corners)


if __name__ == '__main__':
    minx, miny, maxx, maxy = 0, 0, 2., 2.
    domain = create_domain(minx, miny, maxx, maxy)

    PE = 64
    N  = 2**12

    #__points     = np.zeros((N, 2))
    #__points[:int(N/2), :]   = uniform_circle(int(N/2), cx=1,   cy=1,   R=0.5)
    __points = uniform_circle(N, cx=1.0, cy=1.0)

    print(__points)

    #bias         = np.zeros((int(N/2), 2))
    #bias[:, 0]   = np.random.rand(int(N/2)) - 0.5
    #bias[:, 1]   = np.random.rand(int(N/2)) - 0.5

    velocities   = expansion_velocity(__points, [1.0, 1.0]) #np.random.rand(N, 2) #expansion_velocity(__points)
    #velocities[:int(N/2), :] = np.random.rand(int(N/2), 2) + bias
    #velocities[int(N/2):, :] = np.random.rand(int(N/2), 2) + bias

    rcb_velocities   = deepcopy(velocities)
    norcb_velocities = deepcopy(velocities)

    average_vel  = np.mean(velocities, axis=0)

    rcb_points   = deepcopy(__points)
    norcb_points = deepcopy(__points)

    rcb   = partitioning_rcb(None, PE, rcb_points,    velocities, domain=domain, i=1, depth=0)
    norcb =     partitioning(plt, PE, norcb_points,  velocities, domain=domain, i=1, depth=0)
    p0bis = norcb[0][0]

    plt.show()
    dt = 2e-3
    NITER = 100
    rcb_imbalance   = []
    U_rcb = [0] * NITER
    U_norcb = [0] * NITER
    norcb_imbalance = []

    points = __points
    for i in range(NITER):
        W = eval(norcb, points) * 1e-6
        derose = (np.max(W) - np.mean(W))
        U_norcb[i] = (np.max(W) - np.mean(W)) if i == 0 else U_norcb[i-1] + (np.max(W) - np.mean(W))
        norcb_imbalance.append(compute_imbalance(W))

        W = eval(rcb, points) * 1e-6
        derose = (np.max(W) - np.mean(W))
        U_rcb[i] = (np.max(W) - np.mean(W)) if i == 0 else U_rcb[i - 1] + (np.max(W) - np.mean(W))
        rcb_imbalance.append(compute_imbalance(W))
        points, velocities = move(minx, miny, maxx, maxy, points, velocities, dt)

        fig, ax = plt.subplots(2, 1, figsize=(5, 5))
        show_partitioned_points(ax[0], norcb,   points)
        show_partitioned_points(ax[1],   rcb,   points)
        ax[0].set_xlim(minx, maxx)
        ax[1].set_xlim(minx, maxx)
        ax[0].set_ylim(miny, maxy)
        ax[1].set_ylim(miny, maxy)
        ax[0].set_aspect('equal', 'box')
        ax[1].set_aspect('equal', 'box')
        fig.tight_layout()
        plt.savefig(f"img-{i:03d}.png", dpi=300)
        plt.close()

    print(U_rcb, U_norcb)

    plt.figure()
    plt.title("Imbalance over time")
    plt.xlabel("Imbalance metric")
    plt.ylabel("Iteration")
    plt.plot(rcb_imbalance, label='rcb')
    plt.plot(norcb_imbalance, label='norcb')
    plt.legend()
    plt.tight_layout()
    plt.show()
    plt.close()

    plt.figure()
    plt.title("U over time")
    plt.xlabel("U")
    plt.ylabel("Iteration")
    plt.plot(U_rcb, label='rcb')
    plt.plot(U_norcb, label='norcb')
    plt.legend()
    plt.tight_layout()
    plt.show()
    plt.close()


#     intersections = []
#     for i, bis1 in enumerate(p0bis):
#         intersections = intersections + find_intersections_with_bisectors(bis1, p0bis[i+1:])
#     intersections = list(filter(lambda i: is_vertex(i, p0bis), intersections))
#     print(intersections)
#     plt.show()
#
#     eval_partition(norcb, norcb_points)
#     show_partitioned_points(ax[0][0], norcb, norcb_points)
#     show_partitioned_points(ax[0][1], rcb,     rcb_points)
#     dt = 0.01
#
#     for i in range(100):
#         rcb_points = rcb_points + dt * velocities
#         rcb_imbalance.append(compute_imbalance(eval_partition(rcb, rcb_points)))
#
#     for i in range(100):
#         norcb_points = norcb_points + dt * velocities
#         norcb_imbalance.append(compute_imbalance(eval_partition(norcb, norcb_points)))
#
#     ax[2][1].plot(rcb_imbalance, label='RCB')
#     ax[2][1].set_ylabel('imbalance')
#     ax[2][1].set_xlabel('iteration')
#     ax[2][1].legend()
#
#     ax[2][0].plot(norcb_imbalance, label='NoRCB')
#     ax[2][0].set_xlabel('iteration')
#     ax[2][0].set_ylabel('imbalance')
#     ax[2][0].legend()
#
#     ax[0][0].set_title(f"{PE} {N} {average_vel}")
#     ax[2][0].set_ylim(0, 2 + max( max(norcb_imbalance), max(rcb_imbalance) ))
#     ax[2][1].set_ylim(0, 2 + max( max(norcb_imbalance), max(rcb_imbalance) ))
#
#     show_partitioned_points(ax[1][1], rcb, rcb_points)
#     show_partitioned_points(ax[1][0], norcb, norcb_points)
#
#     plt.show()
#
