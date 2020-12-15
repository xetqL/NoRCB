from lib.norcb import *
from skgeom.draw import draw
from numpy import abs


def unique(l, eq):
    r = []
    for i, x in enumerate(l):
        inside = any([eq(x, y) for y in l[i+1:]])
        if not inside:
            r.append(x)
    return r


def divide_segment(segment, pi):
    s1 = sg.Segment2(segment.source(), pi)
    s2 = sg.Segment2(pi, segment.target())

    return s1, s2


def point_to_numpy(point):
    return np.asarray([point.x(), point.y()], dtype=np.float64)


def segment_to_vector(segment):
    return point_to_numpy(segment.target()) - point_to_numpy(segment.source())


def segment_to_rvector(segment):
    return point_to_numpy(segment.source()) - point_to_numpy(segment.target())


def select_bisection(s1, s2, vec, point, p):
    _p = point_to_numpy(p)
    sgn= sign(side1(vec, _p - point_to_numpy(point)))
    if sgn > 0:
        s1.append(p)
    else:
        s2.append(p)


if __name__ == '__main__':
    plt.figure()
    p1 = sg.Point2( 0,  0)
    px = sg.Point2( 4,  0)
    p2 = sg.Point2(10,  5)
    p3 = sg.Point2(10, 10)
    p4 = sg.Point2( 0, 10)

    corners = [p1, px, p2, p3, p4]
    pseg = zip(corners, rotate(corners, -1))
    segments = []
    for b, e in pseg:
        segments.append(sg.Segment2(b, e))
    vp = np.asarray([ 1,  1])
    vm = np.asarray([-1, -1])
    pr1 = sg.Ray2(sg.Point2(1, 1), sg.Vector2(*vp))
    mr1 = sg.Ray2(sg.Point2(1, 1), sg.Vector2(*vm))

    # start with a set of segment and a ray
    # find intersections
    intersections = []
    for s in segments:
        iplus = (sg.intersection(mr1, s))
        iminus= (sg.intersection(pr1, s))
        i = iplus if iplus else iminus
        intersections.append(i)

    intersections = list(filter(lambda i: i, intersections))
    intersections = unique(intersections, lambda x, y: abs(float(x.x())) == abs(float(y.x())) and abs(float(x.y())) == abs(float(y.y())))

    sgbisect = sg.Segment2(*intersections)
    vec = segment_to_rvector(sg.Segment2(*intersections))

    s1 = [] + intersections
    s2 = [] + intersections
    for s in segments:
        select_bisection(s1, s2, vec, sgbisect.source(), s.source())
        select_bisection(s1, s2, vec, sgbisect.source(), s.target())

    chull = sg.convex_hull.graham_andrew(s1)
    phull = sg.Polygon(chull)
    print("1.", list(phull.edges))
    draw(phull)

    chull = sg.convex_hull.graham_andrew(s2)
    phull = sg.Polygon(chull)
    print("2.", list(phull.edges))
    draw(phull)

    plt.show()

#
#     # Generate random points
#     PE = 4
#     N  = 8192
#
#     rcb_imbalance   = []
#     norcb_imbalance = []
#
#     __points     = uniform_circle(N, R=1)
#     bias         = np.zeros_like(__points)
#
#     bias[:, 0]   = np.random.rand(N) * 2.
#     bias[:, 1]   = np.random.rand(N) - 0.5
#
#     velocities   = expansion_velocity(__points)
#     average_vel  = np.mean(velocities, axis=0)
#
#     rcb_points   = deepcopy(__points)
#     norcb_points = deepcopy(__points)
#     f, ax        = plt.subplots(nrows=3, ncols=2)
#
#     ax[0][0].set_xlim((-4, 6))
#     ax[0][0].set_ylim((-4, 6))
#
#     ax[0][1].set_xlim((-4, 6))
#     ax[0][1].set_ylim((-4, 6))
#
#     ax[1][0].set_xlim((-4, 6))
#     ax[1][0].set_ylim((-4, 6))
#
#     ax[1][1].set_xlim((-4, 6))
#     ax[1][1].set_ylim((-4, 6))
#
#     rcb   = partitioning_rcb(None, PE, rcb_points,    velocities, i=1, depth=0)
#     norcb =     partitioning(None, PE, norcb_points,  velocities, i=1, depth=0)
#     p0bis = norcb[0][0]
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
