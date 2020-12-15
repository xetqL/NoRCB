from .utils import *
from .math import *
from .show import *
import skgeom as sg


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

VERY_SMALL=1e-12
def split_domain(segments, p, v):
    vp =  v
    vm = -v
    pr1 = sg.Ray2(sg.Point2(*p), sg.Vector2(*vp))
    mr1 = sg.Ray2(sg.Point2(*p), sg.Vector2(*vm))

    # start with a set of segment and a ray
    # find intersections
    intersections = []
    for s in segments:
        iplus = (sg.intersection(mr1, s))
        iminus = (sg.intersection(pr1, s))
        i = iplus if iplus else iminus
        intersections.append(i)

    intersections = list(filter(lambda i: i, intersections))
    intersections = unique(intersections,
                           lambda x, y: abs(float(x.x())-float(y.x())) < VERY_SMALL and abs(float(x.y()) -
                               float(y.y())) < VERY_SMALL)

    sgbisect = sg.Segment2(*intersections)

    vec = segment_to_rvector(sg.Segment2(*intersections))

    s1 = [] + intersections
    s2 = [] + intersections
    for s in segments:
        select_bisection(s1, s2, vec, sgbisect.source(), s.source())
        select_bisection(s1, s2, vec, sgbisect.source(), s.target())

    chull = sg.convex_hull.graham_andrew(s1)
    phull1 = sg.Polygon(chull)
    #print("1.", list(phull.edges))

    chull = sg.convex_hull.graham_andrew(s2)
    phull2 = sg.Polygon(chull)
    #print("2.", list(phull.edges))
    #draw(phull)
    return phull1, phull2


def belongs_to_partition(partition, p):
    bis, _, _ = partition
    return does_satisfy_all_bisectors(bis, p)


def point2(p):
    return sg.Point2(p[0], p[1])


def vel2(v):
    return point2(v)


def partitioning(ax, P, points, velocities, domain, bisectors=[], i=0, depth=0 ):
    N = len(points)
    width = np.amax(points, axis=0) - np.amin(points, axis=0)

    origin = [0., 1.]
    avg_vel = get_average_velocity(velocities, axis=np.argmin(width))

    if N <= 1 or P <= 1:
        sg.draw.draw(domain)
        if ax:
            show_points(ax, points, s=0.5)
        return [(bisectors, i - 2**depth)]

    theta = get_angle(avg_vel, origin)

    # rotate
    R  = get_rotation_matrix( theta)
    mR = get_rotation_matrix(-theta)

    avg_vel = (np.dot(mR, np.abs(np.dot(R, avg_vel))))
    avg_vel = avg_vel / np.linalg.norm(avg_vel)

    # Define resulting zero array B.
    rpoints = np.zeros(points.shape)

    # Loop over rows and determine rotated vectors.
    for idx, v in enumerate(points):
        rpoints[idx] = np.dot(R, v)

    rmedian = np.median(rpoints[:, 0])
    pmedian = np.dot(mR, [rmedian, 1.0])

    # partition
    lpart_mask = np.where(rpoints[:, 0] <= rmedian)
    rpart_mask = np.where(rpoints[:, 0] > rmedian)

    pleft  = points[lpart_mask]
    vleft  = velocities[lpart_mask]
    sl     = sign(side(np.tile(avg_vel, (len(pleft), 1)),   pleft - pmedian))
    assert np.all(np.isin(sl, [-1, 0]))

    pright = points[rpart_mask]
    vright = velocities[rpart_mask]
    sr     = sign(side(np.tile(avg_vel, (len(pright), 1)), pright - pmedian))
    assert np.all(np.isin(sr, [1]))

    d1, d2 = split_domain(domain, pmedian, avg_vel)

    l = partitioning(ax, P / 2, pleft,   vleft, domain=list(d1.edges), bisectors=bisectors + [(pmedian, avg_vel, [-1, 0])], i=2*i,   depth=depth+1)
    r = partitioning(ax, P / 2, pright, vright, domain=list(d2.edges), bisectors=bisectors + [(pmedian, avg_vel, [1])],     i=2*i+1, depth=depth+1)

    return l + r


def eval_partition(struct, points):
    workloads = np.zeros(len(struct))
    for bisectors, i in struct:
        workloads[i] = count_satisfy_all_bisectors(bisectors, points)
    return workloads


def partitioning_rcb(ax, P, points, velocities, bisectors=[], i=0, depth=0):
    N = len(points)

    width = np.amax(points, axis=0) - np.amin(points, axis=0)

    origin  = [0., 1.]
    avg_vel = [0., 0.]
    avg_vel[np.argmin(width)] = 1

    if N <= 1 or P <= 1:
        if ax:
            show_points(ax, points, s=0.5)
        return [(bisectors, i - 2**depth)]

    theta = get_angle(avg_vel, origin)

    # rotate
    R = get_rotation_matrix(theta)
    mR = get_rotation_matrix(-theta)

    avg_vel = (np.dot(mR, np.abs(np.dot(R, avg_vel))))
    avg_vel = avg_vel / np.linalg.norm(avg_vel)

    # Define resulting zero array B.
    rpoints = np.zeros(points.shape)

    # Loop over rows and determine rotated vectors.
    for idx, v in enumerate(points):
        rpoints[idx] = np.dot(R, v)

    rmedian = np.median(rpoints[:, 0])
    pmedian = np.dot(mR, [rmedian, 1.0])

    # partition
    lpart_mask = np.where(rpoints[:, 0] <= rmedian)
    rpart_mask = np.where(rpoints[:, 0] > rmedian)

    pleft = points[lpart_mask]
    vleft = velocities[lpart_mask]
    sl = sign(side(np.tile(avg_vel, (len(pleft), 1)), pleft - pmedian))
    assert np.all(np.isin(sl, [-1, 0]))

    pright = points[rpart_mask]
    vright = velocities[rpart_mask]
    sr = sign(side(np.tile(avg_vel, (len(pright), 1)), pright - pmedian))
    assert np.all(np.isin(sr, [1]))

    l = partitioning_rcb(ax, P / 2, pleft, vleft,   bisectors=bisectors + [(pmedian, avg_vel, [-1, 0])], i=2*i, depth=depth+1)
    r = partitioning_rcb(ax, P / 2, pright, vright, bisectors=bisectors + [(pmedian, avg_vel, [1])], i=2*i+1, depth=depth+1)

    return l + r
