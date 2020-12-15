import skgeom as sg
from itertools import combinations

from .utils import *



# use convex hull to get the segments
def get_segments_from_vertices(vertices):
    points  = sg.convex_hull.graham_andrew(vertices)
    rpoints = rotate(points, -1)
    spoints = zip(points, rpoints)
    segments= []
    for p1, p2 in spoints:
        segments.append(sg.Segment2(p1, p2))
    return segments


def find_intersections_with_bisectors(target, others):
    # build target point
    medtarget     =  sg.Point2(*target[0])
    raytarget     = [sg.Ray2(medtarget, sg.Point2(*target[1])), sg.Ray2(medtarget, sg.Point2(*(-target[1])))]
    intersections = []
    # check all the bisectors
    for bisector in others:
        median    = sg.Point2(*bisector[0])
        rays      = [sg.Ray2(median, sg.Point2(*(bisector[1]))), sg.Ray2(median, sg.Point2(*(-bisector[1])))]
        # intersections for pos and negative ray
        print(())
        for tar, other in list(combinations(raytarget + rays, 2)):
            i = sg.intersection(tar, other)
            if i is not None:
                intersections.append(i)
                break
    return intersections


def is_vertex(P, constraints):
    return does_satisfy_all_bisectors(constraints, P)


def compute_imbalance(workloads):
    return np.max(workloads) / np.mean(workloads) - 1.


def who_satisfy_all_bisectors(bisectors, points):
    N = len(points)
    s = np.ones(len(points)).astype(bool)
    for pmedian, avg_vel, result in bisectors:
        s = s & np.isin(sign(side(np.tile(avg_vel, (N, 1)), points - pmedian)), result)
    return s


def count_satisfy_all_bisectors(bisectors, points):
    return np.count_nonzero(who_satisfy_all_bisectors(bisectors, points))


def does_satisfy_all_bisectors(bisectors, p):
    N = 1
    s = np.ones(1).astype(bool)
    for pmedian, avg_vel, result in bisectors:
        s = s & np.isin(sign(side(np.tile(avg_vel, (N, 1)), p - pmedian)), result)
    return s


def get_average_velocity(velocities, axis):
    v = np.copy(velocities)
    R = get_rotation_matrix(np.pi)
    # print(np.inner(R, v))
    msk = np.where(v[:, axis] >= 0)
    v[msk] = np.inner(R, v[msk]).T
    return np.mean(v, axis=0)


def get_rotation_matrix(theta):
    c, s = np.cos(theta), np.sin(theta)
    R = np.array(((c, -s), (s, c)))
    return R


def degree(rad):
    return rad * 180. / (np.pi)


def radian(deg):
    return (np.pi * deg) / 180.


def normalize(v):
    norm=np.linalg.norm(v, ord=1)
    if norm==0:
        norm=np.finfo(v.dtype).eps
    return v/norm


def round_to_rotated_axis(rvel):
    rvel = np.abs(rvel)
    if np.allclose(rvel, [0, 1]):
        return np.asarray([1, 0])
    elif np.allclose(rvel, [1, 0]):
        return np.asarray([0, 1])
    else:
        return None


def side1(T, D):
    return D[0] * T[1] - D[1] * T[0]


def side(T, D):
    return D[:, 0] * T[:, 1] - D[:, 1] * T[:, 0]


def sign(T):
    zero = np.zeros_like(T)
    none = -np.ones_like(T)
    pone = np.ones_like(T)
    msk_zero = np.isclose(T, zero)
    T = np.select([msk_zero,  ~msk_zero], [zero, T])
    msk_neg  = T < 0
    msk_pos  = T > 0
    T = np.select([msk_neg, msk_pos], [none, pone], default=0)
    return T


def get_angle(v, origin):
    dot = np.dot(v, origin)
    det = np.linalg.det([v, origin])
    return np.arctan2(det, dot)
