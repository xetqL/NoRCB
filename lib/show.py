import matplotlib.pyplot as plt
from .utils import *
from .math import *


def show_partitioned_points(ax, struct, points):
    colors = get_colors(len(struct))
    for bisectors, i in struct:
        who = who_satisfy_all_bisectors(bisectors, points)
        c   = [colors[i]] * len(points[who])
        show_points(ax, points=points[who], c=c, s=0.5)


def show_points(ax, points, **kwargs):
    return ax.scatter(points[:, 0], points[:, 1], **kwargs)


def show_velocity(ax, p, velocity, **kwargs):
    ax.arrow(x = p[0], y = p[1], dx=velocity[0], dy=velocity[1], **kwargs)


def show_velocities(ax, p, velocities, **kwargs):
    for velocity in velocities:
        show_velocity(ax, p, velocity, **kwargs)


def draw_bisector(ax, p, v, c='black'):
    print(p, v)
    ax.plot((p[0], p[0] + 100 * v[0]), (p[1], p[1] + 100 * v[1]), linewidth=0.6, c=c)
    ax.plot((p[0], p[0] - 100 * v[0]), (p[1], p[1] - 100 * v[1]), linewidth=0.6, c=c)


def draw_bisectors(ax, bisectors):
    for p, v, _ in bisectors:
        draw_bisector(ax, p, v)
