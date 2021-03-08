import matplotlib.pyplot as plt
from .utils import *
from .math import *
from skgeom.draw import draw


def show_partitioned_points(ax, struct, points):
    colors = get_colors(len(struct))
    c = get_polygon_id(struct, points)
    c = [colors[x] for x in c]
    show_points(ax, points=points, c=c, s=0.1)
    for p, i in struct:
        for e in p.edges:
            ax.plot([e.source().x(), e.target().x()], [e.source().y(), e.target().y()], c='b', linewidth=0.1)


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
