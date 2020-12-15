import numpy as np
import matplotlib.cm as cm


def rotate(l, n):
    return l[n:] + l[:n]


def get_colors(N):
    evenly_spaced_interval = np.linspace(0, 1, N)
    colors = [cm.rainbow(x) for x in evenly_spaced_interval]
    return colors


def rand_negpos(N):
    r = np.random.rand(N,2) < 0.5
    return np.select([r, ~r], [np.ones(r.shape), -np.ones(r.shape)])


def gen_points_on_diagonal(n):
    l  =  1.0
    x0 = -0.5
    y0 = -0.5
    dx = l / n
    dy = l / n
    points = [[x0, y0]]
    for i in range(1, n):
        points.append([x0 + i*dx, y0 + i*dy])
    return np.asarray(points)


def uniform_circle(num_samples, R=1):
    p = np.ones((num_samples, 2))
    # generate the points
    r2 = R * np.sqrt(np.random.rand(num_samples, 1))
    theta2 = 2 * np.pi * np.random.rand(num_samples, 1)
    p[:, 0] = (r2 * np.cos(theta2)).flatten()
    p[:, 1] = (r2 * np.sin(theta2)).flatten()
    return p


def expansion_velocity(p, c=np.zeros((1, 2))):
    return p - c


def random_biased_velocity(N):
    V    = np.random.rand(N, 2) - 0.5
    bias = np.random.rand(N, 2)
    return V + bias
