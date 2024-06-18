import numpy as np
import random

from I3Tray import I3Units
from icecube import dataclasses


### Sampling utilities ###


def random_sign(u):
    """
    Sample sign (-1/+1).
    Parameters
    ----------
    u : float
        Random variable in [0,1].
    Returns
    -------
    sign : float
    """
    return 1 if u < 0.5 else -1


def sample_power_law(k_min, k_max, gamma, u):
    """
    Sample from power law distribution k^(-gamma).
    Parameters
    ----------
    k_min : float
        Lower distribution boundary
    k_max : float
        Upper distribution boundary
    gamma : float
        Power law
    u : float
        Random variable in [0,1].
    Returns
    -------
    energy : float
    """
    a = (k_max ** (-gamma + 1.0) - k_min ** (-gamma + 1.0)) * u
    b = k_min ** (-gamma + 1.0)
    c = 1.0 / (-gamma + 1.0)
    return (a + b) ** c


def sample_exp_in_range(a, b, u, rate = 0.01):
    """
    Sample from truncated exponential on domain.
    Parameters
    ----------
    a : float
        Lower distribution boundary
    b : float
        Upper distribution boundary
    u : float
        Random variable in [0,1].
    gamma : float
        Power law
    Returns
    -------
    value : float
    """
    exp_rate_a = np.exp(-rate * a)
    return -np.log(exp_rate_a - u * (exp_rate_a - np.exp(-rate * b))) / rate


### IceCube related functions ###


def direction_vector(zenith, azimuth):
    """
    Calculate direction vector from zenith/azimuth.
    Parameters
    ----------
    zenith : float
        zenith of direction [rad]
    azimuth : float
        azimuth of direction [rad]
    Returns
    -------
    direction : array of flots with shape (3,)
        direction vector in cartesian coordinates
    """
    direction = np.array(
        [
            -np.sin(zenith) * np.cos(azimuth),  # x
            -np.sin(zenith) * np.sin(azimuth),  # y
            -np.cos(zenith),  # z
        ]
    )
    return direction


def direction_angles(x, y, z):
    """
    Calculate direction angles from cartesian vector.
    Parameters
    ----------
    x : float
        x of direction
    y : float
        y of direction
    z : float
        y of direction
    Returns
    -------
    direction : array of floats with shape (2,)
        direction zenith and azimuth
    """
    if x==0.0 and y == 0.0:
        arctan_arg = 0.0  # azimuth does not mean anything here..
    elif x==0.0 and y > 0.0:
        arctan_arg = +np.inf
    elif x==0.0 and y < 0.0:
        arctan_arg = -np.inf
    else:
        arctan_arg = y/x
    azimuth_tmp = np.where(x < 0., np.arctan(arctan_arg), np.arctan(arctan_arg) + np.pi)
    direction = np.array(
        [
            np.arccos(-z),  # zenith
            np.where(azimuth_tmp < 0., azimuth_tmp+2*np.pi, azimuth_tmp),  # azimuth
        ]
    )

    return direction


def make_cascade(position, direction, energy, time):
    """
    Make an EM cascade.
    Parameters
    ----------
    position : list of floats with shape (3,)
        Position (x, y, z) [m].
    direction : list of floats with shape (2,)
        Direction (zenith, azimuth) [rad].
    energy : float
        Energy [GeV].
    time : float
        Time [ns].
    Returns
    -------
    particle : I3Particle
        The created cascade particle.
    """

    p = dataclasses.I3Particle(
        dataclasses.I3Particle.ParticleShape.Cascade,
        dataclasses.I3Particle.ParticleType.EMinus,
    )

    p.pos = dataclasses.I3Position(position[0], position[1], position[2])
    p.dir = dataclasses.I3Direction(direction[0], direction[1])

    p.energy = energy * I3Units.GeV
    p.time = time * I3Units.ns
    p.speed = dataclasses.I3Constants.c

    #     important for clsim
    p.location_type = dataclasses.I3Particle.LocationType.InIce
    return p


def add_cascades_to_tree(frame, tree, positions, directions, energies, times):
    """
    Add cascade(s) to the tree.
    Parameters
    ----------
    frame : frame
        Full frame to add keys directly.
    tree : I3MCTree
        I3MCTree to add the cascade(s) to.
    positions : list of floats with shape (n,3)
        Positions (x, y, z) [m].
    directions : list of floats with shape (n,2)
        Directions (zenith, azimuth) [rad].
    energies : list of floats with shape (n)
        Energies [GeV].
    times : list of floats with shape (n)
        Times [ns].
    Returns
    -------
    tree : I3MCTree
        I3MCTree with added cascade(s).
    """

    primary = tree.get_primaries()[0]

    for count in range(len(energies)):
        particle = make_cascade(
            position=positions[count],
            direction=directions[count],
            energy=energies[count],
            time=times[count],
        )
        tree.append_child(primary, particle)
        frame["true_cascade_{}".format(count)] = particle
    return tree
