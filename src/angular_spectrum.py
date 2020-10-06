# -*- coding: utf-8 -*-
import numpy
import math
from scipy.fftpack import fft2, ifft2, fftshift, ifftshift

"""
17.10.2012 Sergey Ilyin

Module for calculations with Angular Spectrum Approach.

"""


def propagate_cartesian(p_0, dx, dy, f, ro, c, z):
    """
    Calculates pressure and velocity distributions (X and Y) on the plane Z = z
    from initial pressure distribution on plan Z = 0.0 in Cartesian coordinate system.

    """

    # Wave number
    k0 = 2 * math.pi * f / c

    s = numpy.shape(p_0)
    # Building kx, ky vectors
    kx = numpy.linspace(-math.pi / dx, math.pi / dx, s[0])
    ky = numpy.linspace(-math.pi / dy, math.pi / dy, s[1])
    ky, kx = numpy.meshgrid(ky, kx)

    # Calculating Angular spectrum
    p0_kx_ky = fftshift(fft2(p_0))

    # Calculating kz and propation term H
    kz2 = k0 ** 2 - kx ** 2 - ky ** 2
    kz = numpy.sqrt(kz2 * (kz2 > 0))
    H = numpy.exp(-z * kz * 1j) * (kz2 > 0)

    # Propagate anglular spectrum
    p_kx_ky = p0_kx_ky * H

    # Pressure from angular spectrum
    p = ifft2(ifftshift(p_kx_ky))

    # Calculating a multiplier for Velocity
    mult = 1 - (kx ** 2 + ky ** 2) / k0 ** 2
    mult = numpy.sqrt(mult * (mult > 0)) / ro / c

    # Calculating Velocity
    v = ifft2(ifftshift(p_kx_ky * mult))

    return p, v


if __name__ == '__main__':
    pass
