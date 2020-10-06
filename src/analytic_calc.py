# -*- coding: utf-8 -*-

"""
Copyright 2012. Sergey Ilin.
Lab 366, Acoustic Department, Faculty of Physics
Lomonosov Moscow State University

"""

import math
import numpy
import time
import logging
from scipy.special import j1
from multiprocessing import Pool
from datetime import timedelta


def calc_exact_on_axis(field, trans, medium):
    """
    Calculation field on axis from the circle piston transducer using the formula
    p/p0 = -2i*exp(i*k/2*(sqrt(a*a+z*z)+z))*sin(k/2*(sqrt(a*a+z*z)-z))

    a - radius of the trans
    k - wave number
    z - axis coord
    """

    logging.info("-----Exact field on axis calculation has been started-----")
    logging.info("Using transducer '%s' and medium '%s'", trans.name, medium.name)
    start_time = time.clock()

    wave_length = medium.speed_of_sound / trans.frequency
    wave_number = 2.0 * math.pi / wave_length
    # kappa = wave_number + 1j * medium.attenuation
    logging.info("Wave length = %f mm", 1.0e3 * wave_length)
    logging.info("Wave number = %f m-1", wave_number)

    # Precalculations for performance
    wave_number_half = wave_number / 2.0
    a2 = trans.element_radius*trans.element_radius

    # z = field.z[0,0]

    field.p[0,0,:] = -2*1j*numpy.exp(1j*wave_number_half*(numpy.sqrt(a2+field.z*field.z)+field.z))*numpy.sin(wave_number_half*(numpy.sqrt(a2+field.z*field.z)-field.z))
    field.vn[0,0,:] = field.p / (medium.density * medium.speed_of_sound)

    # print(numpy.shape(field.p))

    duration = time.clock() - start_time
    duration = timedelta(seconds=duration)
    logging.info("-----Exact field on axis calculation has been finished. It took %s", duration)


def calc_field_analytically(field, trans, medium):
    logging.info("-----Analitical field calculation has been started-----")
    logging.info("Using transducer '%s' and medium '%s'", trans.name, medium.name)
    start_time = time.clock()

    wave_length = medium.speed_of_sound / trans.frequency
    wave_number = 2 * math.pi / wave_length
    # kappa = wave_number + 1j * medium.attenuation
    logging.info("Wave length = %f mm", 1.0e3 * wave_length)
    logging.info("Wave number = %f m-1", wave_number)

    # Precalculations for performance
    logging.info("Forming sources...")
    sources = form_sources_from_array(trans)

    shape = numpy.shape(field.p)
    points = [(field, sources, trans.element_radius, wave_number, i, j, k, medium.density, medium.speed_of_sound) for i in range(0, shape[0]) for j in range(0, shape[1]) for k in range(0, shape[2])]

    pool = Pool()

    logging.info("Entering main cycle...")
    res = pool.map(process_field_calc, points)

    res = numpy.array(res)
    field.p = numpy.reshape(res[:, 0], shape)
    field.vn = numpy.reshape(res[:, 1], shape)

    duration = time.clock() - start_time
    duration = timedelta(seconds=duration)
    logging.info("-----Analitical field calculation has been finished. It took %s", duration)


def form_sources_from_array(trans):
    sources_Xs = numpy.array([element['center_x'] for element in trans.elements])
    sources_Ys = numpy.array([element['center_y'] for element in trans.elements])
    sources_Zs = numpy.array([element['center_z'] for element in trans.elements])
    sources_normal_Xs = numpy.array([element['normal_x'] for element in trans.elements])
    sources_normal_Ys = numpy.array([element['normal_y'] for element in trans.elements])
    sources_normal_Zs = numpy.array([element['normal_z'] for element in trans.elements])
    sources_phase_shift = numpy.array([element['phase_shift'] for element in trans.elements])

    sources = {'Xs': sources_Xs,
                'Ys': sources_Ys,
                'Zs': sources_Zs,
                'n_Xs': sources_normal_Xs,
                'n_Ys': sources_normal_Ys,
                'n_Zs': sources_normal_Zs,
                'phase_shift': sources_phase_shift}

    return sources


def process_field_calc(point):
    field = point[0]
    sources = point[1]
    element_radius = point[2]
    wave_number = point[3]
    i = point[4]
    j = point[5]
    k = point[6]
    density = point[7]
    speed_of_sound = point[8]

    field_x, field_y, field_z = field.get_cartesian_coords(i, j, k)
    field_normal_x, field_normal_y, field_normal_z = field.normal(i, j, k)

    rx, ry, rz, r = calc_distance_vectors(sources, field_x, field_y, field_z)
    
    sines = calc_sines(sources, rx, ry, rz, r)
    coses = calc_coses(rx, ry, rz, r, field_normal_x, field_normal_y, field_normal_z)

    # check if r crosses the source
    # bad_rays_map = check_crosses(rx, ry, rz, sources)

    p, vn = sum_up_analyt(wave_number, r, sines, element_radius, sources['phase_shift'], coses / (density * speed_of_sound))
    return p, vn


def calc_distance_vectors(sources, field_x, field_y, field_z):
    DX = field_x - sources['Xs']
    DY = field_y - sources['Ys']
    DZ = field_z - sources['Zs']
    return DX, DY, DZ, numpy.sqrt(DX ** 2 + DY ** 2 + DZ ** 2)


def calc_sines(sources, rx, ry, rz, r):
    # Using cross product
    # |AxB| = |A| * |B| * sin(A,B)
    # |AxB| = |(Ay*Bz - Az*By, Az*Bx - Ax*Bz, Ax*By - Ay*Bx)|
    module = numpy.sqrt((sources['n_Ys'] * rz - sources['n_Zs'] * ry)**2 + (sources['n_Zs'] * rx - sources['n_Xs'] * rz)**2 + (sources['n_Xs'] * ry - sources['n_Ys'] * rx)**2)
    sines = module / r
    return sines


def calc_coses(rx, ry, rz, r, normal_x, normal_y, normal_z):
    coses = numpy.abs(rx * normal_x + ry * normal_y + rz * normal_z) / r
    return coses


def sum_up_analyt(wave_number, r, sines, element_radius, phase_shift, coef_vn):
    a = element_radius
    val = numpy.exp(1j * wave_number * r) * (1 / (2.0 * r))
    bessel_val = 2 * j1(wave_number * a * sines) / (wave_number * a * sines + 1.0 * (sines == 0.0)) + 1.0 * (sines == 0.0)
    val *= bessel_val * (-1j * wave_number * element_radius * element_radius)
    p = val * numpy.exp(-1j * wave_number * phase_shift)
    vn = p * coef_vn
    return numpy.sum(p), numpy.sum(vn)
