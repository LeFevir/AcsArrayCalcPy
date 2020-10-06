# -*- coding: utf-8 -*-
"""
Module with functions for Numerical calculation of Rayleigh Integral
Includes cases with Cartesian and cylindrical coordinates

"""
import math
import numpy
import time
import logging
from multiprocessing import Pool
from datetime import timedelta


class Triangle:
    def __init__(self, x1, y1, x2, y2, x3, y3):
        self.x1 = x1
        self.y1 = y1
        self.z1 = 0.0
        self.x2 = x2
        self.y2 = y2
        self.z2 = 0.0
        self.x3 = x3
        self.y3 = y3
        self.z3 = 0.0

    def compute_centroids(self):
        self.xc = (self.x1 + self.x2 + self.x3) / 3.0
        self.yc = (self.y1 + self.y2 + self.y3) / 3.0
        self.zc = (self.z1 + self.z2 + self.z3) / 3.0

    def is_need(self, trans):
        for n in range(0, trans.num_of_elements):
            # Go to local coordinate system of each element of transducer
            l1 = trans.element_length * n
            l2 = trans.element_length * (n + 1)
            if (self.yc >= l1) & (self.yc <= l2):
                return True
        return False

    def compute_area(self):
        """
        Косое произведение векторов
        http://habrahabr.ru/post/147691/
        или вики

        S = abs( (x2-x1)(y3-y1) — (x3-x1)(y2-y1) ) / 2

        """
        self.area = abs((self.x2 - self.x1) * (self.y3 - self.y1) - (self.x3 - self.x1) * (self.y2 - self.y1)) / 2.0

        """
        # По теореме Герона
        a = math.sqrt((self.x1-self.x2)**2 + (self.y1-self.y2)**2 + (self.z1-self.z2)**2)
        b = math.sqrt((self.x1-self.x3)**2 + (self.y1-self.y3)**2 + (self.z1-self.z3)**2)
        c = math.sqrt((self.x2-self.x3)**2 + (self.y2-self.y3)**2 + (self.z2-self.z3)**2)
        p = 0.5 * (a + b + c)
        self.area = math.sqrt(p * (p - a) * (p - b) * (p - c))
        """
        return self.area


def calc_field_from_trans_opt(field, trans, medium, numberOfTrisOnLine):
    logging.info("-----Numerical field calculation has been started-----")
    logging.info("Using transducer '%s' and medium '%s'", trans.name, medium.name)
    start_time = time.clock()

    wave_length = medium.speed_of_sound / trans.frequency
    wave_number = 2 * math.pi / wave_length
    kappa = wave_number + 1j * medium.attenuation
    logging.info("Wave length = %f mm", 1.0e3 * wave_length)
    logging.info("Wave number = %f m-1", wave_number)

    # Precalculations for performance
    logging.info("Forming sources...")
    sources = form_sources_from_trans(trans, numberOfTrisOnLine)

    # Initial pressure on transducer = 1 Pa
    sources['Ss'] *= 1.0 / (medium.density * medium.speed_of_sound)

    shape = numpy.shape(field.p)
    points = [(field, sources, medium, kappa, i, j, k) for i in range(0, shape[0]) for j in range(0, shape[1]) for k in range(0, shape[2])]

    pool = Pool()

    logging.info("Entering main cycle...")
    res = pool.map(process_field_calc, points)

    res = numpy.array(res)
    field.p = numpy.reshape(res[:, 0], shape)
    field.vn = numpy.reshape(res[:, 1], shape)

    duration = time.clock() - start_time
    duration = timedelta(seconds=duration)
    logging.info("-----Numerical field calculation has been finished. It took %s", duration)


def calc_field_from_trans_after_interface_opt(field, trans, medium1, medium2, numberOfTrisOnLine):
    logging.info("-----Numerical field calculation after interface has been started-----")
    logging.info("Using transducer '%s', medium '%s' for propagation and medium '%s' after interface", trans.name, medium1.name, medium2.name)
    start_time = time.clock()

    wave_length = medium1.speed_of_sound / trans.frequency
    wave_number = 2 * math.pi / wave_length
    kappa = wave_number + 1j * medium1.attenuation
    logging.info("Wave length = %f mm", 1.0e3 * wave_length)
    logging.info("Wave number = %f m-1", wave_number)

    # Precalculations for performance
    logging.info("Forming sources...")
    sources = form_sources_from_trans(trans, numberOfTrisOnLine)

    # Initial pressure on transducer = 1 Pa
    sources['Ss'] *= 1.0 / (medium1.density * medium1.speed_of_sound)

    shape = numpy.shape(field.p)
    points = [(field, sources, medium1, medium2, kappa, i, j, k) for i in range(0, shape[0]) for j in range(0, shape[1]) for k in range(0, shape[2])]

    pool = Pool()

    logging.info("Entering main cycle...")
    res = pool.map(process_field_calc_interface, points)

    res = numpy.array(res)
    field.p = numpy.reshape(res[:, 0], shape)
    field.vn = numpy.reshape(res[:, 1], shape)

    duration = time.clock() - start_time
    duration = timedelta(seconds=duration)
    logging.info("-----Numerical field calculation has been finished. It took %s", duration)


def calc_field_from_field_opt(field, src_field, medium, frequency):
    logging.info("-----Numerical field calculation has been started-----")
    logging.info("Using field as a source and medium '%s'", medium.name)
    start_time = time.clock()

    wave_length = medium.speed_of_sound / frequency
    wave_number = 2 * math.pi / wave_length
    kappa = wave_number + 1j * medium.attenuation
    logging.info("Wave length = %f mm", 1.0e3 * wave_length)
    logging.info("Wave number = %f m-1", wave_number)

    # Precalculations for performance
    logging.info("Forming sources...")
    sources = form_sources_from_field(src_field)

    shape = numpy.shape(field.p)
    points = [(field, sources, medium, kappa, i, j, k) for i in range(0, shape[0]) for j in range(0, shape[1]) for k in range(0, shape[2])]

    pool = Pool()

    logging.info("Entering main cycle...")
    res = pool.map(process_field_calc, points)

    res = numpy.array(res)
    field.p = numpy.reshape(res[:, 0], shape)
    field.vn = numpy.reshape(res[:, 1], shape)

    duration = time.clock() - start_time
    duration = timedelta(seconds=duration)
    logging.info("-----Numerical field calculation has been finished. It took %s", duration)


def calc_field_from_field_after_interface_opt(field, src_field, medium1, medium2, frequency):
    logging.info("-----Numerical field calculation has been started-----")
    logging.info("Using field as a source, medium '%s' for propagation and medium '%s' after interface", medium1.name, medium2.name)
    start_time = time.clock()

    wave_length = medium1.speed_of_sound / frequency
    wave_number = 2 * math.pi / wave_length
    kappa = wave_number + 1j * medium1.attenuation
    logging.info("Wave length = %f mm", 1.0e3 * wave_length)
    logging.info("Wave number = %f m-1", wave_number)

    # Precalculations for performance
    logging.info("Forming sources...")
    sources = form_sources_from_field(src_field)

    shape = numpy.shape(field.p)
    points = [(field, sources, medium1, medium2, kappa, i, j, k) for i in range(0, shape[0]) for j in range(0, shape[1]) for k in range(0, shape[2])]

    pool = Pool()

    logging.info("Entering main cycle...")
    res = pool.map(process_field_calc_interface, points)

    res = numpy.array(res)
    field.p = numpy.reshape(res[:, 0], shape)
    field.vn = numpy.reshape(res[:, 1], shape)

    duration = time.clock() - start_time
    duration = timedelta(seconds=duration)
    logging.info("-----Numerical field calculation has been finished. It took %s", duration)


def process_field_calc(point):
    field = point[0]
    sources = point[1]
    medium = point[2]
    kappa = point[3]
    i = point[4]
    j = point[5]
    k = point[6]

    field_x, field_y, field_z = field.get_cartesian_coords(i, j, k)
    field_normal_x, field_normal_y, field_normal_z = field.normal(i, j, k)

    rx, ry, rz, r = calc_distance_vectors(sources, field_x, field_y, field_z)
    coses = calc_coses(rx, ry, rz, r, field_normal_x, field_normal_y, field_normal_z)

    # check if r crosses the source
    bad_rays_map = check_crosses(rx, ry, rz, sources)

    p, vn = sum_up_sources_rayleigh(kappa, r, sources['Ss'] * bad_rays_map, 1.0, coses)

    p *= -1j * kappa * medium.density * medium.speed_of_sound / (2.0 * math.pi)
    vn *= -1j * kappa / (2.0 * math.pi)

    return p, vn


def process_field_calc_interface(point):
    field = point[0]
    sources = point[1]
    medium1 = point[2]
    medium2 = point[3]
    kappa = point[4]
    i = point[5]
    j = point[6]
    k = point[7]

    field_x, field_y, field_z = field.get_cartesian_coords(i, j, k)
    field_normal_x, field_normal_y, field_normal_z = field.normal(i, j, k)

    rx, ry, rz, r = calc_distance_vectors(sources, field_x, field_y, field_z)
    incident_angles = calc_incident_angles(rx, ry, rz, r, field_normal_x, field_normal_y, field_normal_z)
    transmission_coses = calc_transmission_coses(incident_angles, medium1, medium2)

    transmission_coefs_p = calc_transmission_coefs_for_pressure(medium1, medium2, incident_angles, transmission_coses)
    transmission_coefs_vn = calc_transmission_coefs_for_velocity(medium1, medium2, incident_angles, transmission_coses)

    # check if r crosses the source
    bad_rays_map = check_crosses(rx, ry, rz, sources)

    p, vn = sum_up_sources_rayleigh(kappa, r, sources['Ss'] * bad_rays_map, transmission_coefs_p, transmission_coefs_vn * transmission_coses)

    p *= -1j * kappa * medium1.density * medium1.speed_of_sound / (2 * math.pi)
    vn *= -1j * kappa / (2 * math.pi)

    return p, vn


def form_target_surface(field):
    surface_Xs = []
    surface_Ys = []
    surface_Zs = []
    surface_normal_x = []
    surface_normal_y = []
    surface_normal_z = []

    surface_ps = []
    surface_vns = []

    shape = numpy.shape(field.p)
    for i in range(0, shape[0]):
        for j in range(0, shape[1]):
            for k in range(0, shape[2]):
                field_x, field_y, field_z = field.get_cartesian_coords(i, j, k)
                norm_x, norm_y, norm_z = field.normal(i, j, k)

                surface_Xs.append(field_x)
                surface_Ys.append(field_y)
                surface_Zs.append(field_z)
                surface_normal_x.append(norm_x)
                surface_normal_y.append(norm_y)
                surface_normal_z.append(norm_z)
                surface_ps.append(field.p[i, j, k])
                surface_vns.append(field.vn[i, j, k])

    surface_Xs = numpy.array(surface_Xs)
    surface_Ys = numpy.array(surface_Ys)
    surface_Zs = numpy.array(surface_Zs)
    surface_normal_Xs = numpy.array(surface_normal_x)
    surface_normal_Ys = numpy.array(surface_normal_y)
    surface_normal_Zs = numpy.array(surface_normal_z)
    surface_ps = numpy.array(surface_ps)
    surface_vns = numpy.array(surface_vns)

    surface = {'Xs': surface_Xs,
                'Ys': surface_Ys,
                'Zs': surface_Zs,
                'n_Xs': surface_normal_Xs,
                'n_Ys': surface_normal_Ys,
                'n_Zs': surface_normal_Zs,
                'ps': surface_ps,
                'vns': surface_vns}
    return surface


def calc_field_from_trans(field, trans, medium, numberOfTrisOnLine):
    logging.info("-----Numerical field calculation has been started-----")
    logging.info("Using transducer '%s' and medium '%s'", trans.name, medium.name)
    start_time = time.clock()

    wave_length = medium.speed_of_sound / trans.frequency
    wave_number = 2 * math.pi / wave_length
    kappa = wave_number + 1j * medium.attenuation
    logging.info("Wave length = %f mm", 1.0e3 * wave_length)
    logging.info("Wave number = %f m-1", wave_number)

    # Precalculations for performance
    logging.info("Forming sources...")
    sources = form_sources_from_trans(trans, numberOfTrisOnLine)

    # Initial pressure on transducer = 1 Pa
    sources['Ss'] *= 1.0 / (medium.density * medium.speed_of_sound)

    field_calc_main_cycle(field, sources, medium, kappa)

    duration = time.clock() - start_time
    duration = timedelta(seconds=duration)
    logging.info("-----Numerical field calculation has been finished. It took %s", duration)


def calc_field_from_field(field, src_field, medium, frequency):
    logging.info("-----Numerical field calculation has been started-----")
    logging.info("Using field as a source and medium '%s'", medium.name)
    start_time = time.clock()

    wave_length = medium.speed_of_sound / frequency
    wave_number = 2 * math.pi / wave_length
    kappa = wave_number + 1j * medium.attenuation
    logging.info("Wave length = %f mm", 1.0e3 * wave_length)
    logging.info("Wave number = %f m-1", wave_number)

    # Precalculations for performance
    logging.info("Forming sources...")
    sources = form_sources_from_field(src_field)

    field_calc_main_cycle(field, sources, medium, kappa)

    duration = time.clock() - start_time
    duration = timedelta(seconds=duration)
    logging.info("-----Numerical field calculation has been finished. It took %s", duration)


def calc_field_from_trans_after_interface(field, trans, medium1, medium2, numberOfTrisOnLine):

    logging.info("-----Numerical field calculation after interface has been started-----")
    logging.info("Using transducer '%s', medium '%s' for propagation and medium '%s' after interface", trans.name, medium1.name, medium2.name)
    start_time = time.clock()

    wave_length = medium1.speed_of_sound / trans.frequency
    wave_number = 2 * math.pi / wave_length
    kappa = wave_number + 1j * medium1.attenuation
    logging.info("Wave length = %f mm", 1.0e3 * wave_length)
    logging.info("Wave number = %f m-1", wave_number)

    # Precalculations for performance
    logging.info("Forming sources...")
    sources = form_sources_from_trans(trans, numberOfTrisOnLine)

    # Initial pressure on transducer = 1 Pa
    sources['Ss'] *= 1.0 / (medium1.density * medium1.speed_of_sound)

    field_calc_after_interface_main_cycle(field, sources, medium1, medium2, kappa)

    duration = time.clock() - start_time
    duration = timedelta(seconds=duration)
    logging.info("-----Numerical field calculation has been finished. It took %s", duration)


def calc_field_from_field_after_interface(field, src_field, medium1, medium2, frequency):

    logging.info("-----Numerical field calculation has been started-----")
    logging.info("Using field as a source, medium '%s' for propagation and medium '%s' after interface", medium1.name, medium2.name)
    start_time = time.clock()

    wave_length = medium1.speed_of_sound / frequency
    wave_number = 2 * math.pi / wave_length
    kappa = wave_number + 1j * medium1.attenuation
    logging.info("Wave length = %f mm", 1.0e3 * wave_length)
    logging.info("Wave number = %f m-1", wave_number)

    # Precalculations for performance
    logging.info("Forming sources...")
    sources = form_sources_from_field(src_field)

    field_calc_after_interface_main_cycle(field, sources, medium1, medium2, kappa)

    duration = time.clock() - start_time
    duration = timedelta(seconds=duration)
    logging.info("-----Numerical field calculation has been finished. It took %s", duration)


def calc_power(field):
    """Calc power as a sum of Poynting vectors (P = Re(p * conj(v))) in each point and multiplyies it with area"""

    poynt_in_points = numpy.real(field.p * numpy.conj(field.vn))
    power = numpy.sum(poynt_in_points)
    power *= field.one_pixel_area
    return power


def field_calc_main_cycle(field, sources, medium, kappa):
    logging.info("Entering main cycle...")
    shape = numpy.shape(field.p)
    for i in range(0, shape[0]):
        for j in range(0, shape[1]):
            for k in range(0, shape[2]):
                field_x, field_y, field_z = field.get_cartesian_coords(i, j, k)
                field_normal_x, field_normal_y, field_normal_z = field.normal(i, j, k)

                rx, ry, rz, r = calc_distance_vectors(sources, field_x, field_y, field_z)
                coses = calc_coses(rx, ry, rz, r, field_normal_x, field_normal_y, field_normal_z)

                # check if r crosses the source
                bad_rays_map = check_crosses(rx, ry, rz, sources)

                p, vn = sum_up_sources_rayleigh(kappa, r, sources['Ss'] * bad_rays_map, 1.0, coses)

                p *= -1j * kappa * medium.density * medium.speed_of_sound / (2.0 * math.pi)
                vn *= -1j * kappa / (2.0 * math.pi)

                field.p[i, j, k] = p
                field.vn[i, j, k] = vn


def field_calc_after_interface_main_cycle(field, sources, medium1, medium2, kappa):
    logging.info("Entering main cycle...")
    shape = numpy.shape(field.p)
    for i in range(0, shape[0]):
        for j in range(0, shape[1]):
            for k in range(0, shape[2]):
                field_x, field_y, field_z = field.get_cartesian_coords(i, j, k)
                field_normal_x, field_normal_y, field_normal_z = field.normal(i, j, k)

                rx, ry, rz, r = calc_distance_vectors(sources, field_x, field_y, field_z)
                incident_angles = calc_incident_angles(rx, ry, rz, r, field_normal_x, field_normal_y, field_normal_z)
                transmission_coses = calc_transmission_coses(incident_angles, medium1, medium2)

                transmission_coefs_p = calc_transmission_coefs_for_pressure(medium1, medium2, incident_angles, transmission_coses)
                transmission_coefs_vn = calc_transmission_coefs_for_velocity(medium1, medium2, incident_angles, transmission_coses)

                # check if r crosses the source
                bad_rays_map = check_crosses(rx, ry, rz, sources)

                p, vn = sum_up_sources_rayleigh(kappa, r, sources['Ss'] * bad_rays_map, transmission_coefs_p, transmission_coefs_vn * transmission_coses)

                p *= -1j * kappa * medium1.density * medium1.speed_of_sound / (2 * math.pi)
                vn *= -1j * kappa / (2 * math.pi)

                field.p[i, j, k] = p
                field.vn[i, j, k] = vn


def slice_trans_on_tris(trans, numberOfTrisOnLine):
    """ Slice the transducer into triangles """
    logging.info("Slicing of transducer into triangles has been started. Number of triangles on line %d", numberOfTrisOnLine)

    numberOfRectsOnLine = numberOfTrisOnLine / 2
    dx = trans.element_width / numberOfRectsOnLine
    dy = trans.element_length / numberOfRectsOnLine
    a = trans.element_width / 2.0
    tris = []
    triArea = 0.0  # Считаем площадь треугольника только один раз

    for j in range(0, trans.num_of_elements * numberOfRectsOnLine):
        for i in range(0, numberOfRectsOnLine):

            tri1 = Triangle(-a+i*dx, j*dy, -a+i*dx, (j+1)*dy, -a+(i+1)*dx, j*dy)
            tri2 = Triangle(-a+(i+1)*dx, j*dy, -a+i*dx, (j+1)*dy, -a+(i+1)*dx, (j+1)*dy)

            tri1.compute_centroids()

            if tri1.is_need(trans):
                if triArea == 0.0:
                    triArea = tri1.compute_area()
                else:
                    tri1.area = triArea

                tris.append(tri1)

                # Для второго треугольника можно не проверять,
                # так как его центр по Y совпадает с первым треуглольником,
                # а проверка попадания проводится по Y
                tri2.compute_centroids()
                tri2.area = triArea
                tris.append(tri2)

    logging.info("Slicing of transducer into triangles has been finished. Number of active triangles: %d. Sides of triangle are %fx%f mm. Area of each triangle is %f mm^2", len(tris), dx/1.0e-03, dy/1.0e-03,  triArea*1.0e+06)

    return tris


def form_sources_from_tris(tris):
    sources_Xs = numpy.array([tri.xc for tri in tris])
    sources_Ys = numpy.array([tri.yc for tri in tris])
    sources_Zs = numpy.array([tri.zc for tri in tris])
    sources_Ss = numpy.array([tri.area for tri in tris])

    # These values are valid only in this case of flat transducer's surface!
    sources_normal_Xs = numpy.zeros_like(sources_Xs)
    sources_normal_Ys = numpy.zeros_like(sources_Xs)
    sources_normal_Zs = numpy.ones_like(sources_Xs)

    sources = {'Xs': sources_Xs,
                'Ys': sources_Ys,
                'Zs': sources_Zs,
                'Ss': sources_Ss,
                'n_Xs': sources_normal_Xs,
                'n_Ys': sources_normal_Ys,
                'n_Zs': sources_normal_Zs}
    return sources


def form_sources_from_trans(trans, numberOfTrisOnLine):
    tris = slice_trans_on_tris(trans, numberOfTrisOnLine)
    return form_sources_from_tris(tris)


def form_sources_from_field(field):
    sources_Xs = []
    sources_Ys = []
    sources_Zs = []
    sources_Ss = []
    sources_normal_x = []
    sources_normal_y = []
    sources_normal_z = []

    shape = numpy.shape(field.p)
    for i in range(0, shape[0]):
        for j in range(0, shape[1]):
            for k in range(0, shape[2]):
                field_x, field_y, field_z = field.get_cartesian_coords(i, j, k)
                norm_x, norm_y, norm_z = field.normal(i, j, k)

                sources_Xs.append(field_x)
                sources_Ys.append(field_y)
                sources_Zs.append(field_z)
                sources_Ss.append(field.vn[i, j, k])
                sources_normal_x.append(norm_x)
                sources_normal_y.append(norm_y)
                sources_normal_z.append(norm_z)

    sources_Xs = numpy.array(sources_Xs)
    sources_Ys = numpy.array(sources_Ys)
    sources_Zs = numpy.array(sources_Zs)
    sources_Ss = numpy.array(sources_Ss)
    sources_Ss *= field.one_pixel_area

    sources_normal_Xs = numpy.array(sources_normal_x)
    sources_normal_Ys = numpy.array(sources_normal_y)
    sources_normal_Zs = numpy.array(sources_normal_z)

    sources = {'Xs': sources_Xs,
                'Ys': sources_Ys,
                'Zs': sources_Zs,
                'Ss': sources_Ss,
                'n_Xs': sources_normal_Xs,
                'n_Ys': sources_normal_Ys,
                'n_Zs': sources_normal_Zs}
    return sources


def calc_distance_vectors(sources, field_x, field_y, field_z):
    DX = field_x - sources['Xs']
    DY = field_y - sources['Ys']
    DZ = field_z - sources['Zs']
    return DX, DY, DZ, numpy.sqrt(DX ** 2 + DY ** 2 + DZ ** 2)


def calc_coses(rx, ry, rz, r, normal_x, normal_y, normal_z):
    coses = numpy.abs(rx * normal_x + ry * normal_y + rz * normal_z) / r
    return coses


def check_crosses(rx, ry, rz, sources):
    # Checks if r crosses the source
    # Checking Scalar multiplying normal vector on distance vector
    # if (mult >= 0) => (angle <= 90) => OK, unit value, else (angle > 90) => Zero value
    mult = rx * sources['n_Xs'] + ry * sources['n_Ys'] + rz * sources['n_Zs']
    return 1.0 * (mult >= 0)


def calc_incident_angles(rx, ry, rz, r, normal_x, normal_y, normal_z):
    coses = calc_coses(rx, ry, rz, r, normal_x, normal_y, normal_z)
    return numpy.arccos(coses)


def calc_transmission_coses(incident_angles, medium1, medium2):
    val = (1 - ((medium2.speed_of_sound / medium1.speed_of_sound) * numpy.sin(incident_angles))**2)
    return numpy.sqrt(val * (val > 0.0))


def calc_transmission_coefs_for_velocity(medium1, medium2, incident_angles, transmission_coses):
    incident_local_impedances = calc_incident_local_impedances(medium1, incident_angles)
    transmission_local_impedances = calc_transmission_local_impedances(medium2, transmission_coses)

    return ((medium1.density * medium1.speed_of_sound) / (medium2.density * medium2.speed_of_sound)) * 2.0 * transmission_local_impedances / (transmission_local_impedances + incident_local_impedances)


def calc_incident_local_impedances(medium, angles):
    return (medium.density * medium.speed_of_sound) / numpy.cos(angles)


def calc_transmission_local_impedances(medium, transmission_coses):
    return (transmission_coses > 0) * (medium.density * medium.speed_of_sound) / (transmission_coses * (transmission_coses > 0) + 1.0 * (transmission_coses == 0))


def calc_transmission_coefs_for_pressure(medium1, medium2, icident_angles, transmission_coses):
    incident_local_impedances = calc_incident_local_impedances(medium1, icident_angles)
    transmission_local_impedances = calc_transmission_local_impedances(medium2, transmission_coses)

    return 2.0 * transmission_local_impedances / (transmission_local_impedances + incident_local_impedances)


def sum_up_sources_rayleigh(kappa, r, S, tran_coef_p, tran_coef_vn):
    val = numpy.exp(1j * kappa * r) * (1 - 1j / (kappa * r)) * (S / r)
    # val = numpy.exp(1j * kappa * r) * (S / r)
    p = val * tran_coef_p
    vn = val * tran_coef_vn
    return numpy.sum(p), numpy.sum(vn)


if __name__ == '__main__':
    pass
    # test case
    # trans = transducer.Transducer("Probe1", 1, 4.0e-03, 5.0e-03, 6.0e06)
    # tris = sliceTransOnTrisProper(trans, 100)

    # trans = transducer.Transducer(
    #   name = "Probe1",
    #   num_of_elements = 1,
    #   element_width = 4.0e-03,
    #   element_length = 5.0e-03,
    #   # frequency = 6.0e06
    #   # frequency = 12.0e06
    #   # frequency = 4.0e06
    #   frequency = 11.6e06
    #   )

    # water = medium.Medium(
    #   name = "Water",
    #   density = 1000.0,
    #   speed_of_sound = 1488.0
    #   )
