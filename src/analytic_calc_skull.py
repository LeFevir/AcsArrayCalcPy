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
    # PHASE SHIFT HERE
    sources = form_sources_from_array(trans)

    skull_info = load_skull_info_from_matlab()

    shape = numpy.shape(field.p)
    points = [(field, sources, trans.element_radius, wave_number, i, j, k, medium.density, medium.speed_of_sound, skull_info) for i in range(0, shape[0]) for j in range(0, shape[1]) for k in range(0, shape[2])]

    pool = Pool()

    logging.info("Entering main cycle...")
    res = pool.map(process_field_calc, points)

    res = numpy.array(res)
    field.p = numpy.reshape(res[:, 0], shape)
    field.vn = numpy.reshape(res[:, 1], shape)

    duration = time.clock() - start_time
    duration = timedelta(seconds=duration)
    logging.info("-----Analitical field calculation has been finished. It took %s", duration)


def load_skull_info_from_matlab():
    # import codecs
    from scipy import io as scipyio

    mat = scipyio.loadmat(r"d:\yandex_disk\SublimeProjects\SkullWay\test_files\3d_skull_coords\3d.mat")

    x = mat['x']/1000
    y = mat['y']/1000
    z3d_inside = mat['z3d_inside']/1000
    z3d_outside = mat['z3d_outside']/1000
    return [x, y, z3d_inside, z3d_outside]

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
    skull_info = point[9]

    field_x, field_y, field_z = field.get_cartesian_coords(i, j, k)
    field_normal_x, field_normal_y, field_normal_z = field.normal(i, j, k)

    rx, ry, rz, r = calc_distance_vectors(sources, field_x, field_y, field_z)
    shift_skull = propagate_skull(sources, field_x, field_y, field_z, skull_info)
    
    sines = calc_sines(sources, rx, ry, rz, r)
    coses = calc_coses(rx, ry, rz, r, field_normal_x, field_normal_y, field_normal_z)

    # check if r crosses the source
    # bad_rays_map = check_crosses(rx, ry, rz, sources)

    p, vn = sum_up_analyt(wave_number, r, sines, element_radius, sources['phase_shift'], coses / (density * speed_of_sound),shift_skull)
    return p, vn


def propagate_skull_3d(sources, field_x, field_y, field_z, skull_info):
    freq = 650*1.0e3  # Hz
    alphaDB = 22*freq*1.0e-06
    alpha_at = alphaDB*100/8.686  # 1/m

    c_skull = 4100  # m/s
    c_water = 1500  # m/s

    ro_skull = 1800
    ro_water = 1000

    z_skull = c_skull*ro_skull
    z_water = c_water*ro_water

    # print(z_skull)
    # print(z_water)

    T12 = 4*z_skull*z_water / ((z_skull+z_water)*(z_skull+z_water))
    # print(T12)

    dk = 2.0 * math.pi * freq * (1.0/c_skull - 1.0/c_water)

    return 1.0

    x = skull_info[0]
    y = skull_info[1]
    z3d_inside = skull_info[2]
    z3d_outside = skull_info[3]

    # len_of_sources = numpy.shape(sources['Xs'])[0]
    # len_of_sources = numpy.shape(sources['Ys'])[0]
    len_of_sources = numpy.shape(sources['Zs'])[0]

    bone_distance = numpy.zeros(len_of_sources)

    for i in range(0, len_of_sources):
        x0 = sources['Xs'][i]
        y0 = sources['Ys'][i]
        z0 = sources['Zs'][i]

        xf = field_x
        yf = field_y
        zf = field_z

        M0M1_outside_X = xf-x0
        M0M1_outside_Y = yf-y0
        M0M1_outside_Z = zf-z0
        M0M1_module = numpy.sqrt(M0M1_outside_X**2 + M0M1_outside_Y**2 + M0M1_outside_Z**2)

        alpha = M0M1_outside_X/M0M1_module
        betha = M0M1_outside_Y/M0M1_module
        gamma = M0M1_outside_Z/M0M1_module
        s = numpy.array([alpha, betha, gamma])  # Направляющий вектор
        # t = numpy.linspace(0, numpy.sqrt((xf-x0)**2+(yf-y0)**2+(zf-z0)**2),100)  # Параметр

        # %Уравнение прямой%
        # x_line = x0+alpha*t
        # y_line = y0+betha*t
        # z_line = z0+gamma*t

        d_outside = numpy.zeros_like(x)
        d_inside = numpy.zeros_like(x)

        # %Ищем расстояние от всех точек плоскости до прямой
        direction_vector_of_line = s

        vector_to_point_Xs = x - x0
        vector_to_point_Ys = y - y0
        vector_to_point_Zs = z3d_outside - z0
        vector_to_point = numpy.array([vector_to_point_Xs, vector_to_point_Ys, vector_to_point_Zs])

        cross_product = numpy.cross(direction_vector_of_line, vector_to_point, axisb=0, axisc=0)
        distances = numpy.sqrt(cross_product[0,:,:]**2 + cross_product[1,:,:]**2 + cross_product[2,:,:]**2)

        # Индексы минимумов
        ids = numpy.unravel_index(distances.argmin(), distances.shape)
        i_min_outside = ids[0]
        j_min_outside = ids[1]

        if i_min_outside == 0:
            i_min_outside = 1
        if j_min_outside == 0:
            j_min_outside = 1

        if i_min_outside == len(x)-1:
            i_min_outside = len(x)-2
        if j_min_outside == len(x)-1:
            j_min_outside = len(x)-2

        # %по трём точкам вокруг найденной строим плоскость
        x1_outside = x[i_min_outside-1,j_min_outside-1]
        y1_outside = y[i_min_outside-1,j_min_outside-1]
        z1_outside = z3d_outside[i_min_outside-1,j_min_outside-1]

        x2_outside = x[i_min_outside+1,j_min_outside-1]
        y2_outside = y[i_min_outside+1,j_min_outside-1]
        z2_outside = z3d_outside[i_min_outside+1,j_min_outside-1]

        x3_outside = x[i_min_outside,j_min_outside+1]
        y3_outside = y[i_min_outside,j_min_outside+1]
        z3_outside = z3d_outside[i_min_outside,j_min_outside+1]

        A_outside = (y2_outside-y1_outside)*(z3_outside-z1_outside)-(z2_outside-z1_outside)*(y3_outside-y1_outside)
        B_outside = -(x2_outside-x1_outside)*(z3_outside-z1_outside)+(z2_outside-z1_outside)*(x3_outside-x1_outside)
        C_outside = (x2_outside-x1_outside)*(y3_outside-y1_outside)-(y2_outside-y1_outside)*(x3_outside-x1_outside)
        D_outside = -A_outside*x1_outside-B_outside*y1_outside-C_outside*z1_outside

        # %ищем пересечение плоскости с прямой
        t_cross_outside = -(A_outside*x0+B_outside*y0+C_outside*z0+D_outside)/(A_outside*alpha+B_outside*betha+C_outside*gamma)

        # %Результат%
        x_cross_outside = x0+alpha*t_cross_outside
        y_cross_outside = y0+betha*t_cross_outside
        z_cross_outside = z0+gamma*t_cross_outside

        # %всё то же для внутренней поверхности

        vector_to_point_Zs = z3d_inside - z0
        vector_to_point = numpy.array([vector_to_point_Xs, vector_to_point_Ys, vector_to_point_Zs])

        # print(numpy.shape(vector_to_point))
        cross_product = numpy.cross(direction_vector_of_line, vector_to_point, axisb=0, axisc=0)
        distances = numpy.sqrt(cross_product[0,:,:]**2 + cross_product[1,:,:]**2 + cross_product[2,:,:]**2)

        # Индексы минимумов
        ids = numpy.unravel_index(distances.argmin(), distances.shape)
        i_min_inside = ids[0]
        j_min_inside = ids[1]

        if i_min_inside == 0:
            i_min_inside = 1
        if j_min_inside == 0:
            j_min_inside = 1

        if i_min_inside == len(x)-1:
            i_min_inside = len(x)-2
        if j_min_inside == len(x)-1:
            j_min_inside = len(x)-2

        x1_inside = x[i_min_inside-1,j_min_inside-1]
        y1_inside = y[i_min_inside-1,j_min_inside-1]
        z1_inside = z3d_inside[i_min_inside-1,j_min_inside-1]

        x2_inside = x[i_min_inside+1,j_min_inside-1]
        y2_inside = y[i_min_inside+1,j_min_inside-1]
        z2_inside = z3d_inside[i_min_inside+1,j_min_inside-1]

        x3_inside = x[i_min_inside,j_min_inside+1]
        y3_inside = y[i_min_inside,j_min_inside+1]
        z3_inside = z3d_inside[i_min_inside,j_min_inside+1]

        A_inside = (y2_inside-y1_inside)*(z3_inside-z1_inside)-(z2_inside-z1_inside)*(y3_inside-y1_inside)
        B_inside = -(x2_inside-x1_inside)*(z3_inside-z1_inside)+(z2_inside-z1_inside)*(x3_inside-x1_inside)
        C_inside = (x2_inside-x1_inside)*(y3_inside-y1_inside)-(y2_inside-y1_inside)*(x3_inside-x1_inside)
        D_inside = -A_inside*x1_inside-B_inside*y1_inside-C_inside*z1_inside

        t_cross_inside = -(A_inside*x0+B_inside*y0+C_inside*z0+D_inside)/(A_inside*alpha+B_inside*betha+C_inside*gamma)

        x_cross_inside = x0+alpha*t_cross_inside
        y_cross_inside = y0+betha*t_cross_inside
        z_cross_inside = z0+gamma*t_cross_inside

        bone_distance[i] = numpy.sqrt((x_cross_inside - x_cross_outside) ** 2 + (y_cross_inside - y_cross_outside) ** 2 + (z_cross_inside - z_cross_outside) ** 2)
    
    # print("Found!")
    
    # bone_distance = numpy.sqrt((x_cross_inside - x_cross_outside) ** 2 + (y_cross_inside - y_cross_outside) ** 2 + (z_cross_inside - z_cross_outside) ** 2)

    return numpy.exp(-alpha * bone_distance) * numpy.exp(1j * dk * bone_distance) * T12
    # return numpy.exp(1j * dk * bone_distance)
    # return 1.0
    # return numpy.exp(-alpha_at * bone_distance) * T12
    # return T12


def propagate_skull(sources, field_x, field_y, field_z, skull_info):
    freq = 650*1.0e3  # Hz
    alphaDB = 22*freq*1.0e-06
    alpha_at = alphaDB*100/8.686  # 1/m

    c_skull = 4100  # m/s
    c_water = 1500  # m/s

    ro_skull = 1800
    ro_water = 1000

    z_skull = c_skull*ro_skull
    z_water = c_water*ro_water

    # print(z_skull)
    # print(z_water)

    T12 = 4*z_skull*z_water / ((z_skull+z_water)*(z_skull+z_water))

    dk = 2.0 * math.pi * freq * (1.0/c_skull - 1.0/c_water)

    x = skull_info[0]
    y = skull_info[1]
    z3d_inside = skull_info[2]
    z3d_outside = skull_info[3]

    # len_of_sources = numpy.shape(sources['Xs'])[0]
    # len_of_sources = numpy.shape(sources['Ys'])[0]
    len_of_sources = numpy.shape(sources['Zs'])[0]

    bone_distance = numpy.zeros(len_of_sources)

    for i in range(0, len_of_sources):
        x0 = sources['Xs'][i]
        y0 = sources['Ys'][i]
        z0 = sources['Zs'][i]

        xf = field_x
        yf = field_y
        zf = field_z

        alpha = (xf-x0)/numpy.sqrt((xf-x0)**2+(yf-y0)**2+(zf-z0)**2)
        betha = (yf-y0)/numpy.sqrt((xf-x0)**2+(yf-y0)**2+(zf-z0)**2)
        gamma = (zf-z0)/numpy.sqrt((xf-x0)**2+(yf-y0)**2+(zf-z0)**2)
        s = numpy.array([alpha, betha, gamma])  # Направляющий вектор
        # t = numpy.linspace(0, numpy.sqrt((xf-x0)**2+(yf-y0)**2+(zf-z0)**2),100)  # Параметр

        # %Уравнение прямой%
        # x_line = x0+alpha*t
        # y_line = y0+betha*t
        # z_line = z0+gamma*t

        d_outside = numpy.zeros_like(x)
        d_inside = numpy.zeros_like(x)

        # %Ищем расстояние от всех точек плоскости до прямой

        M0M1_outside_X = xf-x0
        M0M1_outside_Y = yf-y0
        M0M1_outside_Z = zf-z0
        M0M1_module = numpy.sqrt(M0M1_outside_X**2 + M0M1_outside_Y**2 + M0M1_outside_Z**2)

        direction_vector_of_line_X = M0M1_outside_X/M0M1_module
        direction_vector_of_line_Y = M0M1_outside_Y/M0M1_module
        direction_vector_of_line_Z = M0M1_outside_Z/M0M1_module
        direction_vector_of_line = numpy.array([direction_vector_of_line_X, direction_vector_of_line_Y, direction_vector_of_line_Z])

        vector_to_point_Xs = x - x0
        vector_to_point_Ys = y - y0
        vector_to_point_Zs = z3d_outside - z0
        vector_to_point = numpy.array([vector_to_point_Xs, vector_to_point_Ys, vector_to_point_Zs])

        cross_product = numpy.cross(direction_vector_of_line, vector_to_point, axisb=0, axisc=0)

        distances = numpy.sqrt(cross_product[0,:,:]**2 + cross_product[1,:,:]**2 + cross_product[2,:,:]**2)
        # distances = numpy.sqrt(numpy.sum(cross_product**2, axis=0))

        # Индексы минимумов
        ids = numpy.unravel_index(distances.argmin(), distances.shape)
        i_min_outside = ids[0]
        j_min_outside = ids[1]

        if i_min_outside == 0:
            i_min_outside = 1
        if j_min_outside == 0:
            j_min_outside = 1

        if i_min_outside == len(x)-1:
            i_min_outside = len(x)-2
        if j_min_outside == len(x)-1:
            j_min_outside = len(x)-2

        # %по трём точкам вокруг найденной строим плоскость
        x1_outside = x[i_min_outside-1,j_min_outside-1]
        y1_outside = y[i_min_outside-1,j_min_outside-1]
        z1_outside = z3d_outside[i_min_outside-1,j_min_outside-1]

        x2_outside = x[i_min_outside+1,j_min_outside-1]
        y2_outside = y[i_min_outside+1,j_min_outside-1]
        z2_outside = z3d_outside[i_min_outside+1,j_min_outside-1]

        x3_outside = x[i_min_outside,j_min_outside+1]
        y3_outside = y[i_min_outside,j_min_outside+1]
        z3_outside = z3d_outside[i_min_outside,j_min_outside+1]

        A_outside = (y2_outside-y1_outside)*(z3_outside-z1_outside)-(z2_outside-z1_outside)*(y3_outside-y1_outside)
        B_outside = -(x2_outside-x1_outside)*(z3_outside-z1_outside)+(z2_outside-z1_outside)*(x3_outside-x1_outside)
        C_outside = (x2_outside-x1_outside)*(y3_outside-y1_outside)-(y2_outside-y1_outside)*(x3_outside-x1_outside)
        D_outside = -A_outside*x1_outside-B_outside*y1_outside-C_outside*z1_outside

        # %ищем пересечение плоскости с прямой
        t_cross_outside = -(A_outside*x0+B_outside*y0+C_outside*z0+D_outside)/(A_outside*alpha+B_outside*betha+C_outside*gamma)

        # %Результат%
        x_cross_outside = x0+alpha*t_cross_outside
        y_cross_outside = y0+betha*t_cross_outside
        z_cross_outside = z0+gamma*t_cross_outside

        # %всё то же для внутренней поверхности

        vector_to_point_Zs = z3d_inside - z0
        vector_to_point = numpy.array([vector_to_point_Xs, vector_to_point_Ys, vector_to_point_Zs])

        # print(numpy.shape(vector_to_point))
        cross_product = numpy.cross(direction_vector_of_line, vector_to_point, axisb=0, axisc=0)
        distances = numpy.sqrt(cross_product[0,:,:]**2 + cross_product[1,:,:]**2 + cross_product[2,:,:]**2)

        # Индексы минимумов
        ids = numpy.unravel_index(distances.argmin(), distances.shape)
        i_min_inside = ids[0]
        j_min_inside = ids[1]

        if i_min_inside == 0:
            i_min_inside = 1
        if j_min_inside == 0:
            j_min_inside = 1

        if i_min_inside == len(x)-1:
            i_min_inside = len(x)-2
        if j_min_inside == len(x)-1:
            j_min_inside = len(x)-2

        x1_inside = x[i_min_inside-1,j_min_inside-1]
        y1_inside = y[i_min_inside-1,j_min_inside-1]
        z1_inside = z3d_inside[i_min_inside-1,j_min_inside-1]

        x2_inside = x[i_min_inside+1,j_min_inside-1]
        y2_inside = y[i_min_inside+1,j_min_inside-1]
        z2_inside = z3d_inside[i_min_inside+1,j_min_inside-1]

        x3_inside = x[i_min_inside,j_min_inside+1]
        y3_inside = y[i_min_inside,j_min_inside+1]
        z3_inside = z3d_inside[i_min_inside,j_min_inside+1]

        A_inside = (y2_inside-y1_inside)*(z3_inside-z1_inside)-(z2_inside-z1_inside)*(y3_inside-y1_inside)
        B_inside = -(x2_inside-x1_inside)*(z3_inside-z1_inside)+(z2_inside-z1_inside)*(x3_inside-x1_inside)
        C_inside = (x2_inside-x1_inside)*(y3_inside-y1_inside)-(y2_inside-y1_inside)*(x3_inside-x1_inside)
        D_inside = -A_inside*x1_inside-B_inside*y1_inside-C_inside*z1_inside

        t_cross_inside = -(A_inside*x0+B_inside*y0+C_inside*z0+D_inside)/(A_inside*alpha+B_inside*betha+C_inside*gamma)

        x_cross_inside = x0+alpha*t_cross_inside
        y_cross_inside = y0+betha*t_cross_inside
        z_cross_inside = z0+gamma*t_cross_inside

        bone_distance[i] = numpy.sqrt((x_cross_inside - x_cross_outside) ** 2 + (y_cross_inside - y_cross_outside) ** 2 + (z_cross_inside - z_cross_outside) ** 2)
    
    # print("Found!")
    
    # bone_distance = numpy.sqrt((x_cross_inside - x_cross_outside) ** 2 + (y_cross_inside - y_cross_outside) ** 2 + (z_cross_inside - z_cross_outside) ** 2)

    # return numpy.exp(-alpha * bone_distance) * numpy.exp(1j * dk * bone_distance) * T12
    # return numpy.exp(1j * dk * bone_distance)
    return 1.0
    # return numpy.exp(-alpha_at * bone_distance) * T12


def propagate_skull_2D(sources, field_x, field_y, field_z, skull_info):
    freq = 650*1.0e3  # Hz
    alphaDB = 22*freq*1.0e-06
    alpha = alphaDB*100/8.686  # 1/m

    c_skull = 4100  # m/s
    c_water = 1500  # m/s

    ro_skull = 1800
    ro_water = 1000

    z_skull = c_skull*ro_skull
    z_water = c_water*ro_water

    # print(z_skull)
    # print(z_water)

    T12 = 4*z_skull*z_water / ((z_skull+z_water)*(z_skull+z_water))

    dk = 2.0 * math.pi * freq * (1.0/c_skull - 1.0/c_water)

    import codecs
    # Adding coordinates
    # 'utf-8-sig' means that it skips BOM in UTF files created in Notepad
    file_path1 = r"d:\yandex_disk\SublimeProjects\SkullWay\test_files\x.txt"
    file_path2 = r"d:\yandex_disk\SublimeProjects\SkullWay\test_files\z_inside.txt"
    file_path3 = r"d:\yandex_disk\SublimeProjects\SkullWay\test_files\z_outside.txt"
    with codecs.open(file_path1, encoding='utf-8-sig') as f:
        x = [float(line.strip())/1000 for line in f]

    with codecs.open(file_path2, encoding='utf-8-sig') as f:
        z_inside = [float(line.strip())/1000 for line in f]

    with codecs.open(file_path3, encoding='utf-8-sig') as f:
        z_outside = [float(line.strip())/1000 for line in f]

    DX = field_x - sources['Xs']
    DZ = field_z - sources['Zs']

    def_i = numpy.where(DX == 0.0)

    k = DZ/DX

    k[def_i[0]] = 100

    b = field_z - k*field_x
    z_ray = k*x + b
    dz_out = z_ray - z_outside
    dz_in = z_ray - z_inside

    d_out = numpy.nonzero(numpy.diff(numpy.sign(dz_out)))
    d_in = numpy.nonzero(numpy.diff(numpy.sign(dz_in)))

    d_out = d_out[0]
    d_in = d_in[0]
    # print(d_out)
    # print(d_in)

    if (d_out.size > 1) or (d_in.size > 1) or (d_out.size == 0) or (d_in.size == 0):
        # print('Out of brain')
        return numpy.ones_like(sources['Xs'])
    else:
        d_out = d_out[0]
        d_in = d_in[0]
        # print(d_out)
        # print(d_in)
        z1_out = z_outside[d_out]
        z2_out = z_outside[d_out+1]
        x1_out = x[d_out]
        x2_out = x[d_out+1]
        x_cross_out = (z1_out-(z2_out-z1_out)/(x2_out-x1_out)*x1_out)/(k-(z2_out-z1_out)/(x2_out-x1_out))
        z_cross_out = k*x_cross_out

        z1_in = z_inside[d_in]
        z2_in = z_inside[d_in+1]
        x1_in = x[d_in]
        x2_in = x[d_in+1]
        x_cross_in = (z1_in-(z2_in-z1_in)/(x2_in-x1_in)*x1_in)/(k-(z2_in-z1_in)/(x2_in-x1_in))
        z_cross_in = k*x_cross_in
    
    bone_distance = numpy.sqrt((x_cross_in - x_cross_out) ** 2 + (z_cross_in - z_cross_out) ** 2)

    return numpy.exp(-alpha * bone_distance) * numpy.exp(1j * dk * bone_distance) * T12
    # return numpy.exp(1j * dk * bone_distance)
    # return numpy.exp(-alpha * bone_distance) * T12

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


def sum_up_analyt(wave_number, r, sines, element_radius, phase_shift, coef_vn, shift_skull):
    a = element_radius
    val = numpy.exp(1j * wave_number * r) * (1 / (2.0 * r))
    val *= shift_skull
    bessel_val = 2 * j1(wave_number * a * sines) / (wave_number * a * sines + 1.0 * (sines == 0.0)) + 1.0 * (sines == 0.0)
    val *= bessel_val * (-1j * wave_number * element_radius * element_radius)
    p = val * numpy.exp(-1j * wave_number * phase_shift)
    vn = p * coef_vn
    return numpy.sum(p), numpy.sum(vn)
