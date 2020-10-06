# -*- coding: utf-8 -*-
"""Module with functions for Numerical calculation of Rayleigh Integral
"""
import math
import numpy
import time
import logging
from multiprocessing import Pool, Manager
from datetime import timedelta
import ribs
import pickle
import scipy


class RayCalc:

    """Основной расчетный класс
    """

    def __init__(self, field, trans, medium, num_of_tris_on_line):
        """Создание задания на расчет
        """
        self.field = field
        self.trans = trans
        self.medium = medium
        self.num_of_tris_on_line = num_of_tris_on_line

    def doit(self):
        """Запуск расчета
        """
        self.start_timer()
        self.calc_params()
        self.form_sources()
        self.multiproc_calc_status()
        self.stop_timer()

    def start_timer(self):
        logging.info("-----Numerical field calculation has been started-----")
        logging.info("Using transducer '%s' and medium '%s'", self.trans.name, self.medium.name)
        self.start_time = time.clock()

    def calc_params(self):
        """Вычисление длины волны и волнового числа
        """
        self.wave_length = self.medium.speed_of_sound / self.trans.frequency
        self.wave_number = 2.0 * math.pi / self.wave_length
        logging.info("Wave length = %f mm", 1.0e3 * self.wave_length)
        logging.info("Wave number = %f m-1", self.wave_number)

    def form_sources(self):
        logging.info("Forming sources...")
        self.tris = slice_trans_plane_on_tris(self.trans.aperture, self.num_of_tris_on_line)
        project_tris_on_trans(self.tris, self.trans.curvature_radius)
        compute_centroids_of_tris(self.tris)
        calc_tris_normals_to_trans_focus(self.tris, self.trans.curvature_radius)
        self.delete_extra_tris()
        calc_tris_areas(self.tris)
        # import draw_transducer
        # draw_transducer.draw_transducer_tris(self.trans, self.tris)
        save_tris_on_disk(self.tris)
        self.sources = form_sources_from_tris(self.tris)
        logging.info("Slicing of transducer into triangles has been finished. Number of active triangles: %d. Side of triangle is about %f mm", len(self.tris), math.sqrt(self.tris[0].area) * 1.0e03)
        self.sources_power_check()

    def multiproc_calc_status(self):
        "Попытка сделать статус асинхр задачи"
        STATUS_CALC_EVERY_S = 10
        self.q = Manager().Queue()
        shape = numpy.shape(self.field.p)
        points = [(i, j, k) for i in range(0, shape[0]) for j in range(0, shape[1]) for k in range(0, shape[2])]
        num_points = len(points)
        pool = Pool()

        logging.info("Entering main cycle...")

        res_async = pool.starmap_async(self.process_field_calc_status, points)
        time0 = time.clock()
        while 1:
            if (res_async.ready()):
                break
            duration = time.clock() - time0
            percent = self.q.qsize()/num_points
            if percent != 0.0:
                logging.info('Done {0:%}, remaining time: {1}'.format((percent), timedelta(seconds=duration/percent-duration)))
            time.sleep(STATUS_CALC_EVERY_S)
        res = res_async.get()
        res = numpy.array(res)
        self.field.p = numpy.reshape(res[:, 0], shape)
        self.field.vn = numpy.reshape(res[:, 1], shape)
        self.field.p *= 1j * self.wave_number / (2.0 * numpy.pi)
        self.field.vn *= 1 / (2.0 * numpy.pi * self.medium.density * self.medium.speed_of_sound)

    def process_field_calc_status(self, i, j, k):
        field_x, field_y, field_z = self.field.get_cartesian_coords(i, j, k)
        field_normal_z = 1.0
        # field_normal_x, field_normal_y, field_normal_z = self.field.normal(i, j, k)

        rx, ry, rz, r = calc_distance_vectors(self.sources, field_x, field_y, field_z)
        coses = calc_coses(rz, r, field_normal_z)
        # coses = calc_coses(rx, ry, rz, r, field_normal_x, field_normal_y, field_normal_z)

        # # check if r crosses the source
        # bad_rays_map = check_crosses(rx, ry, rz, self.sources)
        # p, vn = sum_up_sources_rayleigh(self.wave_number, r, self.sources['Ss'] * bad_rays_map, 1.0, coses)

        # without check if r crosses the source
        p, vn = sum_up_sources_rayleigh(self.wave_number, r, self.sources['Ss'], coses)
        self.q.put(None)
        return p, vn

    def stop_timer(self):
        duration = time.clock() - self.start_time
        duration = timedelta(seconds=duration)
        logging.info("-----Numerical field calculation has been finished. It took %s", duration)

    def sources_power_check(self):
        # Power checking
        sources_power = calc_sources_power(self.sources, self.medium)
        ideal_power = calc_ideal_power(self.trans, self.medium)
        logging.info("Sources power = " + str(sources_power))
        logging.info("Ideal power = " + str(ideal_power))

    def delete_extra_tris(self):
        self.tris = [tri for tri in self.tris if tri.is_need(self.trans)]


class RayCalcFromField(RayCalc):

    def __init__(self, field, trans, medium, field_init):
        self.field = field
        self.trans = trans
        self.medium = medium
        self.field_init = field_init

    def form_sources(self):
        logging.info("Forming sources...")
        self.sources = form_sources_from_field(self.field_init)
        # logging.info("Slicing of transducer into triangles has been finished. Number of active triangles: %d. Side of triangle is about %f mm", len(self.tris), math.sqrt(self.tris[0].area) * 1.0e03)
        # self.sources_power_check()


class RayCalcIdealFull(RayCalc):

    def delete_extra_tris(self):
        self.tris = [tri for tri in self.tris if tri.is_need_full_trans_full(self.trans)]


class RayCalcIdealRibs(RayCalc):

    def __init__(self, field, trans, medium, num_of_tris_on_line, ribs_phantom):
        self.field = field
        self.trans = trans
        self.medium = medium
        self.num_of_tris_on_line = num_of_tris_on_line
        self.ribs_phantom = ribs_phantom

    def delete_extra_tris(self):
        self.tris = [tri for tri in self.tris if tri.is_need_full_trans_ribs(self.trans, self.ribs_phantom)]

    # def delete_extra_tris(self):
    #   import random
    #   # self.tris = [tri for tri in self.tris if tri.is_need_full_trans_ribs(self.trans, self.ribs_phantom)]
    #   tris_f = [tri for tri in self.tris if tri.is_need_full_trans_ribs(self.trans, self.ribs_phantom)]
    #   random.seed()
    #   self.tris = random.sample(tris_f, int(len(tris_f)*0.005))

    # def delete_extra_tris(self):
    #   # import random
    #   # self.tris = [tri for tri in self.tris if tri.is_need_full_trans_ribs(self.trans, self.ribs_phantom)]
    #   tris_f = [tri for tri in self.tris if tri.is_need_full_trans_ribs(self.trans, self.ribs_phantom)]
    #   # random.seed()
    #   # self.tris = random.sample(tris_f, len(tris_f)//2)
    #   # for i in range(0, len(tris_f)-1):
    #       # if i % 2 == 0:
    #           # tris_f.pop(i)
    #   self.tris = tris_f[:2500]

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
        self.nx = 0.0
        self.ny = 0.0
        self.nz = 1.0

    def compute_centroids(self):
        self.xc = (self.x1 + self.x2 + self.x3) / 3.0
        self.yc = (self.y1 + self.y2 + self.y3) / 3.0
        self.zc = (self.z1 + self.z2 + self.z3) / 3.0

    def project_on_trans(self, F):
        """Calculates z coordinate of trinangle vertexes
        """

        x0 = self.x1
        y0 = self.y1
        self.z1 = F - math.sqrt(F * F - x0 * x0 - y0 * y0)

        x0 = self.x2
        y0 = self.y2
        self.z2 = F - math.sqrt(F * F - x0 * x0 - y0 * y0)

        x0 = self.x3
        y0 = self.y3
        self.z3 = F - math.sqrt(F * F - x0 * x0 - y0 * y0)

    def calc_normals(self, F):
        # Normals to the center of curvature of trans = (0, 0, F)
        self.nx = (0.0 - self.xc) / F
        self.ny = (0.0 - self.yc) / F
        self.nz = (F - self.zc) / F

    def compute_area(self):
        """Computes area of triangle according to its vertexes
        """

        """Косое произведение векторов
        A = (x2-x1; y2-y1; z2-z1)
        B = (x3-x1; y3-y1; z3-z1)
        S = 0.5*sqrt((Ay*Bz - Az*By)^2 + (Az*Bx - Ax*Bz)^2 + (Ax*By - Ay*Bx)^2 )
        """
        a_x = self.x2 - self.x1
        a_y = self.y2 - self.y1
        a_z = self.z2 - self.z1

        b_x = self.x3 - self.x1
        b_y = self.y3 - self.y1
        b_z = self.z3 - self.z1

        self.area = 0.5 * math.sqrt((a_y * b_z - a_z * b_y) ** 2 + (a_z * b_x - a_x * b_z) ** 2 + (a_x * b_y - a_y * b_x) ** 2)

        """По теореме Герона"""
        # a = math.sqrt((self.x1-self.x2)**2 + (self.y1-self.y2)**2 + (self.z1-self.z2)**2)
        # b = math.sqrt((self.x1-self.x3)**2 + (self.y1-self.y3)**2 + (self.z1-self.z3)**2)
        # c = math.sqrt((self.x2-self.x3)**2 + (self.y2-self.y3)**2 + (self.z2-self.z3)**2)
        # p = 0.5 * (a + b + c)
        # self.area = math.sqrt(p * (p - a) * (p - b) * (p - c))
        

    def is_need(self, trans):
        # Для скорости
        R2 = trans.element_radius * trans.element_radius
        for el in trans.elements:
            # Go to local coordinate system of each element of transducer
            x0 = self.xc - el['center_x']
            y0 = self.yc - el['center_y']

            # Расстояние от центра треугольника до центра элемента в квадрате
            # для скорости
            hypot2 = x0 * x0 + y0 * y0

            # Сравниваем Расстояние от центра треугольника до центра элемента
            # просто с радиусом элемента
            if hypot2 < R2:
                cos_gamma = math.cos(
                    math.asin(math.hypot(el['center_x'], el['center_y']) / trans.curvature_radius))
                a2 = R2
                b2 = cos_gamma * cos_gamma * R2
 
                # Сравниваем Расстояние от центра треугольника до центра
                # элемента с наименьшим радиусом эллипса (проекция повернутого
                # элемента на плоскость)
                if hypot2 < b2:
                    return True
                else:
                    phi = math.atan2(el['center_y'], el['center_x'])
                    beta = math.atan2(y0, x0)
                    rad2 = a2 * \
                        math.sin(phi - beta) * math.sin(phi - beta) + \
                        b2 * math.cos(phi - beta) * math.cos(phi - beta)

                    # Сравниваем Расстояние от центра треугольника до центра
                    # элемента с остальной частью эллипса (проекция повернутого
                    # элемента на плоскость)
                    if hypot2 < rad2:
                        return True
        return False

    def is_need_full_trans_ribs(self, trans, ribs_phantom):
        """
        Для разбития на треугольники всей поверхности решетки с учетом ребер
        """
        # Precalc
        R2 = trans.aperture * trans.aperture / 4.0
        # H2 = trans.hole_radius*trans.hole_radius

        if ribs.is_point_trans_intersects_ribs(self.xc, self.yc, self.zc, trans.curvature_radius, ribs_phantom['dist_to_ribs_plane'], ribs_phantom):
            return False

        if self.xc ** 2 + self.yc ** 2 > R2:
            return False
                # if self.x1**2 + self.y1**2 >= H2:
                #   if self.x2**2 + self.y2**2 >= H2:
                #       if self.x3**2 + self.y3**2 >= H2:
                #           return True
        return True

    def is_need_full_trans_full(self, trans):
        """
        Для разбития на треугольники всей поверхности решетки
        """
        # Precalc
        R2 = trans.aperture * trans.aperture / 4.0

        if self.xc ** 2 + self.yc ** 2 >= R2:
            return False

        return True


def slice_trans_plane_on_tris(aperture, num_of_tris_on_line):
    """Slice the transducer aperture plane into triangles.
    Each rectangle of calc grid is dividing on 2 triangles by the diagonal.
    """
    logging.info("Slicing of transducer aperture plane into triangles has been started. Number of triangles on line %d", num_of_tris_on_line)

    num_of_rects_on_line = num_of_tris_on_line // 2
    dx = aperture / num_of_rects_on_line
    dy = -dx
    a = aperture / 2.0
    tris = []

    for j in range(0, num_of_rects_on_line):
        for i in range(0, num_of_rects_on_line):
            tri1 = Triangle(-a + i * dx, a + j * dy, -a + i * dx, a + (j + 1) * dy, -a + (i + 1) * dx, a + j * dy)
            # tri2 = Triangle(-a + (i + 1) * dx, a + j * dy, -a + i * dx, a + (j + 1) * dy, -a + (i + 1) * dx, a + (j + 1) * dy)
            tris.append(tri1)
            # tris.append(tri2)

    return tris


def compute_centroids_of_tris(tris):
    for tri in tris:
        tri.compute_centroids()


def project_tris_on_trans(tris, F):
    for tri in tris:
        tri.project_on_trans(F)


def calc_tris_normals_to_trans_focus(tris, F):
    for tri in tris:
        tri.calc_normals(F)


def calc_tris_areas(tris):
    for tri in tris:
        tri.compute_area()


def process_field_calc(field, sources, wave_number, i, j, k):
    field_x, field_y, field_z = field.get_cartesian_coords(i, j, k)
    field_normal_z = 1.0
    # field_normal_x, field_normal_y, field_normal_z = field.normal(i, j, k)

    rx, ry, rz, r = calc_distance_vectors(sources, field_x, field_y, field_z)
    coses = calc_coses(rz, r, field_normal_z)
    # coses = calc_coses(rx, ry, rz, r, field_normal_x, field_normal_y, field_normal_z)

    # # check if r crosses the source
    # bad_rays_map = check_crosses(rx, ry, rz, sources)
    # p, vn = sum_up_sources_rayleigh(wave_number, r, sources['Ss'] * bad_rays_map, 1.0, coses)

    # without check if r crosses the source
    p, vn = sum_up_sources_rayleigh(wave_number, r, sources['Ss'], coses)

    return p, vn


def calc_power(field):
    """Calc power as a sum of Poynting vectors (P = 0.5*Re(p * conj(v))) in each point and multiplyies it with area"""

    poynt_in_points = 0.5*numpy.real(field.p * numpy.conj(field.vn))
    power = numpy.sum(poynt_in_points)
    power *= field.one_pixel_area
    return power


def calc_power_in_main_focus(field, focus_size):
    # 1. Вычленим окно с основным фокусом, используя focus_size
    # 2. Просуммируем в окне все, что больше уровня -6 дБ
    """Calc power as a sum of Poynting vectors (P = 0.5*Re(p * conj(v))) in each point and multiplyies it with area"""
    poynt_in_points = 0.5*numpy.real(field.p * numpy.conj(field.vn))

    # Maximum
    max_focus = poynt_in_points.max()
    print("max_focus = ", max_focus)
    print("max_focus p/p0= ", numpy.absolute(field.p).max())
    max_indices = numpy.argwhere(max_focus == poynt_in_points)
    max_index_x = max_indices[0][0]
    max_index_y = max_indices[0][1]
    print("max_indices = ", max_index_x, max_index_y)
    focus_size_index = focus_size // field.dx
    print("focus_size_index = ", focus_size_index)
    field_window = field
    field_window.p = field.p[(max_index_x-focus_size_index):(max_index_x+focus_size_index),(max_index_y-focus_size_index):(max_index_y+focus_size_index),:]
    field_window.vn = field.vn[(max_index_x-focus_size_index):(max_index_x+focus_size_index),(max_index_y-focus_size_index):(max_index_y+focus_size_index),:]
    field_window.x = field.x[(max_index_x-focus_size_index):(max_index_x+focus_size_index)]
    field_window.y = field.y[(max_index_y-focus_size_index):(max_index_y+focus_size_index)]

    import draw_plane_field
    draw_plane_field.draw_XY(field_window)

    poynt_in_points_window = 0.5*numpy.real(field_window.p * numpy.conj(field_window.vn))

    # print(field_window.p)
    # power_array = numpy.where(poynt_in_points_window>=0.5*max_focus)
    # print(power_array)
    # power = numpy.sum(power_array)
    # power *= field.one_pixel_area
    # print(power)
    # return power

    w = 0.0
    max_focus = poynt_in_points_window.max()
    for j, y in enumerate(field_window.y):
        for i, x in enumerate(field_window.x):
            if poynt_in_points_window[i,j,0] >= 0.25*max_focus:
                w += poynt_in_points_window[i,j,0]

    w *= field_window.one_pixel_area
    print("Power in Focus -6dB = ", w)
            # w = field.p_amp(i, j, 0) * field.p_amp(i, j, 0)
            
    # powerOnPlane *= field.one_pixel_area
    # # powerOnRibs *= field.one_pixel_area

    # print("PowerOnPlane = " + str(powerOnPlane))


def form_sources_from_tris(tris):
    sources_Xs = numpy.array([tri.xc for tri in tris])
    sources_Ys = numpy.array([tri.yc for tri in tris])
    sources_Zs = numpy.array([tri.zc for tri in tris])
    sources_Ss = numpy.array([tri.area for tri in tris])

    sources_normal_Xs = numpy.array([tri.nx for tri in tris])
    sources_normal_Ys = numpy.array([tri.ny for tri in tris])
    sources_normal_Zs = numpy.array([tri.nz for tri in tris])

    sources = {
        'Xs': sources_Xs,
        'Ys': sources_Ys,
        'Zs': sources_Zs,
        'Ss': sources_Ss,
        'n_Xs': sources_normal_Xs,
        'n_Ys': sources_normal_Ys,
        'n_Zs': sources_normal_Zs}
    return sources


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
                # sources_Ss.append(field.p[i, j, k])
                sources_Ss.append(field.vn[i, j, k]*1500*1000)
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

    sources = {
        'Xs': sources_Xs,
        'Ys': sources_Ys,
        'Zs': sources_Zs,
        'Ss': sources_Ss,
        'n_Xs': sources_normal_Xs,
        'n_Ys': sources_normal_Ys,
        'n_Zs': sources_normal_Zs}
    return sources


def save_tris_on_disk(tris, file_name="tris.bin"):
    logging.info("Saving tris on disk...")
    with open(file_name, 'wb') as f:
        pickle.dump(tris, f)


def restore_tris_from_disk(file_name="tris.bin"):
    logging.info("Restoring tris...")
    with open(file_name, 'rb') as f:
        tris = pickle.load(f)
    return tris


def calc_distance_vectors(sources, field_x, field_y, field_z):
    DX = field_x - sources['Xs']
    DY = field_y - sources['Ys']
    DZ = field_z - sources['Zs']
    return DX, DY, DZ, numpy.sqrt(DX*DX + DY*DY + DZ*DZ)


# def calc_coses(rx, ry, rz, r, normal_x, normal_y, normal_z):
#   return numpy.abs(rx * normal_x + ry * normal_y + rz * normal_z) / r


def calc_coses(rz, r, normal_z):
    return numpy.abs(rz * normal_z) / r


def check_crosses(rx, ry, rz, sources):
    # Checks if r crosses the source
    # Checking Scalar multiplying normal vector on distance vector
    # if (mult >= 0) => (angle <= 90) => OK, unit value, else (angle > 90) =>
    # Zero value
    mult = rx * sources['n_Xs'] + ry * sources['n_Ys'] + rz * sources['n_Zs']
    return 1.0 * (mult >= 0)


def sum_up_sources_rayleigh(k, r, S, cos_coef):
    recip_r = numpy.reciprocal(r)
    i_k = 1j * k
    p = scipy.exp(i_k * r) * (S * recip_r)
    vn = p * cos_coef * (i_k - recip_r)
    return numpy.sum(p), numpy.sum(vn)
    # reciprocal(x) = 1/x

# def sum_up_sources_rayleigh(k, r, S, cos_coef):
#     import numexpr as ne
#     # recip_r = numpy.reciprocal(r)
#     i_k = 1j * k
#     p = ne.evaluate("exp(i_k * r) * (S * 1/r)")
#     vn = ne.evaluate("p * cos_coef * (1j * k - 1/r)")
#     return numpy.sum(p), numpy.sum(vn)
#     # reciprocal(x) = 1/x

def reference_anal_calc_on_axis(field, trans, medium):
    """
    Расчет для проверки. Поле сферического источника на его оси

    p(z) = 2*p0/(1 - z/F) * sin(k*(z - Rmax)/2)

    Rmax = F*sqrt(1+(1-z/F)^2-2*(1-z/F)*cos(alpha))
    sin(alpha) = a/(2F)
    a - aperture of trans
    F - focus of trans
    k - wave number
    p0 - initial pressure on trans


    vn(z) = p0/(ikc0*rho) * F/(F-z)*[(exp(ikz) - exp(ikRmax))/(F-z) + ik * (exp(ikz) - exp(ikRmax)*(z-delta)/Rmax)]

    Rmax = sqrt((z-delta)^2 + (a/2)^2)
    delta = F * (1 - sqrt(1-(a/2F)^2))
    a - aperture of trans
    F - focus of trans
    k - wave number
    c0 - speed of sound
    rho - density
    p0 - initial pressure on trans

    """
    logging.info(
        "-----Reference field calculation on axis has been started-----")
    logging.info(
        "Using transducer '%s' and medium '%s'", trans.name, medium.name)
    start_time = time.clock()

    wave_length = medium.speed_of_sound / trans.frequency
    wave_number = 2 * math.pi / wave_length
    logging.info("Wave length = %f mm", 1.0e3 * wave_length)
    logging.info("Wave number = %f m-1", wave_number)

    a = trans.aperture
    F = trans.curvature_radius
    alpha = math.asin(a / (2 * F))
    cos_alpha = math.cos(alpha)
    k = wave_number
    rho = medium.density
    c0 = medium.speed_of_sound
    delta = F * (1 - math.sqrt(1 - (a / 2 * F) ** 2))
    p0 = 1

    z = field.z

    Rmax = F * \
        numpy.sqrt(1 + (1 - z / F) * (1 - z / F) - 2 * (1 - z / F) * cos_alpha)

    p = 2 * p0 / (1 - z / F) * numpy.sin(k * (z - Rmax) / 2)
    vn = p0 / (1j * k * c0 * rho) * F / (F - z) * ((numpy.exp(1j * k * z) - numpy.exp(1j * k * Rmax)) / (F - z) + 1j * k * (numpy.exp(1j * k * z) - numpy.exp(1j * k * Rmax) * (z - delta) / Rmax))

    field.p[0, 0, :] = p
    field.vn[0, 0, :] = vn

    duration = time.clock() - start_time
    duration = timedelta(seconds=duration)
    logging.info(
        "-----Reference field calculation on axis has been finished. It took %s", duration)


def reference_anal_calc_in_focal_plane(field, trans, medium):
    """
    Расчет для проверки. Поле сферического источника в фокальной плоскости

    p(r) = p0*k*delta*(2J1(sin(alpha)*k*r))/(sin(alpha)*k*r)

    sin(alpha) = a/(2F)
    a - aperture of trans
    F - focus of trans
    k - wave number
    p0 - initial pressure on trans
    delta = F - sqrt(F*F - a*a/4)

    """

    from scipy.special import j1, jn

    logging.info(
        "-----Reference field calculation in focal plane has been started-----")
    logging.info(
        "Using transducer '%s' and medium '%s'", trans.name, medium.name)
    start_time = time.clock()

    wave_length = medium.speed_of_sound / trans.frequency
    wave_number = 2 * math.pi / wave_length
    logging.info("Wave length = %f mm", 1.0e3 * wave_length)
    logging.info("Wave number = %f m-1", wave_number)

    a = trans.aperture
    F = trans.curvature_radius
    sin_alpha = a / (2 * F)
    delta = F * (1 - math.sqrt(1 - sin_alpha * sin_alpha))
    k = wave_number
    p0 = 1

    r = field.x

    # Обычное аналитическое решение
    p = p0*k*delta*(2*j1(sin_alpha*k*r))/(sin_alpha*k*r)
    
    # Дополнительная точность аналитического решения с помощью добавления бесселей
    # p = p0 * k * delta * (2) / (sin_alpha * k * r) * (j1(sin_alpha * k * r) - ((delta / a) ** 2) * jn(3, sin_alpha * k * r) + ((delta / a) ** 4) * jn(5, sin_alpha * k * r))
    p[numpy.where(r == 0)] = p0*k*delta

    field.p[:, 0, 0] = p
    duration = time.clock() - start_time
    duration = timedelta(seconds=duration)
    logging.info(
        "-----Reference field calculation in focal plane has been finished. It took %s", duration)


def calc_sources_power(sources, medium):
    """
    Calc Initial Energy on Transducer from sources.
    W = I_0 * sigma(S_i)
    S_i - area of each triangle
    I_0 - initial intensity
    I_0 = p_0^2/(2*rho*c)
    p_0 = 1
    """
    I_0 = 1.0 / (2.0 * medium.density * medium.speed_of_sound)
    S = numpy.sum(sources['Ss'])
    W = I_0 * S
    return W


def calc_ideal_power(trans, medium):
    """
    Calc Initial Energy on the whole Transducer, if it is ideal and spherical.
    W = I_0 * S
    S - area of transducer
    I_0 - initial intensity
    I_0 = p_0^2/(2*rho*c)
    p_0 = 1
    """
    I_0 = 1.0 / (2.0 * medium.density * medium.speed_of_sound)

    # Area of transducer
    # s = 2*pi*F^2*(1-cos(alpha))
    # sin(alpha) = aperture/(2*F)
    cos_alpha = math.sqrt(
        1 - (trans.aperture / (2.0 * trans.curvature_radius)) ** 2)
    S = 2 * math.pi * trans.curvature_radius * \
        trans.curvature_radius * (1 - cos_alpha)
    W = I_0 * S
    return W

if __name__ == '__main__':
    pass
    # profile_main_cycle()
    # test_slicing()
    # test case
    # trans = transducer.Transducer("Probe1", 1, 4.0e-03, 5.0e-03, 6.0e06)
    # tris = sliceTransOnTrisProper(trans, 100)

    # trans = transducer.Transducer(
    #   name = "Probe1",
    #   num_of_elements = 1,
    #   element_width = 4.0e-03,
    #   element_length = 5.0e-03,
    # frequency = 6.0e06
    # frequency = 12.0e06
    # frequency = 4.0e06
    #   frequency = 11.6e06
    #   )

    # water = medium.Medium(
    #   name = "Water",
    #   density = 1000.0,
    #   speed_of_sound = 1488.0
    #   )
