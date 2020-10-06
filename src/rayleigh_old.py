# -*- coding: utf-8 -*-
"""
Module with functions for Numerical calculation of Rayleigh Integral
Includes cases with Cartesian and cylindrical coordinates

"""
import math
import numpy
import scipy
import time
import logging
from multiprocessing import Pool
from datetime import timedelta
import ribs
import pickle
import transducer


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

	def project_on_trans(self, trans):
		"""Calculates z coordinate of trinangle vertex
		"""
		F = trans.curvature_radius

		x0 = self.x1
		y0 = self.y1
		self.z1 = F - math.sqrt(F * F - x0 * x0 - y0 * y0)

		x0 = self.x2
		y0 = self.y2
		self.z2 = F - math.sqrt(F * F - x0 * x0 - y0 * y0)

		x0 = self.x3
		y0 = self.y3
		self.z3 = F - math.sqrt(F * F - x0 * x0 - y0 * y0)

	def calc_normals(self, trans):
		# Normals to the center of curvature of trans = (0, 0, curvature_radius)
		F = trans.curvature_radius
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

		"""По теореме Герона
		a = math.sqrt((self.x1-self.x2)**2 + (self.y1-self.y2)**2 + (self.z1-self.z2)**2)
		b = math.sqrt((self.x1-self.x3)**2 + (self.y1-self.y3)**2 + (self.z1-self.z3)**2)
		c = math.sqrt((self.x2-self.x3)**2 + (self.y2-self.y3)**2 + (self.z2-self.z3)**2)
		p = 0.5 * (a + b + c)
		self.area = math.sqrt(p * (p - a) * (p - b) * (p - c))
		"""

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
				# 	if self.x2**2 + self.y2**2 >= H2:
				# 		if self.x3**2 + self.y3**2 >= H2:
				# 			return True
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



def calc_field_from_trans_opt(field, trans, medium, num_of_tris_on_line):
	logging.info("-----Numerical field calculation has been started-----")
	logging.info(
		"Using transducer '%s' and medium '%s'", trans.name, medium.name)
	start_time = time.clock()

	wave_length = medium.speed_of_sound / trans.frequency
	wave_number = 2 * math.pi / wave_length
	kappa = wave_number + 1j * medium.attenuation
	logging.info("Wave length = %f mm", 1.0e3 * wave_length)
	logging.info("Wave number = %f m-1", wave_number)

	# Precalculations for performance
	logging.info("Forming sources...")
	sources = form_sources_from_trans(trans, num_of_tris_on_line)

	# Initial pressure on transducer = 1 Pa
	sources['Ss'] *= 1.0 / (medium.density * medium.speed_of_sound)

	shape = numpy.shape(field.p)
	points = [(field, sources, medium, kappa, i, j, k) for i in range(
		0, shape[0]) for j in range(0, shape[1]) for k in range(0, shape[2])]

	pool = Pool()

	logging.info("Entering main cycle...")
	res = pool.map(process_field_calc, points)

	res = numpy.array(res)
	field.p = numpy.reshape(res[:, 0], shape)
	field.vn = numpy.reshape(res[:, 1], shape)

	duration = time.clock() - start_time
	duration = timedelta(seconds=duration)
	logging.info(
		"-----Numerical field calculation has been finished. It took %s", duration)


def calc_field_from_trans_opt_ideal_trans(field, trans, medium, num_of_tris_on_line, ribs_phantom):
	"""
	Модификация, чтобы считать "идеальный" источник, чтобы посчитать "полосатый" источник в присутствии ребер
	"""
	logging.info("-----Numerical field calculation has been started-----")
	logging.info(
		"Using transducer '%s' and medium '%s'", trans.name, medium.name)
	start_time = time.clock()

	wave_length = medium.speed_of_sound / trans.frequency
	wave_number = 2 * math.pi / wave_length
	kappa = wave_number + 1j * medium.attenuation
	logging.info("Wave length = %f mm", 1.0e3 * wave_length)
	logging.info("Wave number = %f m-1", wave_number)

	# Precalculations for performance
	logging.info("Forming sources...")
	sources = form_sources_from_trans_ribs(
		trans, num_of_tris_on_line, ribs_phantom)

	# Initial pressure on transducer = 1 Pa
	# sources['Ss'] *= 1.0 / (medium.density * medium.speed_of_sound)

	# Power checking
	sources_power = calc_sources_power(sources, medium)
	ideal_power = calc_ideal_power(trans, medium)
	logging.info("Sources power = " + str(sources_power))
	logging.info("Ideal power = " + str(ideal_power))

	shape = numpy.shape(field.p)
	points = [(field, sources, medium, kappa, i, j, k) for i in range(
		0, shape[0]) for j in range(0, shape[1]) for k in range(0, shape[2])]

	pool = Pool()

	logging.info("Entering main cycle...")
	res = pool.map(process_field_calc, points)
	# res = [process_field_calc(point) for point in points]

	res = numpy.array(res)
	field.p = numpy.reshape(res[:, 0], shape)
	field.vn = numpy.reshape(res[:, 1], shape)

	duration = time.clock() - start_time
	duration = timedelta(seconds=duration)
	logging.info(
		"-----Numerical field calculation has been finished. It took %s", duration)


def calc_field_from_trans_opt_full(field, trans, medium, num_of_tris_on_line):
	"""
	Модификация, чтобы считать "идеальный" источник, полный
	"""
	logging.info("-----Numerical field calculation has been started-----")
	logging.info(
		"Using transducer '%s' and medium '%s'", trans.name, medium.name)
	start_time = time.clock()

	wave_length = medium.speed_of_sound / trans.frequency
	wave_number = 2 * numpy.pi / wave_length
	kappa = wave_number + 1j * medium.attenuation
	# kappa = wave_number
	logging.info("Wave length = %f mm", 1.0e3 * wave_length)
	logging.info("Wave number = %f m-1", wave_number)

	# Precalculations for performance
	logging.info("Forming sources...")
	sources = form_sources_from_trans_full(trans, num_of_tris_on_line)

	# Power checking
	sources_power = calc_sources_power(sources, medium)
	ideal_power = calc_ideal_power(trans, medium)
	logging.info("Sources power = " + str(sources_power))
	logging.info("Ideal power = " + str(ideal_power))

	# Initial pressure on transducer = 1 Pa
	# sources['Ss'] *= 1.0 / (medium.density * medium.speed_of_sound)

	shape = numpy.shape(field.p)
	points = [(field, sources, medium, kappa, i, j, k) for i in range(
		0, shape[0]) for j in range(0, shape[1]) for k in range(0, shape[2])]

	pool = Pool()

	logging.info("Entering main cycle...")
	res = pool.map(process_field_calc, points)

	res = numpy.array(res)
	field.p = numpy.reshape(res[:, 0], shape)
	field.vn = numpy.reshape(res[:, 1], shape)

	duration = time.clock() - start_time
	duration = timedelta(seconds=duration)
	logging.info(
		"-----Numerical field calculation has been finished. It took %s", duration)


def calc_field_from_trans_after_interface_opt(field, trans, medium1, medium2, num_of_tris_on_line):
	logging.info(
		"-----Numerical field calculation after interface has been started-----")
	logging.info("Using transducer '%s', medium '%s' for propagation and medium '%s' after interface",
				 trans.name, medium1.name, medium2.name)
	start_time = time.clock()

	wave_length = medium1.speed_of_sound / trans.frequency
	wave_number = 2 * math.pi / wave_length
	kappa = wave_number + 1j * medium1.attenuation
	logging.info("Wave length = %f mm", 1.0e3 * wave_length)
	logging.info("Wave number = %f m-1", wave_number)

	# Precalculations for performance
	logging.info("Forming sources...")
	sources = form_sources_from_trans(trans, num_of_tris_on_line)

	# Initial pressure on transducer = 1 Pa
	sources['Ss'] *= 1.0 / (medium1.density * medium1.speed_of_sound)

	shape = numpy.shape(field.p)
	points = [(field, sources, medium1, medium2, kappa, i, j, k) for i in range(
		0, shape[0]) for j in range(0, shape[1]) for k in range(0, shape[2])]

	pool = Pool()

	logging.info("Entering main cycle...")
	res = pool.map(process_field_calc_interface, points)

	res = numpy.array(res)
	field.p = numpy.reshape(res[:, 0], shape)
	field.vn = numpy.reshape(res[:, 1], shape)

	duration = time.clock() - start_time
	duration = timedelta(seconds=duration)
	logging.info(
		"-----Numerical field calculation has been finished. It took %s", duration)


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
	points = [(field, sources, medium, kappa, i, j, k) for i in range(
		0, shape[0]) for j in range(0, shape[1]) for k in range(0, shape[2])]

	pool = Pool()

	logging.info("Entering main cycle...")
	res = pool.map(process_field_calc, points)

	res = numpy.array(res)
	field.p = numpy.reshape(res[:, 0], shape)
	field.vn = numpy.reshape(res[:, 1], shape)

	duration = time.clock() - start_time
	duration = timedelta(seconds=duration)
	logging.info(
		"-----Numerical field calculation has been finished. It took %s", duration)


def calc_field_from_field_after_interface_opt(field, src_field, medium1, medium2, frequency):
	logging.info("-----Numerical field calculation has been started-----")
	logging.info(
		"Using field as a source, medium '%s' for propagation and medium '%s' after interface", medium1.name, medium2.name)
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
	points = [(field, sources, medium1, medium2, kappa, i, j, k) for i in range(
		0, shape[0]) for j in range(0, shape[1]) for k in range(0, shape[2])]

	pool = Pool()

	logging.info("Entering main cycle...")
	res = pool.map(process_field_calc_interface, points)

	res = numpy.array(res)
	field.p = numpy.reshape(res[:, 0], shape)
	field.vn = numpy.reshape(res[:, 1], shape)

	duration = time.clock() - start_time
	duration = timedelta(seconds=duration)
	logging.info(
		"-----Numerical field calculation has been finished. It took %s", duration)


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
	coses = calc_coses(
		rx, ry, rz, r, field_normal_x, field_normal_y, field_normal_z)

	# check if r crosses the source
	bad_rays_map = check_crosses(rx, ry, rz, sources)
	p, vn = sum_up_sources_rayleigh(
		kappa, r, sources['Ss'] * bad_rays_map, 1.0, coses)

	# without check if r crosses the source
	# p, vn = sum_up_sources_rayleigh(kappa, r, sources['Ss'], 1.0, coses)

	p *= 1j * kappa / (2.0 * numpy.pi)
	vn *= 1 / (2.0 * numpy.pi * medium.density * medium.speed_of_sound)

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
	incident_angles = calc_incident_angles(
		rx, ry, rz, r, field_normal_x, field_normal_y, field_normal_z)
	transmission_coses = calc_transmission_coses(
		incident_angles, medium1, medium2)

	transmission_coefs_p = calc_transmission_coefs_for_pressure(
		medium1, medium2, incident_angles, transmission_coses)
	transmission_coefs_vn = calc_transmission_coefs_for_velocity(
		medium1, medium2, incident_angles, transmission_coses)

	# check if r crosses the source
	bad_rays_map = check_crosses(rx, ry, rz, sources)

	p, vn = sum_up_sources_rayleigh(kappa, r, sources[
									'Ss'] * bad_rays_map, transmission_coefs_p, transmission_coefs_vn * transmission_coses)

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

	surface = {
		'Xs': surface_Xs,
		'Ys': surface_Ys,
		'Zs': surface_Zs,
		'n_Xs': surface_normal_Xs,
		'n_Ys': surface_normal_Ys,
		'n_Zs': surface_normal_Zs,
		'ps': surface_ps,
		'vns': surface_vns
		}
	return surface


def calc_field_from_trans(field, trans, medium, num_of_tris_on_line):
	logging.info("-----Numerical field calculation has been started-----")
	logging.info(
		"Using transducer '%s' and medium '%s'", trans.name, medium.name)
	start_time = time.clock()

	wave_length = medium.speed_of_sound / trans.frequency
	wave_number = 2 * math.pi / wave_length
	kappa = wave_number + 1j * medium.attenuation
	logging.info("Wave length = %f mm", 1.0e3 * wave_length)
	logging.info("Wave number = %f m-1", wave_number)

	# Precalculations for performance
	logging.info("Forming sources...")
	sources = form_sources_from_trans(trans, num_of_tris_on_line)

	# Initial pressure on transducer = 1 Pa
	sources['Ss'] *= 1.0 / (medium.density * medium.speed_of_sound)

	field_calc_main_cycle(field, sources, medium, kappa)

	duration = time.clock() - start_time
	duration = timedelta(seconds=duration)
	logging.info(
		"-----Numerical field calculation has been finished. It took %s", duration)


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
	logging.info(
		"-----Numerical field calculation has been finished. It took %s", duration)


def calc_field_from_trans_after_interface(field, trans, medium1, medium2, num_of_tris_on_line):

	logging.info(
		"-----Numerical field calculation after interface has been started-----")
	logging.info("Using transducer '%s', medium '%s' for propagation and medium '%s' after interface", trans.name, medium1.name, medium2.name)
	start_time = time.clock()

	wave_length = medium1.speed_of_sound / trans.frequency
	wave_number = 2 * math.pi / wave_length
	kappa = wave_number + 1j * medium1.attenuation
	logging.info("Wave length = %f mm", 1.0e3 * wave_length)
	logging.info("Wave number = %f m-1", wave_number)

	# Precalculations for performance
	logging.info("Forming sources...")
	sources = form_sources_from_trans(trans, num_of_tris_on_line)

	# Initial pressure on transducer = 1 Pa
	sources['Ss'] *= 1.0 / (medium1.density * medium1.speed_of_sound)

	field_calc_after_interface_main_cycle(
		field, sources, medium1, medium2, kappa)

	duration = time.clock() - start_time
	duration = timedelta(seconds=duration)
	logging.info(
		"-----Numerical field calculation has been finished. It took %s", duration)


def calc_field_from_field_after_interface(field, src_field, medium1, medium2, frequency):

	logging.info("-----Numerical field calculation has been started-----")
	logging.info(
		"Using field as a source, medium '%s' for propagation and medium '%s' after interface", medium1.name, medium2.name)
	start_time = time.clock()

	wave_length = medium1.speed_of_sound / frequency
	wave_number = 2 * math.pi / wave_length
	kappa = wave_number + 1j * medium1.attenuation
	logging.info("Wave length = %f mm", 1.0e3 * wave_length)
	logging.info("Wave number = %f m-1", wave_number)

	# Precalculations for performance
	logging.info("Forming sources...")
	sources = form_sources_from_field(src_field)

	field_calc_after_interface_main_cycle(
		field, sources, medium1, medium2, kappa)

	duration = time.clock() - start_time
	duration = timedelta(seconds=duration)
	logging.info(
		"-----Numerical field calculation has been finished. It took %s", duration)


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
				field_normal_x, field_normal_y, field_normal_z = field.normal(
					i, j, k)

				rx, ry, rz, r = calc_distance_vectors(
					sources, field_x, field_y, field_z)
				coses = calc_coses(
					rx, ry, rz, r, field_normal_x, field_normal_y, field_normal_z)

				# check if r crosses the source
				bad_rays_map = check_crosses(rx, ry, rz, sources)

				p, vn = sum_up_sources_rayleigh(
					kappa, r, sources['Ss'] * bad_rays_map, 1.0, coses)

				p *= -1j * kappa * medium.density * \
					medium.speed_of_sound / (2.0 * math.pi)
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
				field_normal_x, field_normal_y, field_normal_z = field.normal(
					i, j, k)

				rx, ry, rz, r = calc_distance_vectors(
					sources, field_x, field_y, field_z)
				incident_angles = calc_incident_angles(
					rx, ry, rz, r, field_normal_x, field_normal_y, field_normal_z)
				transmission_coses = calc_transmission_coses(
					incident_angles, medium1, medium2)

				transmission_coefs_p = calc_transmission_coefs_for_pressure(
					medium1, medium2, incident_angles, transmission_coses)
				transmission_coefs_vn = calc_transmission_coefs_for_velocity(
					medium1, medium2, incident_angles, transmission_coses)

				# check if r crosses the source
				bad_rays_map = check_crosses(rx, ry, rz, sources)

				p, vn = sum_up_sources_rayleigh(kappa, r, sources[
												'Ss'] * bad_rays_map, transmission_coefs_p, transmission_coefs_vn * transmission_coses)

				p *= -1j * kappa * medium1.density * \
					medium1.speed_of_sound / (2 * math.pi)
				vn *= -1j * kappa / (2 * math.pi)

				field.p[i, j, k] = p
				field.vn[i, j, k] = vn


def slice_trans_on_tris(trans, num_of_tris_on_line):
	""" Slice the transducer into triangles """
	logging.info(
		"Slicing of transducer into triangles has been started. Number of triangles on line %d", num_of_tris_on_line)

	num_of_rects_on_line = num_of_tris_on_line // 2
	dx = trans.aperture / num_of_rects_on_line
	dy = -dx
	a = trans.aperture / 2.0
	tris = []
	# triArea = 0.0  # Считаем площадь треугольника только один раз

	for j in range(0, num_of_rects_on_line):
		for i in range(0, num_of_rects_on_line):

			# tri1 = Triangle(-a+i*dx, a+float64(j)*dY, -a+float64(i)*dX, a+float64(j+1)*dY, -a+float64(i+1)*dX, a+float64(j)*dY)
			# tri2 = Triangle(-a+float64(i+1)*dX, a+float64(j)*dY, -a+float64(i)*dX, a+float64(j+1)*dY, -a+float64(i+1)*dX, a+float64(j+1)*dY)

			tri1 = Triangle(-a + i * dx, a + j * dy, -a + i *
							dx, a + (j + 1) * dy, -a + (i + 1) * dx, a + j * dy)
			tri2 = Triangle(-a + (i + 1) * dx, a + j * dy, -a + i *
							dx, a + (j + 1) * dy, -a + (i + 1) * dx, a + (j + 1) * dy)

			tri1.compute_centroids()
			tri2.compute_centroids()
			tri1.project_on_trans(trans)
			tri2.project_on_trans(trans)

			if tri1.is_need(trans):
				tri1.project_on_trans(trans)
				tri1.compute_centroids()
				tri1.compute_area()
				tri1.calc_normals_to_trans_focus(trans)
				tris.append(tri1)

			if tri2.is_need(trans):
				tri2.project_on_trans(trans)_
				tri2.compute_centroids()
				tri2.compute_area()
				tri2.calc_normals_to_trans_focus(trans)
				tris.append(tri2)

	logging.info("Slicing of transducer into triangles has been finished. Number of active triangles: %d. Sides of triangle are %fx%f mm", len(
		tris), dx / 1.0e-03, dy / 1.0e-03)

	return tris


def slice_trans_on_tris_ideal(trans, num_of_tris_on_line, ribs_phantom):
	""" Slice the transducer into triangles

	Версия для отключения при ребрах

	"""
	logging.info(
		"Slicing of transducer into triangles with RIBS has been started. Number of triangles on line %d", num_of_tris_on_line)

	num_of_rects_on_line = num_of_tris_on_line // 2
	dx = trans.aperture / num_of_rects_on_line
	dy = -dx
	a = trans.aperture / 2.0
	tris = []
	# triArea = 0.0  # Считаем площадь треугольника только один раз

	for j in range(0, num_of_rects_on_line):
		for i in range(0, num_of_rects_on_line):

			# tri1 = Triangle(-a+i*dx, a+float64(j)*dY, -a+float64(i)*dX, a+float64(j+1)*dY, -a+float64(i+1)*dX, a+float64(j)*dY)
			# tri2 = Triangle(-a+float64(i+1)*dX, a+float64(j)*dY, -a+float64(i)*dX, a+float64(j+1)*dY, -a+float64(i+1)*dX, a+float64(j+1)*dY)

			tri1 = Triangle(-a + i * dx, a + j * dy, -a + i * dx, a + (j + 1) * dy, -a + (i + 1) * dx, a + j * dy)
			tri2 = Triangle(-a + (i + 1) * dx, a + j * dy, -a + i * dx, a + (j + 1) * dy, -a + (i + 1) * dx, a + (j + 1) * dy)
			tri1.project_on_trans(trans)
			tri2.project_on_trans(trans)
			tri1.compute_centroids()
			tri2.compute_centroids()

			if tri1.is_need_full_trans_ribs(trans, ribs_phantom):
				tri1.compute_area()
				tri1.calc_normals_to_trans_focus(trans)
				tris.append(tri1)

			if tri2.is_need_full_trans_ribs(trans, ribs_phantom):
				tri2.compute_area()
				tri2.calc_normals_to_trans_focus(trans)
				tris.append(tri2)

	logging.info("Slicing of transducer into triangles with RIBS has been finished. Number of active triangles: %d. Sides of triangle are %fx%f mm", len(
		tris), dx / 1.0e-03, dy / 1.0e-03)

	return tris


def slice_trans_on_tris_full(trans, num_of_tris_on_line):
	""" Slice the transducer into triangles

	Версия для всего источника

	"""
	logging.info(
		"Slicing of FULL transducer into triangles has been started. Number of triangles on line %d", num_of_tris_on_line)

	num_of_rects_on_line = num_of_tris_on_line // 2
	dx = trans.aperture / num_of_rects_on_line
	dy = -dx
	a = trans.aperture / 2.0
	tris = []
	# triArea = 0.0  # Считаем площадь треугольника только один раз

	for j in range(0, num_of_rects_on_line):
		for i in range(0, num_of_rects_on_line):

			tri1 = Triangle(-a + i * dx, a + j * dy, -a + i *
							dx, a + (j + 1) * dy, -a + (i + 1) * dx, a + j * dy)
			tri2 = Triangle(-a + (i + 1) * dx, a + (j + 1) * dy, -a +
							i * dx, a + (j + 1) * dy, -a + (i + 1) * dx, a + j * dy)

			tri1.compute_centroids()
			tri2.compute_centroids()
			tri1.project_on_trans(trans)
			tri2.project_on_trans(trans)
			tri1.compute_centroids()
			tri2.compute_centroids()

			if tri1.is_need_full_trans_full(trans):
				tri1.compute_area()
				tri1.calc_normals_to_trans_focus(trans)
				tris.append(tri1)

			if tri2.is_need_full_trans_full(trans):
				tri2.compute_area()
				tri2.calc_normals_to_trans_focus(trans)
				tris.append(tri2)

	logging.info("Slicing of FULL transducer into triangles has been finished. Number of active triangles: %d. Sides of triangle are %fx%f mm", len(
		tris), dx / 1.0e-03, dy / 1.0e-03)

	return tris


def form_sources_from_tris(tris):
	sources_Xs = numpy.array([tri.xc for tri in tris])
	sources_Ys = numpy.array([tri.yc for tri in tris])
	sources_Zs = numpy.array([tri.zc for tri in tris])
	sources_Ss = numpy.array([tri.area for tri in tris])

	# These values are valid only in this case of flat transducer's surface!
	# sources_normal_Xs = numpy.zeros_like(sources_Xs)
	# sources_normal_Ys = numpy.zeros_like(sources_Xs)
	# sources_normal_Zs = numpy.ones_like(sources_Xs)

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


def form_sources_from_trans(trans, num_of_tris_on_line):
	tris = slice_trans_on_tris(trans, num_of_tris_on_line)
	return form_sources_from_tris(tris)


def form_sources_from_trans_ribs(trans, num_of_tris_on_line, ribs_phantom):
	"""
	Модификация для ребер
	"""
	tris = slice_trans_on_tris_ideal(trans, num_of_tris_on_line, ribs_phantom)
	import random
	tris = random.sample(tris, len(tris) // 2)
	logging.info("Sampling tris = now " + str(len(tris)))
	save_tris_on_disk(tris)
	# import draw_transducer
	# draw_transducer.draw_transducer_tris(trans, tris)
	return form_sources_from_tris(tris)
	# return None


def form_sources_from_trans_full(trans, num_of_tris_on_line):
	"""
	Модификация для всего источника
	"""
	# tris = slice_trans_on_tris_ideal(trans, num_of_tris_on_line, ribs_phantom)
	tris = slice_trans_on_tris_full(trans, num_of_tris_on_line)
	save_tris_on_disk(tris)
	# import draw_transducer
	# draw_transducer.draw_transducer_tris(trans, tris)
	return form_sources_from_tris(tris)
	# return None


def save_tris_on_disk(tris, file_name="tris.bin"):
	logging.info("Saving tris on disk...")
	with open(file_name, 'wb') as f:
		pickle.dump(tris, f)


def restore_tris_from_disk(file_name="tris.bin"):
	logging.info("Restoring tris...")
	with open(file_name, 'rb') as f:
		tris = pickle.load(f)
	return tris


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

	sources = {
		'Xs': sources_Xs,
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
	# return DX, DY, DZ, numpy.sqrt(DX*DX + DY*DY + DZ*DZ)
	return DX, DY, DZ, numpy.sqrt(numpy.square(DX) + numpy.square(DY) + numpy.square(DZ))


def calc_coses(rx, ry, rz, r, normal_x, normal_y, normal_z):
	return numpy.abs(rx * normal_x + ry * normal_y + rz * normal_z) / r


def check_crosses(rx, ry, rz, sources):
	# Checks if r crosses the source
	# Checking Scalar multiplying normal vector on distance vector
	# if (mult >= 0) => (angle <= 90) => OK, unit value, else (angle > 90) =>
	# Zero value
	mult = rx * sources['n_Xs'] + ry * sources['n_Ys'] + rz * sources['n_Zs']
	return 1.0 * (mult >= 0)


def calc_incident_angles(rx, ry, rz, r, normal_x, normal_y, normal_z):
	coses = calc_coses(rx, ry, rz, r, normal_x, normal_y, normal_z)
	return scipy.arccos(coses)


def calc_transmission_coses(incident_angles, medium1, medium2):
	val = (1 - ((medium2.speed_of_sound / medium1.speed_of_sound) * numpy.sin(incident_angles)) ** 2)
	return scipy.sqrt(val * (val > 0.0))


def calc_transmission_coefs_for_velocity(medium1, medium2, incident_angles, transmission_coses):
	incident_local_impedances = calc_incident_local_impedances(
		medium1, incident_angles)
	transmission_local_impedances = calc_transmission_local_impedances(
		medium2, transmission_coses)

	return ((medium1.density * medium1.speed_of_sound) / (medium2.density * medium2.speed_of_sound)) * 2.0 * transmission_local_impedances / (transmission_local_impedances + incident_local_impedances)


def calc_incident_local_impedances(medium, angles):
	return (medium.density * medium.speed_of_sound) / numpy.cos(angles)


def calc_transmission_local_impedances(medium, transmission_coses):
	return (transmission_coses > 0) * (medium.density * medium.speed_of_sound) / (transmission_coses * (transmission_coses > 0) + 1.0 * (transmission_coses == 0))


def calc_transmission_coefs_for_pressure(medium1, medium2, icident_angles, transmission_coses):
	incident_local_impedances = calc_incident_local_impedances(
		medium1, icident_angles)
	transmission_local_impedances = calc_transmission_local_impedances(
		medium2, transmission_coses)

	return 2.0 * transmission_local_impedances / (transmission_local_impedances + incident_local_impedances)


def sum_up_sources_rayleigh(kappa, r, S, tran_coef_p, tran_coef_vn):
	val = numpy.exp(1j * kappa * r) * (S * numpy.reciprocal(r))
	p = val * tran_coef_p
	vn = val * tran_coef_vn * (1j * kappa - numpy.reciprocal(r))
	return numpy.sum(p), numpy.sum(vn)
	# reciprocal(x) = 1/x
	# return numpy.add.reduce(p), numpy.add.reduce(vn)


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

	# p = p0*k*delta*(2*j1(sin_alpha*k*r))/(sin_alpha*k*r)
	p = p0 * k * delta * (2) / (sin_alpha * k * r) * (j1(sin_alpha * k * r) - (
		(delta / a) ** 2) * jn(3, sin_alpha * k * r) + ((delta / a) ** 4) * jn(5, sin_alpha * k * r))

	field.p[:, 0, 0] = p
	duration = time.clock() - start_time
	duration = timedelta(seconds=duration)
	logging.info(
		"-----Reference field calculation in focal plane has been finished. It took %s", duration)


def test_slicing():
	trans = transducer.Transducer("Test1", 4.0, 3.0, 1.0, 6.0e06, 0.0)
	trans.add_element(0, 0.0, 0.0)
	tris = slice_trans_on_tris(trans, 8)
	# tris = array.sliceOnTrisProper(8)
	assert len(tris) == 7
	print("Wrong number of tris after slicing", len(tris))

	print("Number of tris after slicing", len(tris))


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
	W = I_0 * numpy.sum(sources['Ss'])
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
