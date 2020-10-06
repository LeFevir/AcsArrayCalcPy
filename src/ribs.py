# -*- coding: utf-8 -*-

"""
Copyright 2012. Sergey Ilin.
Lab 366, Acoustic Department, Faculty of Physics
Lomonosov Moscow State University

"""

import logging
import math
import numpy

def create_switched_off_elements_file(ribs, trans, file_path):
	ids = get_switched_off_elements_ids(ribs, trans)
	with open(file_path, 'w') as f:
		for id in ids:
			f.write(str(id) + '\n')
		logging.info("File with switched off elements has been generated")

def is_point_trans_intersects_ribs(src_x, src_y, src_z, focus, proj_plane_z, ribs_phantom):
	y = calc_projected_y(src_x, src_y, src_z, focus, proj_plane_z)
	return is_point_on_ribs(ribs_phantom, y) 

def calc_projected_y(src_x, src_y, src_z, focus, proj_plane_z):
	"""
	Вычисление координаты проекции на плоскость ребер из координаты на решетке
	"""

	# focus = math.sqrt(focus*focus - src_x*src_x - - src_y*src_y)
	# focuslf.z1 = F - math.sqrt(F*F-x0*x0-y0*y0)

	res = src_y * (proj_plane_z - focus) / (src_z - focus)

	# res = src_y * (focus - proj_plane_z) / (focus - src_z)
	# logging.info('focus = '+ str(focus))
	# logging.info('projplane = ' + str(proj_plane_z))
	# logging.info('src_z = ' + str(src_z))
	# res = math.fabs(src_y) * (focus - proj_plane_z) / (focus - math.fabs(src_z))
	# logging.info('res = ' + str(res))
	return res

def is_point_on_ribs(ribs, coord_y):
	# Checks if the point with Y coordinate is on the ribs
	# Works only with symmetrical ribs position and if there are 7 ribs in total
	# Calculates the distance from the center corresponding to ribs position
	y = math.fabs(coord_y)

	dist = ribs['rib_width'] / 2.0
	if y < dist:
		return True
	else:
		dist = dist + ribs['gap_width']
		if y <= dist:
			return False
		else:
			dist = dist + ribs['rib_width']
			if y < dist:
				return True
			else:
				dist = dist + ribs['gap_width']
				if y <= dist:
					return False
				else:
					dist = dist + ribs['rib_width']
					if y < dist:
						return True
					else:
						dist = dist + ribs['gap_width']
						if y <= dist:
							return False
						else:
							dist = dist + ribs['rib_width']
							if y < dist:
								return True
	return False
	# return True

def is_point_on_ribs_old(ribs, coord_y):
	# Checks if the point with Y coordinate is on the ribs
	# Works only with symmetrical ribs position and if there are 5 ribs in total
	# Calculates the distance from the center corresponding to ribs position
	y = math.fabs(coord_y)

	dist = ribs['rib_width'] / 2.0
	if y <= dist:
		return True

	dist += ribs['gap_width']
	if y <= dist:
		return False

	dist += ribs['rib_width']
	if y <= dist:
		return True

	dist += ribs['gap_width']
	if y <= dist:
		return False

	dist += ribs['rib_width']
	if y <= dist:
		return True

	return False

def get_switched_off_elements_ids(ribs, trans):
	ids = []
	for element in trans.elements:
		# Вычисление координаты проекции элемента решетки на плоскость ребер
		coord_y = calc_projected_y(element['center_x'],element['center_y'], element['center_z'], trans.curvature_radius, ribs['dist_to_ribs_plane'])
		if is_point_on_ribs(ribs, coord_y):
			# Добавление ID элемента в список на отключение
			ids.append(element['id'])

	return ids

def calc_poynting_in_point(field, i, j, k):
	# Calculating of Poynting Vector in point i,j,k of field field
	# P = ( P * V.conjugate ).real
	val = 0.5 * (field.p[i,j,k] * field.vn[i,j,k].conjugate()).real 
	return val

def CalcPowerOnRibs(field, ribs):
	powerOnPlane = 0.0
	powerOnRibs = 0.0
	w = 0.0

	for j, y in enumerate(field.y):
		pointOnRibs = is_point_on_ribs(ribs, y)
		for i, x in enumerate(field.x):
			# w = field.p_amp(i, j, 0) * field.p_amp(i, j, 0)
			w = calc_poynting_in_point(field, i, j, 0)
			powerOnPlane += w
			if pointOnRibs:
				powerOnRibs += w
	
	powerOnPlane *= field.one_pixel_area
	powerOnRibs *= field.one_pixel_area

	logging.info("PowerOnPlane = " + str(powerOnPlane))
	logging.info("PowerOnRibs = " + str(powerOnRibs))
	logging.info("Ratio = " + str(powerOnRibs/powerOnPlane))


def CalcPowerOnPlane(field):
	powerOnPlane = 0.0
	w = 0.0

	for j, y in enumerate(field.y):
		for i, x in enumerate(field.x):
			# w = field.p_amp(i, j, 0) * field.p_amp(i, j, 0)
			w = calc_poynting_in_point(field, i, j, 0)
			powerOnPlane += w
			
	powerOnPlane *= field.one_pixel_area
	# powerOnRibs *= field.one_pixel_area

	print("PowerOnPlane = " + str(powerOnPlane))
	# logging.info("PowerOnRibs = " + str(powerOnRibs))
	# logging.info("Ratio = " + str(powerOnRibs/powerOnPlane))


def CalcPowerOnRibs_sveta(sveta_file, ribs):
	with open(sveta_file, 'r') as f:
		Z = numpy.loadtxt(f)
	Z = numpy.transpose(Z)

	y_ar = numpy.linspace(-100, 100, num=2001)
	x_ar = numpy.linspace(-100, 100, num=2001)

	powerOnPlane = 0.0
	powerOnRibs = 0.0
	w = 0.0

	for j, y in enumerate(y_ar):
		pointOnRibs = is_point_on_ribs(ribs, y)
		for i, x in enumerate(x_ar):
			# w = field.p_amp(i, j, 0) * field.p_amp(i, j, 0)
			# w = calc_poynting_in_point(field, i, j, 0)
			# w = Z[i,j]
			w = Z[i,j]
			powerOnPlane += w
			if pointOnRibs:
				powerOnRibs += w
	
	powerOnPlane *= 0.1*0.1*1.0e-06
	powerOnRibs *= 0.1*0.1*1.0e-06

	# logging.info("PowerOnPlane = " + str(powerOnPlane))
	# logging.info("PowerOnRibs = " + str(powerOnRibs))
	# logging.info("Ratio = " + str(powerOnRibs/powerOnPlane))
	print("PowerOnPlane = " + str(powerOnPlane))
	print("PowerOnRibs = " + str(powerOnRibs))
	print("Ratio = " + str(powerOnRibs/powerOnPlane))

def rib_test():
	import numpy

	ribs_phantom = {
	# Ребра - плоские горизонтальные полоски одинаковой толщины
	
	# Расстояние до плоскости ребер от решетки
	'dist_to_ribs_plane': 45.0,

	# Количество ребер
	'ribs_count': 5,

	# Ширина ребра
	'rib_width': 18.0,

	# Ширина щели между ребрами
	'gap_width': 14.0,
	}

	F = 130
	rib_width = 18
	gap_width = 14
	aperture = 170

	y_s = numpy.linspace(-85.0, 85.0, 340)
	z_s = F - numpy.sqrt(F*F-y_s*y_s)
	x_s = numpy.ones_like(y_s)

	for i in range(0, len(y_s)):
		if is_point_trans_intersects_ribs(0.0, y_s[i], z_s[i], F, 45, ribs_phantom):
			x_s[i] = 0.0

	import matplotlib.pyplot as plt
	plt.plot(y_s, x_s)
	plt.axis([-85, 85, -1, 2])
	# plt.plot(y_s, z_s)
	# plt.gca().set_aspect('equal')
	plt.show()

	

if __name__ == '__main__':
	pass
	# rib_test()
	# ribs_phantom = {
	# # Ребра - плоские горизонтальные полоски одинаковой толщины
	
	# # Расстояние до плоскости ребер от решетки
	# 'dist_to_ribs_plane': 45.0,

	# # Количество ребер
	# 'ribs_count': 5,

	# # Ширина ребра
	# 'rib_width': 18.0,

	# # Ширина щели между ребрами
	# 'gap_width': 14.0,
	# }
	# CalcPowerOnRibs_sveta(r"c:\Downloads\Calc\Sveta_res\intens.txt",ribs_phantom)

	# calc
	# # test case
	# ribs_phantom = {
	# # Ребра - плоские горизонтальные полоски одинаковой толщины
	
	# # Расстояние до плоскости ребер от решетки
	# 'dist_to_ribs_plane': 45.0e-03,

	# # Количество ребер
	# 'ribs_count': 5,

	# # Ширина ребра
	# 'rib_width': 18.0e-03,

	# # Ширина щели между ребрами
	# 'gap_width': 14.0e-03,

	# # Координата Y нижней грани нижнего ребра
	# 'bottom_coord': 14.0e-03
	# }

	# y = 84.0e-03
	# z = 10.0e-03
	# x=0.0
	# res_y = calc_projected_y(x, y, z, 130e-03, 45.0e-03)
	# print('Coord = ',y,z,', so RES = ',res_y)	

	# # res_y = 70.1e-03
	# res = is_point_on_ribs(ribs_phantom, res_y)
	# print('Coord = ',res_y,', so RES = ',res)

	# # res_y_s = range(0.0, 85e-03,0.1e-03)

	# # res_s = [is_point_on_ribs(ribs_phantom, res_y) for res_y in res_y_s]

	# # plt.plot(res_y_s, res_s)
	# # plt.show()
