# -*- coding: utf-8 -*-

"""
Copyright 2014. Sergey Ilin.
Lab 366, Acoustics Department, Faculty of Physics
Lomonosov Moscow State University

"""

import matplotlib.pyplot as plt
import math
import numpy
import random


def write_elements_to_file(elements, file_path):
	# Сохраняем координаты элементов в файл
	with open(file_path, 'w') as f:
		for el in elements:
			f.write('{0: 7.5f}\t{1: 7.5f}\n'.format(el['x'], el['y']))


def write_trans_params_to_file(name, aperture, curvature_radius, element_radius, frequency, hole_radius, file_path):
	# Сохраняем параметры решетки
	with open(file_path, 'w') as f:
		f.write(name+"\n")
		f.write('{}\n'.format(aperture))
		f.write('{}\n'.format(curvature_radius))
		f.write('{}\n'.format(element_radius))
		f.write('{}\n'.format(frequency))
		f.write('{}'.format(hole_radius))


def generate_elements_rand():
	random.seed()
	elements = []
	external_radius = 85e-03
	internal_radius = 20e-03
	# element_radius = 3.5e-03
	# number_of_elements = 256
	element_radius = 0.875e-03
	number_of_elements = 4096
	# element_radius = 0.875e-03
	# number_of_elements = 4096

	possible_points = []

	possible_area = math.pi*(external_radius*external_radius - internal_radius*internal_radius)
	elements_area = number_of_elements*math.pi*element_radius*element_radius

	# Check if it is possible to position elements due to its size
	print("Spareness is {0:.1f}%.".format(100*elements_area/possible_area))
	if possible_area <= elements_area:
		print("It's impossible to position elements due to its size. Possible area", possible_area, ", elements area", elements_area)
		return

	# Ищем возможные координаты элемента
	ro = numpy.linspace(internal_radius+element_radius, external_radius-element_radius, 400)
	phi = numpy.linspace(0.0, 2*math.pi, 400)

	fail = False
	for r in ro:
		for f in phi:
			fail = False
			for p in possible_points:
				# Проверяем, что элемент не пересекает другие элементы
				if math.hypot(r*math.cos(f)-p['x'], r*math.sin(f)-p['y']) <= (2.0 + 0.1*random.random())*element_radius:
					fail = True
					break
			if fail is False:
				possible_points.append({'x': r*math.cos(f), 'y': r*math.sin(f)})

	if len(possible_points) < number_of_elements:
		print("Unable to position all elements, only have", len(possible_points))
		return

	# Выбираем из возможных только нужное количество элементов случайным образом
	print("Choosing good ones from ", len(possible_points), "elements")
	elements = random.sample(possible_points, number_of_elements)
	print("Now there are ", len(elements), "elements")

	write_elements_to_file(elements, r"e:\Downloads\New_calc\array_elements_gen.txt")

	write_trans_params_to_file(
		name="Generated Mine",
		aperture=external_radius*2.0,
		curvature_radius=130.0e-03,
		element_radius=element_radius,
		frequency=1.0e06,
		hole_radius=internal_radius,
		file_path=r"e:\Downloads\New_calc\array_gen.txt")

	print("Drawing array...")
	# Рисуем решетку
	plt.figure(figsize=(5.0, 5.0), dpi=100, facecolor='w')

	# Draw aperture of trans as a circle
	circle = plt.Circle((0.0, 0.0), radius=external_radius*1000, linewidth=2, facecolor='w', edgecolor='k')
	plt.gca().add_patch(circle)
		
	# Draw central hole
	circle = plt.Circle((0.0, 0.0), radius=internal_radius*1000, linewidth=1, facecolor='w', edgecolor='k')
	plt.gca().add_patch(circle)

	for el in elements:
		circle = plt.Circle((el['x']*1000, el['y']*1000), radius=element_radius*1000, linewidth=1, facecolor='w', edgecolor='k')
		plt.gca().add_patch(circle)

	plt.axis('scaled') # to scale circles properly
	plt.show()


def generate_elements_rand_new():
	# Не активный!!!
	# Хотел сделать быстрый алгоритм в сферических координатах, этот медленный!!!
	random.seed()
	elements = []
	array_radius = 130e-03
	array_aperture = 170e-03
	external_radius = 85e-03
	internal_radius = 20e-03
	element_radius = 3.5e-03
	number_of_elements = 256
	# element_radius = 1.75
	# number_of_elements = 1024
	# element_radius = 0.875e-03
	# number_of_elements = 4096

	possible_points = []

	possible_area = math.pi*(external_radius*external_radius - internal_radius*internal_radius)
	elements_area = number_of_elements*math.pi*element_radius*element_radius

	# Check if it is possible to position elements due to its size
	if possible_area <= elements_area:
		print("It's impossible to position elements due to its size")
		return

	# Ищем возможные координаты элемента
	# Используется сферическая СК с центром в центре сферы решетки, радиус - радиус кривизны решетки

	theta_min = math.asin(internal_radius/array_radius)
	theta_max = math.asin(external_radius/array_radius)

	phi_min = 0.0
	phi_max = 2.0*math.pi

	# Размер элемента по угловым координатам
	element_size_in_theta = math.asin(element_radius/array_radius)

	# Массив углов
	print("Creating arrays")
	thetaS = numpy.linspace(theta_min, theta_max, 40)
	phiS = numpy.linspace(phi_min, phi_max, 360)

	print("Looking for possible points")
	fail = False
	for theta in thetaS:
		for phi in phiS:
			for p in possible_points:
				# Проверяем, что элемент не пересекает другие элементы
				if math.hypot((theta-p['theta'])*math.cos((phi+p['phi'])/2), phi-p['phi']) <= (2.0 + random.random())*element_size_in_theta:
					fail = True
					break
			if fail == False:
				possible_points.append({'phi': phi, 'theta': theta})
			fail = False

	if len(possible_points) < number_of_elements:
		print("Unable to position all elements")
		return
	
	# Выбираем из возможных только нужное количество элементов случайным образом
	print("Choosing good ones from ", len(possible_points), "elements")
	elements = random.sample(possible_points, number_of_elements)
	print("Now there are ", len(elements), "elements")

	# Переделываем из сферических координат в декартовые
	for el in elements:
		el['x'] = array_radius*math.sin(el['theta'])*math.cos(el['phi'])
		el['y'] = array_radius*math.sin(el['theta'])*math.sin(el['phi'])

	write_elements_to_file(elements, "elements_gen.txt")

	print("Drawing array...")
	# Рисуем решетку
	plt.figure(figsize=(5.0, 5.0), dpi=100, facecolor='w')

	# Draw aperture of trans as a circle
	circle = plt.Circle((0.0, 0.0), radius=external_radius*1000, linewidth=2, facecolor='w', edgecolor='k')
	plt.gca().add_patch(circle)
		
	# Draw central hole
	circle = plt.Circle((0.0, 0.0), radius=internal_radius*1000, linewidth=1, facecolor='w', edgecolor='k')
	plt.gca().add_patch(circle)

	for el in elements:
		circle = plt.Circle((el['x']*1000, el['y']*1000), radius=element_radius*1000, linewidth=1, facecolor='w', edgecolor='k')
		plt.gca().add_patch(circle)

	plt.axis('scaled') # to scale circles properly
	plt.show()

if __name__ == '__main__':
	# test case
	generate_elements_rand()
