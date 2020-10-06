# -*- coding: utf-8 -*-
"""
Main Module
"""
import transducer
import field
import medium
import output
import logging
import sys
import os
import time
import analytic_calc
import ribs
import rayleigh

def main_test_num():
	new_work_dir_with_time_stamp()
	prepare_logging()

	array = transducer.transducer_from_file(r"..\array.txt")
	array.add_elements_from_file(r"..\array_elements.txt")
	# array.switch_off_elements_from_file(ur"d:\Educ\АКУСТИКА\.Current\Расчет ак поля решетки\CalcPy\HandOffElements.txt")
	# array.switch_off_elements_from_file(ur"d:\Downloads Projects\Расчет для статьи\switched_off.txt")

	med = medium.medium_from_file(r"..\medium.txt")

	calc_field = field.PressureFieldCartesian()
	calc_field.set_nodes_num(341, 341, 1)
	calc_field.set_grid_bottom(-100.0e-03, -100.0e-03, 45.0e-03)
	calc_field.set_grid_top(100.0e-03, 100.0e-03, 45.0e-03)
	# calc_field.set_nodes_num(1, 1, 131)
	# calc_field.set_grid_bottom(0.0e-03, 0.0e-03, 1.0e-03)
	# calc_field.set_grid_top(0.0e-03, 0.0e-03, 130.0e-03)
	calc_field.prepare_grid()

	# # calc_field = field.PressureFieldCartesian()
	# # calc_field.set_nodes_num(1, 121, 121)
	# # calc_field.set_grid_bottom(0.0e-03, -60.0e-03, 70.0e-03)
	# # calc_field.set_grid_top(0.0e-03, 60.0e-03, 190.0e-03)
	# # calc_field.prepare_grid()

	# array.set_focus_from_file(r"d:\Educ\АКУСТИКА\.Current\Расчет ак поля решетки\CalcPy\focus.txt")

	# -------------------------CALCULATIONS-----------------------------------------
	# analytic_calc.calc_field_analytically(calc_field, array, med)
	rayleigh.calc_field_from_trans_opt(calc_field, array, med, 300)
	output.save_field_on_disk(calc_field, "field_plane")
	# output.PrintCylAmpPhaseZ_Binary(calc_field, 1, 1)

def main():
	new_work_dir_with_time_stamp()
	prepare_logging()

	array = transducer.transducer_from_file(r"d:\Downloads Projects\Расчет для статьи\array.txt")
	array.add_elements_from_file(r"d:\Downloads Projects\Расчет для статьи\array_elements.txt")
	# array.switch_off_elements_from_file(ur"d:\Educ\АКУСТИКА\.Current\Расчет ак поля решетки\CalcPy\HandOffElements.txt")
	# array.switch_off_elements_from_file(ur"d:\Downloads Projects\Расчет для статьи\switched_off.txt")

	med = medium.medium_from_file(r"d:\Downloads Projects\Расчет для статьи\medium.txt")

	calc_field = field.PressureFieldCartesian()
	# calc_field.set_nodes_num(341, 341, 1)
	# calc_field.set_grid_bottom(-85.0e-03, -85.0e-03, 45.0e-03)
	# calc_field.set_grid_top(85.0e-03, 85.0e-03, 45.0e-03)
	calc_field.set_nodes_num(1, 1, 131)
	calc_field.set_grid_bottom(0.0e-03, 0.0e-03, 1.0e-03)
	calc_field.set_grid_top(0.0e-03, 0.0e-03, 130.0e-03)
	calc_field.prepare_grid()

	# # calc_field = field.PressureFieldCartesian()
	# # calc_field.set_nodes_num(1, 121, 121)
	# # calc_field.set_grid_bottom(0.0e-03, -60.0e-03, 70.0e-03)
	# # calc_field.set_grid_top(0.0e-03, 60.0e-03, 190.0e-03)
	# # calc_field.prepare_grid()

	array.set_focus_from_file(r"d:\Educ\АКУСТИКА\.Current\Расчет ак поля решетки\CalcPy\focus.txt")

	# -------------------------CALCULATIONS-----------------------------------------
	analytic_calc.calc_field_analytically(calc_field, array, med)
	output.save_field_on_disk(calc_field, "field_plane")
	# output.PrintCylAmpPhaseZ_Binary(calc_field, 1, 1)

def main_one_el():
	new_work_dir_with_time_stamp()
	prepare_logging()

	array = transducer.transducer_from_file(r"..\\array.txt")
	array.add_element("0", 0.0, 0.0)
	med = medium.medium_from_file(r"..\water_medium.txt")
	# array.add_elements_from_file(r"array_elements.txt")
	# array.switch_off_elements_from_file(ur"d:\Educ\АКУСТИКА\.Current\Расчет ак поля решетки\CalcPy\HandOffElements.txt")
	# array.switch_off_elements_from_file(ur"d:\Downloads Projects\Расчет для статьи\switched_off.txt")

	calc_field = field.PressureFieldCartesian()
	# calc_field.set_nodes_num(341, 341, 1)
	# calc_field.set_grid_bottom(-85.0e-03, -85.0e-03, 45.0e-03)
	# calc_field.set_grid_top(85.0e-03, 85.0e-03, 45.0e-03)
	calc_field.set_nodes_num(1, 1, 2601)
	calc_field.set_grid_bottom(0.0e-03, 0.0e-03, 0.1e-03)
	calc_field.set_grid_top(0.0e-03, 0.0e-03, 150.0e-03)
	calc_field.prepare_grid()

	# # calc_field = field.PressureFieldCartesian()
	# # calc_field.set_nodes_num(1, 121, 121)
	# # calc_field.set_grid_bottom(0.0e-03, -60.0e-03, 70.0e-03)
	# # calc_field.set_grid_top(0.0e-03, 60.0e-03, 190.0e-03)
	# # calc_field.prepare_grid()

	# array.set_focus_from_file(r"d:\Educ\АКУСТИКА\.Current\Расчет ак поля решетки\CalcPy\focus.txt")

	# # -------------------------CALCULATIONS-----------------------------------------
	# analytic_calc.calc_field_analytically(calc_field, array, med)
	analytic_calc.calc_exact_on_axis(calc_field, array, med)

	output.save_field_on_disk(calc_field, "field_Z")
	# output.PrintCylAmpPhaseZ_Binary(calc_field, 1, 1)


def main_1d():
	new_work_dir_with_time_stamp()
	prepare_logging()

	array = transducer.transducer_from_file(r"..\\array.txt")
	array.add_elements_from_file(r"..\\array_elements.txt")
	med = medium.medium_from_file(r"..\water_medium.txt")
	# array.add_elements_from_file(r"array_elements.txt")
	# array.switch_off_elements_from_file(ur"d:\Educ\АКУСТИКА\.Current\Расчет ак поля решетки\CalcPy\HandOffElements.txt")
	# array.switch_off_elements_from_file(ur"d:\Downloads Projects\Расчет для статьи\switched_off.txt")

	calc_field = field.PressureFieldCartesian()
	# calc_field.set_nodes_num(341, 341, 1)
	# calc_field.set_grid_bottom(-85.0e-03, -85.0e-03, 45.0e-03)
	# calc_field.set_grid_top(85.0e-03, 85.0e-03, 45.0e-03)
	calc_field.set_nodes_num(1, 1, 1001)
	calc_field.set_grid_bottom(0.0e-03, 0.0e-03, 100.0e-03)
	calc_field.set_grid_top(0.0e-03, 0.0e-03, 150.0e-03)
	calc_field.prepare_grid()

	# # calc_field = field.PressureFieldCartesian()
	# # calc_field.set_nodes_num(1, 121, 121)
	# # calc_field.set_grid_bottom(0.0e-03, -60.0e-03, 70.0e-03)
	# # calc_field.set_grid_top(0.0e-03, 60.0e-03, 190.0e-03)
	# # calc_field.prepare_grid()

	# array.set_focus_from_file(r"d:\Educ\АКУСТИКА\.Current\Расчет ак поля решетки\CalcPy\focus.txt")

	# # -------------------------CALCULATIONS-----------------------------------------
	# analytic_calc.calc_field_analytically(calc_field, array, med)
	rayleigh.calc_field_from_trans_opt(calc_field, array, med, 1000)
	output.save_field_on_disk(calc_field, "field_Z")
	# output.PrintCylAmpPhaseZ_Binary(calc_field, 1, 1)


def main_y():
	prepare_logging()

	array = transducer.transducer_from_file(r"array.txt")
	array.add_elements_from_file(r"array_elements.txt")
	med = medium.medium_from_file(r"water_medium.txt")
	calc_field = field.PressureFieldCartesian()
	calc_field.set_nodes_num(1, 401, 1)
	calc_field.set_grid_bottom(0.0e-03, -100.0e-03, 45.0e-03)
	calc_field.set_grid_top(0.0e-03, 100.0e-03, 45.0e-03)
	calc_field.prepare_grid()

	ribs_phantom = {
	# Ребра - плоские горизонтальные полоски одинаковой толщины
	
	# Расстояние до плоскости ребер от решетки
	'dist_to_ribs_plane': 45.0e-03,

	# Количество ребер
	'ribs_count': 5,

	# Ширина ребра
	'rib_width': 18.0e-03,

	# Ширина щели между ребрами
	'gap_width': 14.0e-03,

	# Координата Y нижней грани нижнего ребра
	'bottom_coord': 14.0e-03
	}
	
	# ribs.create_switched_off_elements_file(ribs_phantom, array, "switched_off.txt")
	array.switch_off_elements_from_file(r"switched_off.txt")
	analytic_calc.calc_field_analytically(calc_field, array, med)
	output.save_field_on_disk(calc_field, "field_2")

def ribs_calc():
	new_work_dir_with_time_stamp()
	prepare_logging()

	array = transducer.transducer_from_file(r"..\array.txt")
	array.add_elements_from_file(r"..\array_elements.txt")
	med = medium.medium_from_file(r"..\water_medium.txt")

	ribs_phantom = {
	# Ребра - плоские горизонтальные полоски одинаковой толщины
	
	# Расстояние до плоскости ребер от решетки
	'dist_to_ribs_plane': 45.0e-03,
	# 'dist_to_ribs_plane': 65.0e-03,

	# Количество ребер
	'ribs_count': 5,

	# Ширина ребра
	'rib_width': 18.0e-03,

	# Ширина щели между ребрами
	'gap_width': 14.0e-03,

	# Координата Y нижней грани нижнего ребра
	'bottom_coord': 14.0e-03
	}

	calc_field = field.PressureFieldCartesian()
	calc_field.set_nodes_num(201, 201, 1)
	calc_field.set_grid_bottom(-100.0e-03, -100.0e-03, 45.0e-03)
	calc_field.set_grid_top(100.0e-03, 100.0e-03, 45.0e-03)
	# calc_field.set_nodes_num(341, 341, 1)
	# calc_field.set_grid_bottom(-85.0e-03, -85.0e-03, 45.0e-03)
	# calc_field.set_grid_top(85.0e-03, 85.0e-03, 45.0e-03)
	# calc_field.set_nodes_num(1000, 1000, 1)
	# calc_field.set_grid_bottom(-100.0e-03, -100.0e-03, 45.0e-03)
	# calc_field.set_grid_top(100.0e-03, 100.0e-03, 45.0e-03)
	calc_field.prepare_grid()

	# ribs.create_switched_off_elements_file(ribs_phantom, array, "switched_off.txt")
	# array.switch_off_elements_from_file(r"switched_off.txt")
	
	# analytic_calc.calc_field_analytically(calc_field, array, med)
	# rayleigh.calc_field_from_trans_opt(calc_field, array, med, 300)
	
	rayleigh.calc_field_from_trans_opt_ideal_trans(calc_field, array, med, 2000, ribs_phantom)
	# rayleigh.calc_field_from_trans_opt_full(calc_field, array, med, 300)

	
	# calc_field = output.restore_field_from_disk("field_plane")
	output.save_field_on_disk(calc_field, "field_plane")
	ribs.CalcPowerOnRibs(calc_field, ribs_phantom)

def main_num_test():
	new_work_dir_with_time_stamp()
	prepare_logging()

	array = transducer.transducer_from_file(r"..\array.txt")
	med = medium.medium_from_file(r"..\water_medium.txt")

	calc_field = field.PressureFieldCartesian()
	# calc_field.set_nodes_num(501, 1, 1)
	# calc_field.set_grid_bottom(-30.0e-03, 0.0e-03, 130.0e-03)
	# calc_field.set_grid_top(30.0e-03, 0.0e-03, 130.0e-03)
	calc_field.set_nodes_num(1, 1, 1623)
	calc_field.set_grid_bottom(0.0e-03, 0.0e-03, 30.0e-03)
	calc_field.set_grid_top(0.0e-03, 0.0e-03, 232.0e-03)
	calc_field.prepare_grid()

	rayleigh.calc_field_from_trans_opt_full(calc_field, array, med, 300)
	# rayleigh.reference_anal_calc_in_focal_plane(calc_field, array, med)
	# rayleigh.reference_anal_calc_on_axis(calc_field, array, med)
	output.save_field_on_disk(calc_field, "field_plane")


def prepare_logging():
	logging.basicConfig(
		filename='log.txt',
		filemode='w',
		format='%(asctime)s %(message)s',
		datefmt='%Y/%m/%d %H:%M:%S',
		level=logging.INFO
		)
	console = logging.StreamHandler(stream=sys.stdout)
	console.setLevel(logging.INFO)
	console.setFormatter(logging.Formatter('%(message)s'))
	logging.getLogger('').addHandler(console)


def new_work_dir_with_time_stamp():
	path = os.getcwd()
	path += "\\"+ time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
	os.mkdir(path)
	os.chdir(path)


if __name__ == '__main__':
	# import cProfile
	# cProfile.run('test_opt()')
	# test_opt()
	# tpx_new_scheme_interfaces()
	# rexolite_old_scheme_interfaces()
	# ellipse_interfaces()
	# main()
	ribs_calc()
	# main_num_test()
	# main_y()
	# main_one_el()
	# main_1d()
	# main_test_num()
