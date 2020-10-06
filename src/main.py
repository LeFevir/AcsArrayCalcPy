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
import ribs
import rayleigh
import analytic_calc



def main():
    new_work_dir_with_time_stamp()
    prepare_logging()

    array = transducer.transducer_from_file(r"..\array.txt")
    array.add_elements_from_file(r"..\array_elements.txt")
    med = medium.medium_from_file(r"..\water_medium.txt")
    dist_s = array.calc_distances_betweeen_elements()
    print(min(dist_s))
    print(max(dist_s))

    import json
    with open(r"..\ribs.json",'r') as f:
        ribs_phantom = json.loads(f.read())

    calc_field = field.PressureFieldCartesian()
    calc_field.set_nodes_num(241, 1, 301)
    calc_field.set_grid_bottom(-100.0e-03, 0.0e-03, 5.0e-03)
    calc_field.set_grid_top(100.0e-03, 0.0e-03, 130.0e-03)
    calc_field.prepare_grid()

    # ribs.create_switched_off_elements_file(ribs_phantom, array, "switched_off.txt")
    # array.switch_off_elements_from_file(r"..\switched_off.txt")

    # analytic_calc.calc_field_analytically(calc_field, array, med)
    # rayleigh.calc_field_from_trans_opt(calc_field, array, med, 300)

    rayleigh.RayCalc(calc_field, array, med, 500).doit()
    # rayleigh.RayCalcIdealRibs(calc_field, array, med, 300, ribs_phantom).doit()
    # rayleigh.RayCalcIdealFull(calc_field, array, med, 500).doit()
    # rayleigh.calc_field_from_trans_opt_ideal_trans(calc_field, array, med, 2000, ribs_phantom)
    # rayleigh.calc_field_from_trans_opt_full(calc_field, array, med, 300)

    # calc_field = output.restore_field_from_disk("field_plane")
    output.save_field_on_disk(calc_field, "field_plane")
    # ribs.CalcPowerOnRibs(calc_field, ribs_phantom)

def main_skull():
    import analytic_calc_skull
    new_work_dir_with_time_stamp()
    prepare_logging()

    array = transducer.transducer_from_file(r"..\array.txt")
    array.add_elements_from_file(r"..\array_elements.txt")
    med = medium.medium_from_file(r"..\water_medium.txt")
    
    calc_field = field.PressureFieldCartesian()
    # calc_field.set_nodes_num(201, 1, 130)
    # calc_field.set_grid_bottom(-100.0e-03, 0.0e-03, 5.0e-03)
    # calc_field.set_grid_top(100.0e-03, 0.0e-03, 130.0e-03)
    calc_field.set_nodes_num(1, 1, 11)
    calc_field.set_grid_bottom(0.0e-03, 0.0e-03, 125.0e-03)
    calc_field.set_grid_top(0.0e-03, 0.0e-03, 135.0e-03)
    calc_field.prepare_grid()



    # ribs.create_switched_off_elements_file(ribs_phantom, array, "switched_off.txt")
    # array.switch_off_elements_from_file(r"..\switched_off.txt")

    # analytic_calc.calc_field_analytically(calc_field, array, med)
    analytic_calc_skull.calc_field_analytically(calc_field, array, med)
    # analytic_calc.calc_field_analytically(calc_field, array, med)
    # rayleigh.calc_field_from_trans_opt(calc_field, array, med, 300)

    # rayleigh.RayCalc(calc_field, array, med, 500).doit()
    # rayleigh.RayCalcIdealRibs(calc_field, array, med, 300, ribs_phantom).doit()
    # rayleigh.RayCalcIdealFull(calc_field, array, med, 500).doit()
    # rayleigh.calc_field_from_trans_opt_ideal_trans(calc_field, array, med, 2000, ribs_phantom)
    # rayleigh.calc_field_from_trans_opt_full(calc_field, array, med, 300)

    # calc_field = output.restore_field_from_disk("field_plane")
    output.save_field_on_disk(calc_field, "field_plane")
    # ribs.CalcPowerOnRibs(calc_field, ribs_phantom)



def main_one_el_test():
    new_work_dir_with_time_stamp()
    prepare_logging()

    array = transducer.transducer_from_file(r"..\array.txt")
    array.add_elements_from_file(r"..\array_elements.txt")
    med = medium.medium_from_file(r"..\water_medium.txt")
    dist_s = array.calc_distances_betweeen_elements()
    print(min(dist_s))
    print(max(dist_s))

    import json
    with open(r"..\ribs.json",'r') as f:
        ribs_phantom = json.loads(f.read())

    calc_field = field.PressureFieldCartesian()
    calc_field.set_nodes_num(101, 101, 1)
    calc_field.set_grid_bottom(-25.0e-03, -25.0e-03, 130.0e-03)
    calc_field.set_grid_top(25.0e-03, 25.0e-03, 130.0e-03)
    calc_field.prepare_grid()

    # ribs.create_switched_off_elements_file(ribs_phantom, array, "switched_off.txt")
    # array.switch_off_elements_from_file(r"..\switched_off.txt")

    # analytic_calc.calc_field_analytically(calc_field, array, med)
    # rayleigh.calc_field_from_trans_opt(calc_field, array, med, 300)

    rayleigh.RayCalc(calc_field, array, med, 500).doit()
    # rayleigh.RayCalcIdealRibs(calc_field, array, med, 300, ribs_phantom).doit()
    # rayleigh.RayCalcIdealFull(calc_field, array, med, 500).doit()
    # rayleigh.calc_field_from_trans_opt_ideal_trans(calc_field, array, med, 2000, ribs_phantom)
    # rayleigh.calc_field_from_trans_opt_full(calc_field, array, med, 300)

    # calc_field = output.restore_field_from_disk("field_plane")
    output.save_field_on_disk(calc_field, "field_plane")
    # ribs.CalcPowerOnRibs(calc_field, ribs_phantom)


def main_from_field():
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

    field_init = output.restore_field_from_disk(r"..\field_plane_ribs")
    clean_field_on_ribs(field_init, ribs_phantom)
    output.save_field_on_disk(field_init, r"..\field_plane_ribs_clean")

    calc_field = field.PressureFieldCartesian()
    # calc_field.set_nodes_num(121, 121, 1)
    # calc_field.set_grid_bottom(-15.0e-03, -15.0e-03, 130.0e-03)
    # calc_field.set_grid_top(15.0e-03, 15.0e-03, 130.0e-03)
    calc_field.set_nodes_num(241, 241, 1)
    calc_field.set_grid_bottom(-120.0e-03, -120.0e-03, 130.0e-03)
    calc_field.set_grid_top(120.0e-03, 120.0e-03, 130.0e-03)
    # calc_field.set_nodes_num(1, 201, 1)
    # calc_field.set_grid_bottom(0.0e-03, -25.0e-03, 130.0e-03)
    # calc_field.set_grid_top(0.0e-03, 25.0e-03, 130.0e-03)
    # calc_field.set_nodes_num(1000, 1000, 1)
    # calc_field.set_grid_bottom(-100.0e-03, -100.0e-03, 45.0e-03)
    # calc_field.set_grid_top(100.0e-03, 100.0e-03, 45.0e-03)
    calc_field.prepare_grid()

    # ribs.create_switched_off_elements_file(ribs_phantom, array, "switched_off.txt")
    # array.switch_off_elements_from_file(r"..\switched_off.txt")

    # analytic_calc.calc_field_analytically(calc_field, array, med)
    # rayleigh.calc_field_from_trans_opt(calc_field, array, med, 300)

    # rayleigh.RayCalc(calc_field, array, med, 500).doit()
    # rayleigh.RayCalcIdealRibs(calc_field, array, med, 500, ribs_phantom).doit()
    rayleigh.RayCalcFromField(calc_field, array, med, field_init).doit()
    # rayleigh.RayCalcIdealFull(calc_field, array, med, 300).doit()
    # rayleigh.calc_field_from_trans_opt_ideal_trans(calc_field, array, med, 2000, ribs_phantom)
    # rayleigh.calc_field_from_trans_opt_full(calc_field, array, med, 300)

    output.save_field_on_disk(calc_field, "field_plane")
    # ribs.CalcPowerOnRibs(calc_field, ribs_phantom)


def clean_field_on_ribs(field, ribs_phantom):
    import numpy
    import ribs
    sh = numpy.shape(field.p)
    for j in range(0, sh[1]):
        field_x, field_y, field_z = field.get_cartesian_coords(0, j, 0)
        if ribs.is_point_on_ribs(ribs_phantom, field_y):
            field.p[:, j, :] = 0 + 1j*0
            field.vn[:, j, :] = 0 + 1j*0


def main_spectrum_focus():
    import angular_spectrum
    new_work_dir_with_time_stamp()
    prepare_logging()

    array = transducer.transducer_from_file(r"..\array.txt")
    array.add_elements_from_file(r"..\array_elements.txt")
    med = medium.medium_from_file(r"..\water_medium.txt")


    calc_field = field.PressureFieldCartesian()
    calc_field.set_nodes_num(101, 101, 1)
    calc_field.set_grid_bottom(-100.0e-03, -100.0e-03, 130.0e-03)
    calc_field.set_grid_top(100.0e-03, 100.0e-03, 130.0e-03)
    calc_field.prepare_grid()

    # ribs.create_switched_off_elements_file(ribs_phantom, array, "switched_off.txt")
    # array.switch_off_elements_from_file(r"..\switched_off.txt")

    # analytic_calc.calc_field_analytically(calc_field, array, med)
    # rayleigh.calc_field_from_trans_opt(calc_field, array, med, 300)

    # rayleigh.RayCalc(calc_field, array, med, 500).doit()
    # rayleigh.RayCalcIdealRibs(calc_field, array, med, 500, ribs_phantom).doit()
    field_ribs = output.restore_field_from_disk(r"..\field_plane_ribs")
    # Где считать поле
    z = -(130.0e-03 - 45.0e-03)
    p_f, v_f = angular_spectrum.propagate_cartesian(field_ribs.p[:, :, 0], 0.66e-03, 0.66e-03, array.frequency, med.density, med.speed_of_sound, z)

    calc_field.p[:,:,0] = p_f
    calc_field.vn[:,:,0] = v_f
    # rayleigh.RayCalcIdealFull(calc_field, array, med, 300).doit()
    # rayleigh.calc_field_from_trans_opt_ideal_trans(calc_field, array, med, 2000, ribs_phantom)
    # rayleigh.calc_field_from_trans_opt_full(calc_field, array, med, 300)

    output.save_field_on_disk(calc_field, "field_plane")
    # ribs.CalcPowerOnRibs(calc_field, ribs_phantom)


def main_one_el():
    new_work_dir_with_time_stamp()
    prepare_logging()

    # array = transducer.transducer_from_file(r"..\array.txt")
    array = transducer.Transducer("Test", 185.0e-03, 130.0e-03, 6.0e-03, 1.0e06, 10.0e-03)
    # array.add_element("0", 0.0, 0.0)
    med = medium.medium_from_file(r"..\water_medium.txt")

    calc_field = field.PressureFieldCartesian()
    calc_field.set_nodes_num(241, 1, 1)
    calc_field.set_grid_bottom(-3.0e-03, 0.0e-03, 130.0e-03)
    calc_field.set_grid_top(3.0e-03, 0.0e-03, 130.0e-03)
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

    # rayleigh.RayCalc(calc_field, array, med, 1000).doit()
    rayleigh.reference_anal_calc_in_focal_plane(calc_field, array, med)
    # rayleigh.RayCalcIdealRibs(calc_field, array, med, 1300, ribs_phantom).doit()
    # rayleigh.RayCalcIdealFull(calc_field, array, med, 300).doit()
    # rayleigh.calc_field_from_trans_opt_ideal_trans(calc_field, array, med, 2000, ribs_phantom)
    # rayleigh.calc_field_from_trans_opt_full(calc_field, array, med, 300)

    # calc_field = output.restore_field_from_disk("field_plane")
    output.save_field_on_disk(calc_field, "field_plane")
    # ribs.CalcPowerOnRibs(calc_field, ribs_phantom)

def main_ideal_ref():
    new_work_dir_with_time_stamp()
    prepare_logging()

    array = transducer.transducer_from_file(r"..\array.txt")
    # array.add_element("0", 0.0, 0.0)
    med = medium.medium_from_file(r"..\water_medium.txt")

    calc_field = field.PressureFieldCartesian()
    # calc_field.set_nodes_num(1, 1, 801)
    # calc_field.set_grid_bottom(0.0e-03, 0.0e-03, 110.0e-03)
    # calc_field.set_grid_top(0.0e-03, 0.0e-03, 150.0e-03)
    calc_field.set_nodes_num(201, 1, 1)
    calc_field.set_grid_bottom(-5.0e-03, 0.0e-03, 130.0e-03)
    calc_field.set_grid_top(5.0e-03, 0.0e-03, 130.0e-03)
    # calc_field.set_nodes_num(1000, 1000, 1)
    # calc_field.set_grid_bottom(-100.0e-03, -100.0e-03, 45.0e-03)
    # calc_field.set_grid_top(100.0e-03, 100.0e-03, 45.0e-03)
    calc_field.prepare_grid()

    # ribs.create_switched_off_elements_file(ribs_phantom, array, "switched_off.txt")
    # array.switch_off_elements_from_file(r"switched_off.txt")

    # analytic_calc.calc_field_analytically(calc_field, array, med)
    # rayleigh.calc_field_from_trans_opt(calc_field, array, med, 300)

    # rayleigh.RayCalc(calc_field, array, med, 1000).doit()
    # rayleigh.reference_anal_calc_on_axis(calc_field, array, med)
    # rayleigh.reference_anal_calc_in_focal_plane(calc_field, array, med)
    # rayleigh.RayCalcIdealRibs(calc_field, array, med, 1300, ribs_phantom).doit()
    rayleigh.RayCalcIdealFull(calc_field, array, med, 1500).doit()
    # rayleigh.calc_field_from_trans_opt_ideal_trans(calc_field, array, med, 2000, ribs_phantom)
    # rayleigh.calc_field_from_trans_opt_full(calc_field, array, med, 300)

    # calc_field = output.restore_field_from_disk("field_plane")
    output.save_field_on_disk(calc_field, "field_plane")
    # ribs.CalcPowerOnRibs(calc_field, ribs_phantom)

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
    path += "\\" + time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
    os.mkdir(path)
    os.chdir(path)


if __name__ == '__main__':
    # main()
    main_skull()
    # main_one_el()
    # main_spectrum_focus()
    # main_from_field()
    # main_test()
    # main_ideal_ref()
