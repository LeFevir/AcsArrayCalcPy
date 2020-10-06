# -*- coding: utf-8 -*-

"""
Copyright 2012. Sergey Ilin.
Lab 366, Acoustic Department, Faculty of Physics
Lomonosov Moscow State University

"""

import logging
import math
import codecs


class Transducer:
    # Acoustic transducer - array
    def __init__(self, name, aperture, curvature_radius, element_radius, frequency, hole_radius):
        self.name = name
        self.aperture = aperture
        self.curvature_radius = curvature_radius
        self.element_radius = element_radius
        self.frequency = frequency
        self.elements = []
        self.hole_radius = hole_radius

        logging.info('Transducer has been set: ' + str(self))

    def add_element(self, id, center_x, center_y):
        if (math.hypot(center_x, center_y) > self.aperture / 2.0):
            logging.info("Wrong coordinates of adding element with ID "+str(id)+", not corresponds to aperture")
            return

        if (math.hypot(center_x, center_y) > self.curvature_radius):
            logging.info("Wrong coordinates of adding element with ID "+str(id)+", not corresponds to curvature radius")
            return

        # Calculating Z coordinate corresponding to curvature radius of array
        center_z = self.curvature_radius - math.sqrt(self.curvature_radius ** 2 - center_x ** 2 - center_y ** 2)

        """Calculating normal of element's surface corresponding to curvature radius of array
           Calculating vector:
                start point = coordinates of the center of the element (center_x, center_y, center_z),
                end point = coordinates of center of curvature of array (0.0, 0.0, curvature_radius)
            norm of this vector = curvature_radius
        """
        normal_x = (0.0 - center_x) / self.curvature_radius
        normal_y = (0.0 - center_y) / self.curvature_radius
        normal_z = (self.curvature_radius - center_z) / self.curvature_radius

        # Adding new element to the Array of Elements
        element = {'id': id,
                'center_x': center_x,
                'center_y': center_y,
                'center_z': center_z,
                'normal_x': normal_x,
                'normal_y': normal_y,
                'normal_z': normal_z,
                'phase_shift': 0.0
                }
        self.elements.append(element)

        # logging.info('Element of array with ID {} has been added. Center in ({}, {}, {}) mm, normal ({}, {}, {})'.format(
        #     id,
        #     center_x,
        #     center_y,
        #     center_z,
        #     normal_x,
        #     normal_y,
        #     normal_z
        #     ))

    def add_elements_from_file(self, file_path):
        with codecs.open(file_path, encoding='utf-8-sig') as f:
            lines = [line.strip().split() for line in f]

        counter = 0
        for line in lines:
            self.add_element(counter, float(line[0]), float(line[1]))
            counter += 1
        logging.info("Elements of the Array have been read from file, now total: " + str(self.elements_count()))

    def switch_off_element_by_id(self, id):
        if self.elements_count == 0: 
            return

        for i in range(0, self.elements_count()):
            if self.elements[i]['id'] == id:
                self.elements.pop(i)
                return

    def switch_off_elements_from_file(self, file_path):
        with codecs.open(file_path, encoding='utf-8-sig') as f:
            lines = [line.strip() for line in f]

        for line in lines:
            id = int(line)
            self.switch_off_element_by_id(id)

        logging.info("Elements of the Array have been switched off from file, now total: " + str(self.elements_count()))

    def set_focus_to(self, focus_x, focus_y, focus_z):
        # This Function sets the Electronic focus into defined point with fX,fY,fZ coordinated
        for element in self.elements:
            element['phase_shift'] = math.sqrt((focus_x - element['center_x']) ** 2 + (focus_y - element['center_y']) ** 2 + (focus_z - element['center_z']) ** 2)

        logging.info("Phase focus of the array has been put into ({}, {}, {}) mm".format(
            focus_x/1.0e-03,
            focus_y/1.0e-03,
            focus_z/1.0e-03
            ))

    def set_focus_from_file(self, file_path):
        with codecs.open(file_path, encoding='utf-8-sig') as f:
            lines = [line.strip() for line in f]

        if len(lines) != 3:
            logging.info("Couldn't read focus parameters from file " + file_path)
            return
        self.set_focus_to(float(lines[0]), float(lines[1]), float(lines[2]))

    def elements_count(self):
        return len(self.elements)

    def calc_distances_betweeen_elements(self):
        import numpy
        if self.elements_count() == 0:
            return

        x_s = [element['center_x'] for element in self.elements]
        y_s = [element['center_y'] for element in self.elements]
        x_s = numpy.array(x_s)
        y_s = numpy.array(y_s)
        dist_s = []
        
        for i in range(len(self.elements)):
            x = x_s[i]
            y = y_s[i]
            D = (x-x_s)**2 + (y-y_s)**2
            dist = list(numpy.sqrt(D))
            dist.pop(i)
            dist_s.append(min(dist))

        return dist_s
            
    def __str__(self):
        return 'Transducer "{}" has aperture = {} mm, curvature radius = {} mm, each element radius = {} mm, operating frequency = {} MHz, central hole radius = {} mm'.format(
            str(self.name),
            str(self.aperture/1.0e-03),
            str(self.curvature_radius/1.0e-03),
            str(self.element_radius/1.0e-03),
            str(self.frequency/1.0e6),
            str(self.hole_radius/1.0e-03)
            )


def transducer_from_file(file_path):
    # Adding new transducer from file
    # 'utf-8-sig' means that it skips BOM in UTF files created in Notepad
    with codecs.open(file_path, encoding='utf-8-sig') as f:
        lines = [line.strip() for line in f]

    if len(lines) < 5:
        logging.info("Couldn't read array parameters from file " + file_path)
        return

    trans = Transducer(
        name = lines[0],
        aperture = float(lines[1]),
        curvature_radius = float(lines[2]),
        element_radius = float(lines[3]),
        frequency = float(lines[4]),
        hole_radius = float(lines[5])
        )
    return trans


if __name__ == '__main__':
    # test case
    # trans = Transducer("Hand Gavrilov", 100.0, 120.0, 5.0, 6.0e06)
    # print trans

    trans = transducer_from_file(r"..\test_files\array.txt")
    print(trans)
    trans.add_elements_from_file(r"..\test_files\array_elements.txt")
    print(trans.elements_count())
    print(trans.elements[0])
    trans.set_focus_to(0.0, -20.0e-03, 120.0e-03)
    print(trans.elements[0])
