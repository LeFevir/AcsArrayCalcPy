# -*- coding: utf-8 -*-
"""
Module for field class
"""
import numpy
import math
import logging


class PressureFieldCartesian:

    def set_nodes_num(self, n_x, n_y, n_z):
        self.n_x = n_x
        self.n_y = n_y
        self.n_z = n_z

    def set_grid_bottom(self, min_x, min_y, min_z):
        self.min_x = min_x
        self.min_y = min_y
        self.min_z = min_z

    def set_grid_top(self, max_x, max_y, max_z):
        self.max_x = max_x
        self.max_y = max_y
        self.max_z = max_z

    def prepare_grid(self):
        if self.n_x == 1:
            self.x = [self.min_x]
            self.dx = 0.0
        else:
            self.x, self.dx = numpy.linspace(self.min_x, self.max_x, self.n_x, retstep=True)

        if self.n_y == 1:
            self.y = [self.min_y]
            self.dy = 0.0
        else:
            self.y, self.dy = numpy.linspace(self.min_y, self.max_y, self.n_y, retstep=True)

        if self.n_z == 1:
            self.z = [self.min_z]
            self.dz = 0.0
        else:
            self.z, self.dz = numpy.linspace(self.min_z, self.max_z, self.n_z, retstep=True)

        # Area of one element of grid
        self.one_pixel_area = self.dx * self.dy

        self.p = numpy.zeros((self.n_x, self.n_y, self.n_z), dtype=complex)
        self.vn = numpy.zeros((self.n_x, self.n_y, self.n_z), dtype=complex)

        logging.info('Grid has been set: ' + str(self))

    def p_amp(self, i, j, k):
        return numpy.absolute(self.p[i, j, k])

    def p_phase(self, i, j, k):
        return numpy.angle(self.p[i, j, k])

    def vn_amp(self, i, j, k):
        return numpy.absolute(self.vn[i, j, k])

    def vn_phase(self, i, j, k):
        return numpy.angle(self.vn[i, j, k])

    def __str__(self):
        return 'Grid in cartesian coordinates (x, y, z) has {}x{}x{} nodes with steps ({}, {}, {} mm), from ({}, {}, {} mm) to ({}, {}, {} mm). One pixel area is {} mm^2'.format(
            str(self.n_x),
            str(self.n_y),
            str(self.n_z),
            self.dx * 1.0e03,
            self.dy * 1.0e03,
            self.dz * 1.0e03,
            self.min_x * 1.0e03,
            self.min_y * 1.0e03,
            self.min_z * 1.0e03,
            self.max_x * 1.0e03,
            self.max_y * 1.0e03,
            self.max_z * 1.0e03,
            self.one_pixel_area * 1.0e+06
            )

    def get_cartesian_coords(self, i, j, k):
        """i, j, k - point in cartesian coordinates (X, Y, Z)"""
        return self.x[i], self.y[j], self.z[k]

    def normal(self, i, j, k):
        """Generates a normal vector in point (i, j, k)"""
        return 0.0, 0.0, 1.0


class PressureFieldCyl:

    def set_nodes_num(self, n_ro, n_phi, n_z):
        self.n_ro = n_ro
        self.n_phi = n_phi
        self.n_z = n_z

    def set_grid_bottom(self, min_ro, min_phi, min_z):
        self.min_ro = min_ro
        self.min_phi = min_phi
        self.min_z = min_z

    def set_grid_top(self, max_ro, max_phi, max_z):
        self.max_ro = max_ro
        self.max_phi = max_phi
        self.max_z = max_z

    def prepare_grid(self):
        if self.n_ro == 1:
            self.ro = [self.min_ro]
            self.dro = 0.0
        else:
            self.ro, self.dro = numpy.linspace(self.min_ro, self.max_ro, self.n_ro, retstep=True)

        if self.n_phi == 1:
            self.phi = [self.min_phi]
            self.dphi = 0.0
        else:
            self.phi, self.dphi = numpy.linspace(self.min_phi, self.max_phi, self.n_phi, retstep=True)

        if self.n_z == 1:
            self.z = [self.min_z]
            self.dz = 0.0
        else:
            self.z, self.dz = numpy.linspace(self.min_z, self.max_z, self.n_z, retstep=True)

        # Area of one element of grid
        self.one_pixel_area = self.dphi * self.ro[0] * self.dz

        self.p = numpy.zeros((self.n_ro, self.n_phi, self.n_z), dtype=complex)
        self.vn = numpy.zeros((self.n_ro, self.n_phi, self.n_z), dtype=complex)

        logging.info('Grid has been set: ' + str(self))

    def __str__(self):
        return 'Grid in cylindrical coordinates (ro, phi, z) has {}x{}x{} nodes with steps ({} mm, {} rad, {} mm), from ({} mm, {} rad, {} mm) to ({} mm, {} rad, {} mm). One pixel area is {} mm^2'.format(
            str(self.n_ro),
            str(self.n_phi),
            str(self.n_z),
            self.dro * 1.0e03,
            self.dphi,
            self.dz * 1.0e03,
            self.min_ro * 1.0e03,
            self.min_phi,
            self.min_z * 1.0e03,
            self.max_ro * 1.0e03,
            self.max_phi,
            self.max_z * 1.0e03,
            self.one_pixel_area * 1.0e+06
            )

    def p_amp(self, i, j, k):
        return numpy.absolute(self.p[i, j, k])

    def p_phase(self, i, j, k):
        return numpy.angle(self.p[i, j, k])

    def vn_amp(self, i, j, k):
        return numpy.absolute(self.vn[i, j, k])

    def vn_phase(self, i, j, k):
        return numpy.angle(self.vn[i, j, k])

    # TODO Make normal shift, not only on Z
    def set_shift_cartesian_z(self, dz):
        self.shift_z = dz

    def get_cartesian_coords(self, i, j, k):
        """i, j, k - point in cylindrical coordinates (Ro, Phi, Z)"""

        x = self.ro[i] * math.sin(self.phi[j])
        y = self.z[k]
        z = self.ro[i] * math.cos(self.phi[j]) - self.shift_z
        return x, y, z

    def normal(self, i, j, k):
        """Generates a normal vector in point (i, j, k)"""
        normal_x = math.sin(self.phi[j])
        normal_y = 0.0
        normal_z = math.cos(self.phi[j])
        return normal_x, normal_y, normal_z




class FieldElliptic:
    """Grid for elliptic surface.
        b = b
        Theta - eccentric angle
        Y - coordinate, orthogonal to the ellipse's plane

        a - semi-major axis
        b - semi-minor axis
        c - focus coordinate
        e - eccentricity of an ellipse

        x = a * cos(theta)
        y = b * sin(theta)
        e = c/a

        a = sqrt(b**2 + c**2)
    """

    def set_nodes_num(self, n_theta, n_y):
        self.n_b = 1
        self.n_theta = n_theta
        self.n_y = n_y

    def set_grid_bottom(self, min_b, min_theta, min_y):
        self.min_b = min_b
        self.min_theta = min_theta
        self.min_y = min_y

    def set_grid_top(self, max_b, max_theta, max_y):
        self.max_b = max_b
        self.max_theta = max_theta
        self.max_y = max_y

    def set_focus(self, c):
        self.c = c

    def prepare_grid(self):
        self.b = [self.min_b]
        self.db = 0.0

        if self.n_theta == 1:
            self.theta = [self.min_theta]
            self.dtheta = 0.0
        else:
            self.theta, self.dtheta = numpy.linspace(self.min_theta, self.max_theta, self.n_theta, retstep=True)

        if self.n_y == 1:
            self.y = [self.min_y]
            self.dy = 0.0
        else:
            self.y, self.dy = numpy.linspace(self.min_y, self.max_y, self.n_y, retstep=True)

        # Calc parameter of ellipse
        self.a = math.sqrt(self.b[0]**2 + self.c**2)
        # Perimeter of ellipse apprx
        L = 4 * (math.pi*self.a*self.b[0] + (self.a - self.b[0])**2) / (self.a + self.b[0])
        logging.info('Perimeter of ellipse: ' + str(L))
        # Area of one element of grid
        self.one_pixel_area = self.dy * L * abs(self.max_theta - self.min_theta) / (2.0 * math.pi * self.dtheta)

        self.p = numpy.zeros((self.n_b, self.n_theta, self.n_y), dtype=complex)
        self.vn = numpy.zeros((self.n_b, self.n_theta, self.n_y), dtype=complex)

        logging.info('Grid has been set: ' + str(self))

    def __str__(self):
        return 'Grid in elliptic coordinate system (b, theta, y) has {}x{}x{} nodes with steps ({} mm, {} rad, {} mm), from ({} mm, {} rad, {} mm) to ({} mm, {} rad, {} mm). One pixel area is {} mm^2. Parameters of ellipse: a = {} mm, b = {}mm, c = {} mm, e = {} mm'.format(
            str(self.n_b),
            str(self.n_theta),
            str(self.n_y),
            self.db * 1.0e03,
            self.dtheta,
            self.dy * 1.0e03,
            self.min_b * 1.0e03,
            self.min_theta,
            self.min_y * 1.0e03,
            self.max_b * 1.0e03,
            self.max_theta,
            self.max_y * 1.0e03,
            self.one_pixel_area * 1.0e+06,
            self.a * 1.0e03,
            self.min_b * 1.0e03,
            self.c * 1.0e03,
            (self.c/self.a)
            )


    def p_amp(self, i, j, k):
        return numpy.absolute(self.p[i, j, k])

    def p_phase(self, i, j, k):
        return numpy.angle(self.p[i, j, k])

    def vn_amp(self, i, j, k):
        return numpy.absolute(self.vn[i, j, k])

    def vn_phase(self, i, j, k):
        return numpy.angle(self.vn[i, j, k])

    def get_cartesian_coords(self, i, j, k):
        a = self.a
        b = self.b[i]
        r = a * b / math.sqrt((b*math.cos(self.theta[j]))**2 + (a*math.sin(self.theta[j]))**2)
        x = r * math.cos(self.theta[j])
        y = self.y[k]
        z = r * math.sin(self.theta[j])
        return x, y, z

    def normal(self, i, j, k):
        """Generates a normal vector in point (i, j, k)
            nz = (a^2 * z0 / (b^2 * x0)) * nx
            nx^2 + nz^2 = 1

        """
        x0, y0, z0 = self.get_cartesian_coords(i, j, k)

        # Special cases
        if x0 == 0.0:
            return 0.0, 0.0, 1.0
        if z0 == 0.0:
            return 1.0, 0.0, 0.0

        k = (self.a**2) * z0 / ((self.b[i]**2) * x0)
        normal_x = math.copysign(1.0, k) / math.sqrt(1 + k**2)
        normal_y = 0.0
        normal_z = k * normal_x
        return normal_x, normal_y, normal_z


if __name__ == '__main__':
    # test case
    presCylField = PressureFieldCyl()
    presCylField.set_nodes_num(1, 101, 101)
    presCylField.set_grid_bottom(3.07e-03, -0.3*3.14, -2.0e-03)
    presCylField.set_grid_top(3.07e-03, 0.3*3.14, 9.0e-03)
    presCylField.prepare_grid()
