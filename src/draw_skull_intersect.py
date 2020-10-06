# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import output
import numpy
import codecs

import math


def main():
    draw_skull()
    # propagate_skull()


def propagate_skull():
    import analytic_calc_skull
    x, y, z3d_inside, z3d_outside = analytic_calc_skull.load_skull_info_from_matlab()

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

    # len_of_sources = numpy.shape(sources['Xs'])[0]
    # len_of_sources = numpy.shape(sources['Ys'])[0]

    x0 = 0
    y0 = 0
    z0 = 0

    xf = 50.0e-03
    yf = 10.0e-03
    zf = 130.0e-03

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
    cross_product2 = cross_product**2
    distances = numpy.sqrt(cross_product2[0,:,:] + cross_product2[1,:,:] + cross_product2[2,:,:])

    # Индексы минимумов
    ids = numpy.unravel_index(distances.argmin(), distances.shape)
    i_min_outside = ids[0]
    j_min_outside = ids[1]
    print(ids)

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
    print(ids)

    # i_min_array_inside = d_inside.argmin(0)
    # j_min_inside = (d_inside.argmin(0)).argmin(0)
    # i_min_inside = i_min_array_inside[j_min_inside]

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

    bone_distance = numpy.sqrt((x_cross_inside - x_cross_outside) ** 2 + (y_cross_inside - y_cross_outside) ** 2 + (z_cross_inside - z_cross_outside) ** 2)
    

    inside_coords = [x_cross_inside, y_cross_inside, z_cross_inside]
    outside_coords = [x_cross_outside, y_cross_outside, z_cross_outside]
    return inside_coords, outside_coords
    
    # bone_distance = numpy.sqrt((x_cross_inside - x_cross_outside) ** 2 + (y_cross_inside - y_cross_outside) ** 2 + (z_cross_inside - z_cross_outside) ** 2)

    # return numpy.exp(-alpha * bone_distance) * numpy.exp(1j * dk * bone_distance) * T12
    # return numpy.exp(1j * dk * bone_distance)
    # return numpy.exp(-alpha_at * bone_distance) * T12


def draw_skull():
    import analytic_calc_skull
    x, y, z3d_inside, z3d_outside = analytic_calc_skull.load_skull_info_from_matlab()
    inside_coords, outside_coords = propagate_skull()

    x0 = 0
    y0 = 0
    z0 = 0

    xf = 50.0e-03
    yf = 10.0e-03
    zf = 130.0e-03

    alpha = (xf-x0)/numpy.sqrt((xf-x0)**2+(yf-y0)**2+(zf-z0)**2)
    betha = (yf-y0)/numpy.sqrt((xf-x0)**2+(yf-y0)**2+(zf-z0)**2)
    gamma = (zf-z0)/numpy.sqrt((xf-x0)**2+(yf-y0)**2+(zf-z0)**2)
    s = numpy.array([alpha, betha, gamma])  # Направляющий вектор
    # print(numpy.shape(s))
    t = numpy.linspace(0, numpy.sqrt((xf-x0)**2+(yf-y0)**2+(zf-z0)**2),100)  # Параметр

    # %Уравнение прямой%
    x_line = x0+alpha*t
    y_line = y0+betha*t
    z_line = z0+gamma*t

    # Set Font
    plt.rc('font', **{'family':'serif','serif':['Times New Roman'],'size':18})
    # Set figure size and dpi
    fig = plt.figure(figsize=(6.0, 8.0), dpi=96, facecolor='w', edgecolor='k')

    # plt.subplot(121)
    # plt.contourf(X, Y, Z, 100, cmap=plt.cm.jet)
    # plt.plot(x, z_inside, '-w', lw=1.5)
    # plt.plot(x, z_outside, '-w', lw=1.5)
    # plt.colorbar(orientation = 'vertical')

    from mpl_toolkits.mplot3d import Axes3D
    # fig = plt.figure()
    ax = fig.gca(projection='3d')


    # matr = numpy.meshgrid(x,y,z3d_inside)
    # x,y,z3d_outside = numpy.meshgrid(x,y,z3d_outside)

    
    ax.plot_surface(x,y,z3d_inside)
    # ax.plot_wireframe(x,y,z3d_inside)
    # ax.plot_wireframe(x,y,z3d_inside)
    # ax.plot_trisurf(x,y,z3d_inside, 'g')
    ax.plot_surface(x,y,z3d_outside, 'b')
    ax.plot(x_line,y_line,z_line, 'r')

    ax.scatter(inside_coords[0], inside_coords[1], inside_coords[2], 'r')
    ax.scatter(outside_coords[0], outside_coords[1], outside_coords[2], 'r')
    # ax.plot(outside_coords, 'r')
    # plot3(x_cross_outside,y_cross_outside,z_cross_outside,'LineStyle','o','Color','r');
    # plot3(x_cross_inside,y_cross_inside,z_cross_inside,'LineStyle','o','Color','r');
    # axis('equal');

    # plt.xlabel(r'$\it{x}$, mm')
    # plt.ylabel(r'$\it{z}$, mm')

    # plt.ylim(0, 50)
    # plt.clim(0, 1.80)

    # plt.title('Amplitude')
    # plt.gca().set_aspect('equal')

    # plt.subplot(122)
    # plt.contourf(X,Y,Z2,20,cmap=plt.cm.jet)
    # plt.title('Phase')
    # plt.xlabel(r'$\it{x}$, mm')
    # # plt.xlim(-2.5, 7.5)
    # plt.gca().set_aspect('equal')

    # plt.suptitle("Distance from trans = 10 mm")

    # plt.savefig(r"fig_XZ.png")
    plt.show()


if __name__ == '__main__':
    main()
