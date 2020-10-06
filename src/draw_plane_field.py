# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import output
import numpy
import codecs


def main():
    # draw_XY_sveta_with_ribs(r"c:\Downloads\Calc\Sveta_res\intens.txt")
    # field = output.restore_field_from_disk(r"c:\Downloads\Calc\2014-06-30_00-06-31 good\field_plane")
    # draw_Y_compare_sveta(r"c:\Downloads\Calc\Sveta_res\intens.txt", field)

    # field = output.restore_field_from_disk(r"d:\Downloads\Calc Skull\2015-06-15_15-09-40\field_plane")
    # draw_X(field)
    # field1 = output.restore_field_from_disk(r"d:\Downloads\Calc Skull\2015-06-15_15-01-23 x orig\field_plane")
    # field2 = output.restore_field_from_disk(r"d:\Downloads\Calc Skull\2015-06-15_15-15-35 x all\field_plane")
    # draw_X_compare(field1, field2)
    field = output.restore_field_from_disk(r"d:\Downloads\Calc Skull\2015-06-15_16-50-12\field_plane")
    # Z1 = numpy.absolute(field.p[:, 0, :])
    # Z2 = numpy.absolute(field2.p[:, 0, :])
    # print(Z1-Z2)

    # draw_XZ(field)
    draw_XZ_with_skull(field)
    # field = output.restore_field_from_disk(r"d:\Downloads\Calc\2014-09-15_16-53-47 one el 4.95\field_plane")
    # field = output.restore_field_from_disk(r"d:\Downloads\Calc\2014-09-15_16-48-12 one el 3.5\field_plane")
    # field = output.restore_field_from_disk(r"d:\Downloads\Calc\2014-09-17_17-51-10\field_plane")
    # draw_XY_with_ribs(field)

    # field = output.restore_field_from_disk(r"d:\Downloads\Calc\2014-09-22_16-17-08\field_plane")
    # field = output.restore_field_from_disk(r"d:\Downloads\Calc\2014-10-01_09-54-02\field_plane")
    # field = output.restore_field_from_disk(r"d:\Downloads\Calc\2014-10-01_09-48-05\field_plane")
    # field = output.restore_field_from_disk(r"d:\Downloads\Calc\field_plane_ribs_clean")
    # import ribs
    # ribs.CalcPowerOnPlane(field)
    # import rayleigh
    # rayleigh.calc_power_in_main_focus(field, 4.0e-03)
    # draw_XY(field)

    # field1 = output.restore_field_from_disk(r"d:\Downloads\Calc Skull\2015-05-12_11-06-00\field_plane")
    # field2 = output.restore_field_from_disk(r"d:\Downloads\Calc Skull\2015-05-12_11-07-50\field_plane")
    # draw_X_compare(field1, field2)

    # field1 = output.restore_field_from_disk(r"c:\Downloads\Calc\Comp_x\anal2\field_plane")
    # field2 = output.restore_field_from_disk(r"c:\Downloads\Calc\Comp_x\ray\field_plane")
    # field2 = output.restore_field_from_disk(r"c:\Downloads\Calc\2014-06-26_14-44-09\field_plane")
    # field2 = output.restore_field_from_disk(r"c:\Downloads\Calc\field_plane")
    # draw_X_compare(field1, field2)
    # field2 = output.restore_field_from_disk(r"c:\Downloads\Calc\Comp_z\ray\field_plane")
    # field1 = output.restore_field_from_disk(r"c:\Downloads\Calc\Comp_z\anal\field_plane")
    # draw_Z_compare(field1, field2)


    # field = output.restore_field_from_disk(r"c:\Downloads\Calc\field_plane")
    # draw_XY_with_ribs(field)

    # draw_Y(field)
    # draw_PhiZ(field)

    # if field.n_z == 1:
    #   draw_XY_with_ribs(field)
    #   return

    # if field.n_y == 1:
    #   draw_XZ(field)
    #   return

    # if field.n_x == 1:
    #   draw_YZ(field)
    #   return

    # field1 = output.restore_field_from_disk(ur'd:\Educ\АКУСТИКА\Philips R&D\AcField Cylindrical shell\Py\2012-11-08_11-41-00 Direct Z\field')
    # field2 = output.restore_field_from_disk(ur'd:\Educ\АКУСТИКА\Philips R&D\AcField Cylindrical shell\Py\2012-11-08_12-35-10\field')

    # compare_Z(field1, field2)

def draw_XY_with_ribs(field):
    Y, X = numpy.meshgrid(field.y, field.x)
    X *= 1.0e03
    Y *= 1.0e03
    # Z = numpy.absolute(field.p[:, :, 0]) # * numpy.absolute(field.p[:, :, 0])
    Z =  0.5 * numpy.real(field.p[:,:,0] * numpy.conj(field.vn[:,:,0]*1000.0 * 1500.0))

    # Z *= 1000.0 * 1500.0
    print("Full Energy = ", numpy.sum(Z)*0.1*0.1*1.0e-06)
    print(numpy.max(Z))

    #-----------------------------------

    # Set Font
    plt.rc('font', **{'family':'serif','serif':['Times New Roman'],'size':12})
    # Set figure size and dpi
    plt.figure(figsize=(6.0, 6.0), dpi=96, facecolor='w', edgecolor='k')

    # plt.subplot(121)
    cs = plt.contourf(X, Y, Z, 100, cmap=plt.cm.jet)
    # plt.colorbar(orientation = 'vertical')
    # plt.xlabel('x, mm')
    plt.xlabel(r'$\it{x}$, mm')
    plt.ylabel(r'$\it{y}$, mm')
    # plt.ylabel('y, mm')
    plt.locator_params(nbins=6)
    cb = plt.colorbar(cs, orientation='vertical', format='%.1f')
    cb.set_label(r'$\it{p/p_{0}}$')
    cb.locator = plt.MaxNLocator(nbins=6)
    cb.update_ticks()
    plt.gca().set_aspect('equal')

    rib_x = numpy.arange(-100,100,1)
    rib1_top_y = 73.0
    rib1_bot_y = 55.0
    rib2_top_y = 41.0
    rib2_bot_y = rib2_top_y - 18.0
    rib3_top_y = rib2_bot_y - 14.0
    rib3_bot_y = rib3_top_y - 18.0
    rib4_top_y = rib3_bot_y - 14.0
    rib4_bot_y = rib4_top_y - 18.0
    rib5_top_y = rib4_bot_y - 14.0
    rib5_bot_y = rib5_top_y - 18.0
    plt.plot(rib_x, rib1_top_y*numpy.ones(200), '-.w')
    plt.plot(rib_x, rib1_bot_y*numpy.ones(200), '-.w')
    plt.plot(rib_x, rib2_top_y*numpy.ones(200), '-.w')
    plt.plot(rib_x, rib2_bot_y*numpy.ones(200), '-.w')
    plt.plot(rib_x, rib3_top_y*numpy.ones(200), '-.w')
    plt.plot(rib_x, rib3_bot_y*numpy.ones(200), '-.w')
    plt.plot(rib_x, rib4_top_y*numpy.ones(200), '-.w')
    plt.plot(rib_x, rib4_bot_y*numpy.ones(200), '-.w')
    plt.plot(rib_x, rib5_top_y*numpy.ones(200), '-.w')
    plt.plot(rib_x, rib5_bot_y*numpy.ones(200), '-.w')

    plt.savefig(r"d:\Downloads\Calc\fig_XY.png")
    plt.show()

def draw_XY_sveta_with_ribs(sveta_file):
    with open(sveta_file, 'r') as f:
        Z = numpy.loadtxt(f)

    Z = numpy.transpose(Z)

    print("Full Energy = ", numpy.sum(Z)*0.1*0.1*1.0e-06)

    print(numpy.max(Z))

    x = numpy.linspace(-100, 100, num=2001)
    y = x
    Y, X = numpy.meshgrid(x, y)
    # X *= 1.0e03
    # Y *= 1.0e03
    # Z = numpy.absolute(field.p[:, :, 0]) # * numpy.absolute(field.p[:, :, 0])
    # Z2 = numpy.angle(field.p[:, :, 0])

    # #-----------------------------------

    # Set Font
    plt.rc('font', **{'family':'serif','serif':['Times New Roman'],'size':12})
    # Set figure size and dpi
    plt.figure(figsize=(6.0, 6.0), dpi=96, facecolor='w', edgecolor='k')

    # plt.subplot(121)
    cs = plt.contourf(X, Y, Z, 100, cmap=plt.cm.jet)
    # plt.colorbar(orientation = 'vertical')
    # plt.xlabel('x, mm')
    plt.xlabel(r'$\it{x}$, mm')
    plt.ylabel(r'$\it{y}$, mm')
    # plt.ylabel('y, mm')
    plt.locator_params(nbins=6)
    cb = plt.colorbar(cs, orientation='vertical', format='%.1f')
    cb.set_label(r'$\it{p/p_{0}}$')
    cb.locator = plt.MaxNLocator(nbins=6)
    cb.update_ticks()
    plt.gca().set_aspect('equal')

    rib_x = numpy.arange(-100,100,1)
    rib1_top_y = 73.0
    rib1_bot_y = 55.0
    rib2_top_y = 41.0
    rib2_bot_y = rib2_top_y - 18.0
    rib3_top_y = rib2_bot_y - 14.0
    rib3_bot_y = rib3_top_y - 18.0
    rib4_top_y = rib3_bot_y - 14.0
    rib4_bot_y = rib4_top_y - 18.0
    rib5_top_y = rib4_bot_y - 14.0
    rib5_bot_y = rib5_top_y - 18.0
    plt.plot(rib_x, rib1_top_y*numpy.ones(200), '-.w')
    plt.plot(rib_x, rib1_bot_y*numpy.ones(200), '-.w')
    plt.plot(rib_x, rib2_top_y*numpy.ones(200), '-.w')
    plt.plot(rib_x, rib2_bot_y*numpy.ones(200), '-.w')
    plt.plot(rib_x, rib3_top_y*numpy.ones(200), '-.w')
    plt.plot(rib_x, rib3_bot_y*numpy.ones(200), '-.w')
    plt.plot(rib_x, rib4_top_y*numpy.ones(200), '-.w')
    plt.plot(rib_x, rib4_bot_y*numpy.ones(200), '-.w')
    plt.plot(rib_x, rib5_top_y*numpy.ones(200), '-.w')
    plt.plot(rib_x, rib5_bot_y*numpy.ones(200), '-.w')

    # plt.savefig(r"fig_XY.png")
    plt.show()

def draw_Y_compare_sveta(sveta_file, field):
    with open(sveta_file, 'r') as f:
        Z = numpy.loadtxt(f)

    Z = numpy.transpose(Z)
    Z1 = Z[1001,:]
    Y1 = numpy.linspace(-100, 100, num=2001)

    Z2 =  0.5 * numpy.real(field.p[1001,:,0] * numpy.conj(field.vn[1001,:,0]*1000.0 * 1500.0))
    #-----------------------------------

    # Set Font
    plt.rc('font', **{'family':'serif','serif':['Times New Roman'],'size':12})
    # Set figure size and dpi
    plt.figure(figsize=(11.0, 4.0), dpi=96, facecolor='w', edgecolor='k')

    # plt.subplot(121)
    plt.plot(Y1, Z1, 'k')
    plt.plot(Y1, Z2, 'r')
    plt.xlabel(r'$\it{y}$, mm')
    plt.axis([-100,100,0,3])
    plt.legend(['Sveta', 'Sergey'])

    # plt.savefig(r"fig_XY.png")
    plt.show()

def draw_XY(field):
    # X, Y = numpy.meshgrid(field.y, field.x)
    Y, X = numpy.meshgrid(field.y, field.x)
    X *= 1.0e03
    Y *= 1.0e03
    Z = numpy.absolute(field.p[:, :, 0]) * numpy.absolute(field.p[:, :, 0])
    Z2 = numpy.angle(field.p[:, :, 0])

    #-----------------------------------

    # print(numpy.max(Z))
    print(numpy.sum(Z)*field.one_pixel_area)

    # print numpy.shape(X)

    # Set Font
    plt.rc('font', **{'family':'serif','serif':['Times New Roman'],'size':12})
    # Set figure size and dpi
    plt.figure(figsize=(8.0, 8.0), dpi=96, facecolor='w', edgecolor='k')

    # plt.subplot(121)
    cs = plt.contourf(X, Y, Z, 100, cmap=plt.cm.jet)
    # plt.colorbar(orientation = 'vertical')
    # plt.xlabel('x, mm')
    plt.xlabel(r'$\it{x}$, mm')
    plt.ylabel(r'$\it{y}$, mm')
    # plt.ylabel('y, mm')
    plt.locator_params(nbins=6)
    cb = plt.colorbar(cs, orientation='vertical', format='%.1f')
    cb.set_label(r'$\it{p/p_{0}}$')
    cb.locator = plt.MaxNLocator(nbins=6)
    cb.update_ticks()
    plt.gca().set_aspect('equal')

    # plt.plot(X[:,0], Y[:, 12], 'w--')
    # plt.plot(X[:,0], Y[:, 30], 'w--')
    # plt.plot(X[:,0], Y[:, 44], 'w--')
    # plt.plot(X[:,0], Y[:, 62], 'w--')
    # plt.plot(X[:,0], Y[:, 76], 'w--')
    # plt.plot(X[:,0], Y[:, 94], 'w--')
    # plt.plot(X[:,0], Y[:, 108], 'w--')
    # plt.plot(X[:,0], Y[:, 126], 'w--')
    # plt.plot(X[:,0], Y[:, 140], 'w--')
    # plt.plot(X[:,0], Y[:, 158], 'w--')
    # plt.clim(37.0, 50.0)
    # plt.title('Amplitude')

    # plt.subplot(122)
    # cs = plt.contourf(Y, X, Z2, 30, cmap=plt.cm.jet)
    # plt.locator_params(nbins=6)
    # cb = plt.colorbar(orientation='horizontal', format='%.2f')
    # cb.set_label('radians')
    # cb.locator = plt.MaxNLocator(nbins=6)
    # cb.update_ticks()
    # plt.gca().set_aspect('equal')
    # plt.xlabel('x, mm')
    # # plt.xlabel(r'$\it{y}$, mm')
    # plt.yticks([])
    # # plt.ylabel('points')
    # plt.title('Phase')

    # plt.suptitle("Distance from trans = 10 mm")

    plt.savefig(r"d:\Downloads\Calc\fig_XY.png")
    plt.show()


def draw_XZ(field):
    Y, X = numpy.meshgrid(field.z, field.x)
    X *= 1.0e03
    Y *= 1.0e03
    Z = numpy.absolute(field.p[:, 0, :])
    # Z2 = numpy.angle(field.p[:, 0, :])

    #-----------------------------------

    # print(numpy.max(Z))

    # Set Font
    plt.rc('font', **{'family':'serif','serif':['Times New Roman'],'size':18})
    # Set figure size and dpi
    plt.figure(figsize=(6.0, 8.0), dpi=96, facecolor='w', edgecolor='k')

    # plt.subplot(121)
    plt.contourf(X, Y, Z, 100, cmap=plt.cm.jet)
    plt.colorbar(orientation = 'vertical')
    plt.xlabel(r'$\it{x}$, mm')
    plt.ylabel(r'$\it{z}$, mm')

    # plt.ylim(0, 50)
    # plt.clim(0, 1.80)

    # plt.title('Amplitude')
    plt.gca().set_aspect('equal')

    # plt.subplot(122)
    # plt.contourf(X,Y,Z2,20,cmap=plt.cm.jet)
    # plt.title('Phase')
    # plt.xlabel(r'$\it{x}$, mm')
    # # plt.xlim(-2.5, 7.5)
    # plt.gca().set_aspect('equal')

    # plt.suptitle("Distance from trans = 10 mm")

    # plt.savefig(r"fig_XZ.png")
    plt.show()


def draw_YZ(field):
    Y, X = numpy.meshgrid(field.z, field.y)
    # X, Y = numpy.meshgrid(field.y, field.z)
    X *= 1.0e03
    Y *= 1.0e03
    Z = numpy.absolute(field.p[0, :, :])
    # Z2 = numpy.angle(field.p[:, 0, :])

    #-----------------------------------

    # Set Font
    plt.rc('font', **{'family':'serif','serif':['Times New Roman'],'size':18})
    # Set figure size and dpi
    plt.figure(figsize=(10.0, 8.0), dpi=96, facecolor='w', edgecolor='k')

    # plt.subplot(121)
    plt.contourf(Y, X, Z, 100, cmap=plt.cm.jet)
    plt.colorbar(orientation = 'vertical')
    plt.xlabel(r'$\it{z}$, mm')
    plt.ylabel(r'$\it{y}$, mm')
    # plt.ylim(0, 50)
    # plt.title('Amplitude')
    plt.gca().set_aspect('equal')

    # plt.savefig(ur"d:\Educ\АКУСТИКА\Philips R&D\AcField Cylindrical shell\Py\fig_YZ.png")
    plt.show()

def draw_YZ_contour(field):
    Y, X = numpy.meshgrid(field.z, field.y)
    # X, Y = numpy.meshgrid(field.y, field.z)
    X *= 1.0e03
    Y *= 1.0e03
    Z = numpy.absolute(field.p[0, :, :])

    is1, is2 = find_isoline(Z)
    m = numpy.max(Z)
    l = m*0.5

    y1_s = 1.0e03*numpy.take(field.y, is1)
    y2_s = 1.0e03*numpy.take(field.y, is2)

    #-----------------------------------

    # Set Font
    plt.rc('font', **{'family':'serif','serif':['Times New Roman'],'size':18})
    # Set figure size and dpi
    plt.figure(figsize=(10.0, 4.0), dpi=96, facecolor='w', edgecolor='k')

    # plt.subplot(121)
    plt.contourf(Y, X, Z, 100, cmap=plt.cm.jet)
    # plt.set_dashes(cs)
    plt.colorbar(orientation = 'vertical')
    plt.xlabel(r'$\it{z}$, mm')
    plt.ylabel(r'$\it{y}$, mm')
    # plt.ylim(0, 50)
    # plt.title('Amplitude')
    plt.contour(Y, X, Z, 1, levels=[0.0, l])
    plt.gca().set_aspect('equal')

    plane_y = numpy.arange(-20,20,1)
    plane_z = 45.0
    plt.plot(plane_z*numpy.ones(40), plane_y, '--w')
    plt.plot(field.z*1000, y1_s, '--w', lw=2)
    plt.plot(field.z*1000, y2_s, '--w', lw=2)


    plt.savefig(r"d:\Downloads\Calc\fig_one.png")
    plt.show()

def find_isoline(z):
    # На каждой линии по Х
    # print(numpy.shape(z))
    (x_len, y_len) = numpy.shape(z)
    # isoline = numpy.zeros(101, 101, 2)
    isoline_1 = []
    isoline_2 = []
    for i in range(x_len):
        ar = z[:, i]
        m = numpy.max(ar)
        l = m*0.5
        m_x = numpy.argmax(ar)
        if len(numpy.abs(ar[:m_x]-l)) != 0:
            l1 = numpy.abs(ar[:m_x]-l).argmin()
        if len(numpy.abs(ar[m_x:]-l)) != 0:
            l2 = numpy.abs(ar[m_x:]-l).argmin()
        isoline_1.append(l1)
        isoline_2.append(l2+m_x)
    return isoline_1, isoline_2

def draw_PhiZ(field):
    Y, X = numpy.meshgrid(field.z, field.phi)
    # X *= 1.0e03
    Y *= 1.0e03
    Z = numpy.absolute(field.p[0, :, :])
    # Z2 = numpy.angle(field.p[:, 0, :])

    #-----------------------------------

    # Set Font
    plt.rc('font', **{'family':'serif','serif':['Times New Roman'],'size':18})
    # Set figure size and dpi
    plt.figure(figsize=(8.0, 7.0), dpi=96, facecolor='w', edgecolor='k')

    # plt.subplot(121)
    plt.contourf(X, Y, Z, 100, cmap=plt.cm.jet)
    plt.colorbar(orientation = 'vertical')
    plt.xlabel(r'$\it{phi}$, rad')
    plt.ylabel(r'$\it{z}$, mm')
    # plt.title('Amplitude')
    # plt.gca().set_aspect('equal')

    # plt.subplot(122)
    # plt.contourf(X,Y,Z2,20,cmap=plt.cm.jet)
    # plt.title('Phase')
    # plt.xlabel(r'$\it{x}$, mm')
    # # plt.xlim(-2.5, 7.5)
    # plt.gca().set_aspect('equal')

    # plt.suptitle("Distance from trans = 10 mm")

    plt.savefig(r"fig_PhiZ.png")
    plt.show()


def compare_Z(field1, field2):
    X1 = field1.z * 1.0e03
    Z1 = numpy.absolute(field1.p[0, 0, :])
    X2 = field2.z * 1.0e03
    Z2 = numpy.absolute(field2.p[0, 0, :])
    # Z2 = numpy.angle(field.p[:, 0, :])

    #-----------------------------------

    # Set Font
    # plt.rc('font', **{'family':'serif','serif':['Times New Roman'],'size':18})
    # Set figure size and dpi
    plt.figure(figsize=(10.0, 4.0), dpi=96, facecolor='w', edgecolor='k')

    # plt.subplot(121)
    plt.plot(X1, Z1, 'k')
    plt.plot(X2, Z2, 'r')
    plt.xlabel(r'$\it{z}$, mm')
    # plt.ylabel(r'$\it{z}$, mm')
    # plt.title('Amplitude')
    # plt.gca().set_aspect('equal')

    # plt.subplot(122)
    # plt.contourf(X,Y,Z2,20,cmap=plt.cm.jet)
    # plt.title('Phase')
    # plt.xlabel(r'$\it{x}$, mm')
    # # plt.xlim(-2.5, 7.5)
    # plt.gca().set_aspect('equal')

    # plt.suptitle("Distance from trans = 10 mm")

    plt.savefig(r"fig__comp_Z.png")
    plt.show()

def draw_Y(field):
    Y = field.y * 1.0e03
    Z = numpy.absolute(field.p[0, :, 0])
    #-----------------------------------

    # Set Font
    # plt.rc('font', **{'family':'serif','serif':['Times New Roman'],'size':18})
    # Set figure size and dpi
    plt.figure(figsize=(10.0, 4.0), dpi=96, facecolor='w', edgecolor='k')

    # plt.subplot(121)
    plt.plot(Y, Z, 'k')
    plt.xlabel(r'$\it{z}$, mm')
    # plt.ylabel(r'$\it{z}$, mm')
    # plt.title('Amplitude')
    # plt.gca().set_aspect('equal')

    # plt.subplot(122)
    # plt.contourf(X,Y,Z2,20,cmap=plt.cm.jet)
    # plt.title('Phase')
    # plt.xlabel(r'$\it{x}$, mm')
    # # plt.xlim(-2.5, 7.5)
    # plt.gca().set_aspect('equal')

    # plt.suptitle("Distance from trans = 10 mm")

    plt.savefig(r"e:\Downloads\New_calc\fig__Y.png")
    plt.show()

def draw_X(field):
    Y = field.x * 1.0e03
    Z = numpy.absolute(field.p[:, 0, 0])
    #-----------------------------------

    # Set Font
    # plt.rc('font', **{'family':'serif','serif':['Times New Roman'],'size':18})
    # Set figure size and dpi
    plt.figure(figsize=(10.0, 4.0), dpi=96, facecolor='w', edgecolor='k')

    # plt.subplot(121)
    plt.plot(Y, Z, 'k')
    plt.xlabel(r'$\it{x}$, mm')
    # plt.ylabel(r'$\it{z}$, mm')
    # plt.title('Amplitude')
    # plt.gca().set_aspect('equal')

    # plt.subplot(122)
    # plt.contourf(X,Y,Z2,20,cmap=plt.cm.jet)
    # plt.title('Phase')
    # plt.xlabel(r'$\it{x}$, mm')
    # # plt.xlim(-2.5, 7.5)
    # plt.gca().set_aspect('equal')

    # plt.suptitle("Distance from trans = 10 mm")

    # plt.savefig(r"e:\Downloads\New_calc\fig__Y.png")
    plt.show()

def draw_Y_compare(field1, field2):
    Y1 = field1.y * 1.0e03
    Z1 = numpy.absolute(field1.p[0, :, 0])

    Y2 = field2.y * 1.0e03
    Z2 = numpy.absolute(field2.p[0, :, 0])
    #-----------------------------------

    # Set Font
    # plt.rc('font', **{'family':'serif','serif':['Times New Roman'],'size':18})
    # Set figure size and dpi
    plt.figure(figsize=(10.0, 4.0), dpi=96, facecolor='w', edgecolor='k')

    # plt.subplot(121)
    plt.plot(Y1, Z1, '--k')
    plt.plot(Y2, Z2, 'k')
    plt.xlabel(r'$\it{z}$, mm')
    # plt.ylabel(r'$\it{z}$, mm')
    # plt.title('Amplitude')
    # plt.gca().set_aspect('equal')

    # plt.subplot(122)
    # plt.contourf(X,Y,Z2,20,cmap=plt.cm.jet)
    # plt.title('Phase')
    # plt.xlabel(r'$\it{x}$, mm')
    # # plt.xlim(-2.5, 7.5)
    # plt.gca().set_aspect('equal')

    # plt.suptitle("Distance from trans = 10 mm")


    rib_x = numpy.arange(0,2,0.01)
    rib1_top_y = 73.0
    rib1_bot_y = 55.0
    rib2_top_y = 41.0
    rib2_bot_y = rib2_top_y - 18.0
    rib3_top_y = rib2_bot_y - 14.0
    rib3_bot_y = rib3_top_y - 18.0
    rib4_top_y = rib3_bot_y - 14.0
    rib4_bot_y = rib4_top_y - 18.0
    rib5_top_y = rib4_bot_y - 14.0
    rib5_bot_y = rib5_top_y - 18.0
    plt.plot(rib1_top_y*numpy.ones(200), rib_x, '-.w')
    plt.plot(rib1_bot_y*numpy.ones(200), rib_x, '-.w')
    plt.plot(rib_x, rib2_top_y*numpy.ones(200), rib_x, '-.w')
    plt.plot(rib_x, rib2_bot_y*numpy.ones(200), rib_x, '-.w')
    plt.plot(rib_x, rib3_top_y*numpy.ones(200), '-.w')
    plt.plot(rib_x, rib3_bot_y*numpy.ones(200), '-.w')
    plt.plot(rib_x, rib4_top_y*numpy.ones(200), '-.w')
    plt.plot(rib_x, rib4_bot_y*numpy.ones(200), '-.w')
    plt.plot(rib_x, rib5_top_y*numpy.ones(200), '-.w')
    plt.plot(rib_x, rib5_bot_y*numpy.ones(200), '-.w')



    plt.savefig(r"e:\Downloads\New_calc\fig__Y.png")
    plt.show()


def draw_X_compare(field1, field2):
    Y1 = field1.x * 1.0e03
    Z1 = numpy.absolute(field1.p[:, 0, -1])

    Y2 = field2.x * 1.0e03
    Z2 = numpy.absolute(field2.p[:, 0, -1])
    #-----------------------------------

    # Set Font
    # plt.rc('font', **{'family':'serif','serif':['Times New Roman'],'size':18})
    # Set figure size and dpi
    plt.figure(figsize=(10.0, 4.0), dpi=96, facecolor='w', edgecolor='k')

    # plt.subplot(121)
    plt.plot(Y1, Z1, 'k')
    plt.plot(Y2, Z2, 'r')
    plt.xlabel(r'$\it{x}$, mm')
    # plt.legend(['Anal', 'Num'])
    plt.legend(['No corrections', 'All corrections'])
    plt.title('Distribution in focal plane z = 130 mm')
    plt.show()

def draw_Z(field):
    X1 = field.z * 1.0e03
    Z1 = numpy.absolute(field.p[0, 0, :])
    #-----------------------------------

    # Set Font
    plt.rc('font', **{'family':'serif','serif':['Times New Roman'],'size':12})
    
    # Set figure size and dpi
    plt.figure(figsize=(10.0, 4.0), dpi=96, facecolor='w', edgecolor='k')

    plt.plot(X1, Z1, 'k')
    plt.xlabel(r'$\it{z}$, mm')
    # plt.legend(['Anal', 'Num'])

    plt.show()


def draw_Z_compare(field1, field2):
    X1 = field1.z * 1.0e03
    Z1 = numpy.absolute(field1.p[0, 0, :])
    V1 = numpy.absolute(field1.vn[0, 0, :])

    X2 = field2.z * 1.0e03
    Z2 = numpy.absolute(field2.p[0, 0, :])
    V2 = numpy.absolute(field2.vn[0, 0, :])
    #-----------------------------------

    # Set Font
    plt.rc('font', **{'family':'serif','serif':['Times New Roman'],'size':12})
    
    # Set figure size and dpi
    plt.figure(figsize=(10.0, 4.0), dpi=96, facecolor='w', edgecolor='k')

    plt.plot(X1, Z1, 'k')
    plt.plot(X2, Z2, 'r')
    plt.xlabel(r'$\it{z}$, mm')
    plt.legend(['Orig', 'Phase'])

    # # Set figure size and dpi
    # plt.figure(figsize=(10.0, 4.0), dpi=96, facecolor='w', edgecolor='k')

    # plt.plot(X1, V1, 'k')
    # plt.plot(X2, V2, 'r')
    # plt.xlabel(r'$\it{z}$, mm')
    # plt.legend(['Anal', 'Num'])
    
    plt.show()


def draw_XZ_with_skull(field):
    # Adding new medium from file
    # 'utf-8-sig' means that it skips BOM in UTF files created in Notepad
    file_path1 = r"d:\yandex_disk\SublimeProjects\SkullWay\test_files\x.txt"
    file_path2 = r"d:\yandex_disk\SublimeProjects\SkullWay\test_files\z_inside.txt"
    file_path3 = r"d:\yandex_disk\SublimeProjects\SkullWay\test_files\z_outside.txt"
    with codecs.open(file_path1, encoding='utf-8-sig') as f:
        x = [float(line.strip()) for line in f]

    with codecs.open(file_path2, encoding='utf-8-sig') as f:
        z_inside = [float(line.strip()) for line in f]

    with codecs.open(file_path3, encoding='utf-8-sig') as f:
        z_outside = [float(line.strip()) for line in f]

    # DX = field_x - sources['Xs']
    # DZ = field_z - sources['Zs']
    # k = DZ/DX
    # b = field_z - k*field_x
    x = numpy.array(x)
    z_inside = numpy.array(z_inside)
    z_outside = numpy.array(z_outside)
    k = -2
    b = 0
    z_ray = k*x + b
    dz_out = z_ray - z_outside
    dz_in = z_ray - z_inside

    d_out = numpy.nonzero(numpy.diff(numpy.sign(dz_out)))
    d_in = numpy.nonzero(numpy.diff(numpy.sign(dz_in)))

    if ((len(d_out) > 1) or (len(d_in) > 1)):
        print('Out of brain')
    else:
        d_out = d_out[0]
        d_in = d_in[0]
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


    Y, X = numpy.meshgrid(field.z, field.x)
    X *= 1.0e03
    Y *= 1.0e03
    Z = numpy.absolute(field.p[:, 0, :])
    # Z2 = numpy.angle(field.p[:, 0, :])

    #-----------------------------------

    # print(numpy.max(Z))

    # Set Font
    plt.rc('font', **{'family':'serif','serif':['Times New Roman'],'size':18})
    # Set figure size and dpi
    plt.figure(figsize=(6.0, 8.0), dpi=96, facecolor='w', edgecolor='k')

    # plt.subplot(121)

    import matplotlib as mpl

    norm = mpl.colors.Normalize(vmin=0, vmax=0.405)

    plt.contourf(X, Y, Z, 500, cmap=plt.cm.jet, vmin=0, vmax = 0.405)
    plt.plot(x, z_inside, '-w', lw=1.5)
    plt.plot(x, z_outside, '-w', lw=1.5)
    cb = plt.colorbar(orientation = 'vertical')
    cb.set_clim(0, 0.405)
    # plt.vmin(0, 0.405)

    # plt.plot(X[:,-1], Z[:,-1], '-b', lw=1.0)
    # plt.plot(X, Z, , lw=1.0)

    plt.xlabel(r'$\it{x}$, mm')
    plt.ylabel(r'$\it{z}$, mm')

    # plt.ylim(50, 130)
    # plt.clim(0, 1.80)

    # plt.title('Amplitude')
    plt.gca().set_aspect('equal')

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
