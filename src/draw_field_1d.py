# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import output
import numpy

def main():
	# field1 = output.restore_field_from_disk(r"d:\Downloads\Calc\Фокусир излуч\X_axis\Number 300 on line\field_plane")
	# field2 = output.restore_field_from_disk(r"d:\Downloads\Calc\Фокусир излуч\X_axis\Reference\field_plane")
	# field1 = output.restore_field_from_disk(r"d:\Downloads\Calc\2014-12-25_14-11-26 500\field_plane")
	field1 = output.restore_field_from_disk(r"d:\Downloads\Calc\2014-12-25_14-45-00 через 1\field_plane")
	field2 = output.restore_field_from_disk(r"d:\Downloads\Calc\2014-12-25_14-12-51 1000\field_plane")
	field3 = output.restore_field_from_disk(r"d:\Downloads\Calc\2014-12-25_14-20-55 1200\field_plane")
	# draw_Y_compare(field1, field2, field3)
	# draw_Z(field1)
	# draw_X_compare(field1,field2)
	# draw_Z_compare(field1,field2)
	draw_X_compare(field1, field2, field3)


def main2():
	# field1 = output.restore_field_from_disk(r"d:\Downloads\Calc\Фокусир излуч\X_axis\Number 300 on line\field_plane")
	# field2 = output.restore_field_from_disk(r"d:\Downloads\Calc\Фокусир излуч\X_axis\Reference\field_plane")
	field1 = output.restore_field_from_disk(r"d:\Downloads\Calc\Фокусир излуч\Z _axis\Number 300 on line\field_plane")
	field2 = output.restore_field_from_disk(r"d:\Downloads\Calc\Фокусир излуч\Z _axis\Reference\field_plane")
	# draw_Y_compare(field1, field2, field3)
	# draw_Z(field1)
	# draw_X_compare(field1,field2)
	draw_Z_compare(field1,field2)

def main_1_el():
	field1 = output.restore_field_from_disk(r"c:\Downloads\Calc\2014-06-23_08-03-42 - Exact\field_Z")
	field2 = output.restore_field_from_disk(r"c:\Downloads\Calc\2014-06-23_08-03-15 - Anal\field_Z")
	draw_Z_compare(field1,field2)

def main_compare():
	# import matplotlib.font_manager

	# matplotlib.font_manager.FontProperties(fname=r'c:\Windows\Fonts\Times New Roman Cyr Regular.ttf')
	# print(matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext='ttf'))
	# print([i for i in matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext='ttf') if 'times' in i.lower()])
	field1 = output.restore_field_from_disk(r"E:\YandexDisk\Python test\New_calc_scan\array_256_num_Z_field\field_Z")
	field2 = output.restore_field_from_disk(r"E:\YandexDisk\Python test\New_calc_scan\array_256_Z_field\field_Z")
	field3 = output.restore_field_from_disk(r"E:\YandexDisk\Python test\New_calc_scan\array_256_reg_Z_field\field_Z")
	field4 = output.restore_field_from_disk(r"E:\YandexDisk\Python test\New_calc_scan\array_1024_Z_field\field_Z")
	# draw_Z(field)
	# draw_Z_compare(field1, field2)
	# draw_Z_compare3(field1, field2, field3)
	draw_Z_compare3_subplot(field1, field2, field3,field4)

	# draw_Y(field)
	# draw_PhiZ(field)

	# if field.n_z == 1:
	# 	draw_XY_with_ribs(field)
	# 	return

	# if field.n_y == 1:
	# 	draw_XZ(field)
	# 	return

	# if field.n_x == 1:
	# 	draw_YZ(field)
	# 	return

	# field1 = output.restore_field_from_disk(ur'd:\Educ\АКУСТИКА\Philips R&D\AcField Cylindrical shell\Py\2012-11-08_11-41-00 Direct Z\field')
	# field2 = output.restore_field_from_disk(ur'd:\Educ\АКУСТИКА\Philips R&D\AcField Cylindrical shell\Py\2012-11-08_12-35-10\field')

	# compare_Z(field1, field2)

def draw_X(field):
	X = field.x * 1.0e03
	Z = numpy.absolute(field.p[:, 0, 0]) * numpy.absolute(field.p[:, 0, 0])
	#-----------------------------------

	# Set Font
	# plt.rc('font', **{'family':'Times New Roman','size':16, 'weight':'bold', 'style':'normal'})
	# plt.rc('font', family='serif')
	plt.rc('font', serif='Times New Roman')
	plt.rc('font', family='serif')
	plt.rc('font', size='11')
	# plt.rc('text', usetex=True)
	# plt.rc('font', serif='Times')	
	# plt.rc('font', weight='bold')
	# plt.rc('text', usetex=True)
	# plt.rc('axes', labelweight='light')
	# plt.rc('xtick', labelsize=10)
	
	
	# plt.rc('font', family='Times New Roman') 
	# plt.rc('font', **{'family':'serif','serif':['Times New Roman'],'size':16, 'weight':20, 'style':'normal'})
	# Set figure size and dpi
	plt.figure(figsize=(5.0, 5.0), dpi=96, facecolor='w', edgecolor='k')

	# plt.subplot(121)
	plt.plot(X, Z, 'k')
	# plt.gca().axis([0,130,0,2])
	# plt.xlabel(r'$\it{z}$, mm')
	# plt.ylabel(r'$\it{p/p}_0$')
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


def draw_Y(field):
	X = field.y * 1.0e03
	Z = numpy.absolute(field.p[0, :, 0]) * numpy.absolute(field.p[0, :, 0])
	#-----------------------------------

	# Set Font
	# plt.rc('font', **{'family':'Times New Roman','size':16, 'weight':'bold', 'style':'normal'})
	# plt.rc('font', family='serif')
	plt.rc('font', serif='Times New Roman')
	plt.rc('font', family='serif')
	plt.rc('font', size='11')
	# plt.rc('text', usetex=True)
	# plt.rc('font', serif='Times')	
	# plt.rc('font', weight='bold')
	# plt.rc('text', usetex=True)
	# plt.rc('axes', labelweight='light')
	# plt.rc('xtick', labelsize=10)
	
	
	# plt.rc('font', family='Times New Roman') 
	# plt.rc('font', **{'family':'serif','serif':['Times New Roman'],'size':16, 'weight':20, 'style':'normal'})
	# Set figure size and dpi
	plt.figure(figsize=(5.0, 5.0), dpi=96, facecolor='w', edgecolor='k')

	# plt.subplot(121)
	plt.plot(X, Z, 'k')
	# plt.gca().axis([0,130,0,2])
	# plt.xlabel(r'$\it{z}$, mm')
	# plt.ylabel(r'$\it{p/p}_0$')
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


def draw_Y_compare(field1, field2, field3):
	X1 = field1.y * 1.0e03
	Z1 = numpy.absolute(field1.p[0, :, 0]) * numpy.absolute(field1.p[0, :, 0])
	X2 = field2.y * 1.0e03
	Z2 = numpy.absolute(field2.p[0, :, 0]) * numpy.absolute(field2.p[0, :, 0])
	X3 = field3.y * 1.0e03
	Z3 = numpy.absolute(field3.p[0, :, 0]) * numpy.absolute(field3.p[0, :, 0])
	
	#-----------------------------------

	# Set Font
	# plt.rc('font', **{'family':'Times New Roman','size':16, 'weight':'bold', 'style':'normal'})
	# plt.rc('font', family='serif')
	plt.rc('font', serif='Times New Roman')
	plt.rc('font', family='serif')
	plt.rc('font', size='11')
	# plt.rc('text', usetex=True)
	# plt.rc('font', serif='Times')	
	# plt.rc('font', weight='bold')
	# plt.rc('text', usetex=True)
	# plt.rc('axes', labelweight='light')
	# plt.rc('xtick', labelsize=10)
	
	
	# plt.rc('font', family='Times New Roman') 
	# plt.rc('font', **{'family':'serif','serif':['Times New Roman'],'size':16, 'weight':20, 'style':'normal'})
	# Set figure size and dpi
	plt.figure(figsize=(5.0, 5.0), dpi=96, facecolor='w', edgecolor='k')

	# plt.subplot(121)
	plt.plot(X1, Z1, 'k')
	plt.plot(X2, Z2, 'g')
	plt.plot(X3, Z3, 'r')
	# plt.gca().axis([0,130,0,2])
	# plt.xlabel(r'$\it{z}$, mm')
	# plt.ylabel(r'$\it{p/p}_0$')
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


def draw_X_compare(field1, field2, field3):
	X1 = field1.x * 1.0e03
	Z1 = numpy.absolute(field1.p[:, 0, 0])
	X2 = field2.x * 1.0e03
	Z2 = numpy.absolute(field2.p[:, 0, 0])
	X3 = field3.x * 1.0e03
	Z3 = numpy.absolute(field3.p[:, 0, 0])
	#-----------------------------------

	# Set Font
	# plt.rc('font', **{'family':'Times New Roman','size':16, 'weight':'bold', 'style':'normal'})
	# plt.rc('font', family='serif')
	plt.rc('font', serif='Times New Roman')
	plt.rc('font', family='serif')
	plt.rc('font', size='12')
	
	# Set figure size and dpi
	plt.figure(figsize=(7.0, 7.0), dpi=96, facecolor='w', edgecolor='w')

	# plt.subplot(121)
	plt.plot(X1, Z1, 'k', linewidth=1)
	plt.plot(X2, Z2, 'r', linewidth=2)
	plt.plot(X3, Z3, 'g', linewidth=1)
	# plt.gca().yaxis.set_ticklabels([])
	# plt.gca().axis([0,150,0,3])
	# plt.locator_params(axis = 'x', nbins = 5)
	# plt.locator_params(axis = 'y', nbins = 5)
	# plt.gca().set_aspect(50.0)
	# plt.grid()
	plt.xlabel(r'$\it{x}$, mm')
	plt.ylabel(r'$\it{p/p}_0$')
	plt.legend(['Через 1 выключ', '0.3 mm', '0.2 mm'])

	# plt.savefig(r"d:\Downloads\Calc\fig3.tiff", format='tiff', dpi=300)
	plt.savefig(r"d:\Downloads\Calc\fig_ref_focal.png", dpi=300)
	print(max(Z1))
	print(max(Z2))
	print(abs(100-100*max(Z2)/max(Z1)))
	
	plt.show()


def draw_Z(field):
	X = field.z * 1.0e03
	Z = numpy.absolute(field.p[0, 0, :])
	#-----------------------------------

	# Set Font
	# plt.rc('font', **{'family':'Times New Roman','size':16, 'weight':'bold', 'style':'normal'})
	# plt.rc('font', family='serif')
	plt.rc('font', serif='Times New Roman')
	plt.rc('font', family='serif')
	plt.rc('font', size='18')
	# plt.rc('text', usetex=True)
	# plt.rc('font', serif='Times')	
	# plt.rc('font', weight='bold')
	# plt.rc('text', usetex=True)
	# plt.rc('axes', labelweight='light')
	# plt.rc('xtick', labelsize=10)
	
	
	# plt.rc('font', family='Times New Roman') 
	# plt.rc('font', **{'family':'serif','serif':['Times New Roman'],'size':16, 'weight':20, 'style':'normal'})
	# Set figure size and dpi
	plt.figure(figsize=(10.0, 7.0), dpi=100, facecolor='w', edgecolor='k')

	# plt.subplot(121)
	plt.plot(X, Z, 'k')
	# plt.gca().axis([0,130,0,2])
	# plt.xlabel(r'$\it{z}$, mm')
	# plt.ylabel(r'$\it{p/p}_0$')
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

def draw_Z_compare(field1, field2):
	X1 = field1.z * 1.0e03
	Z1 = numpy.absolute(field1.p[0, 0, :])
	X2 = field2.z * 1.0e03
	Z2 = numpy.absolute(field2.p[0, 0, :])
	#-----------------------------------

	# Set Font
	# plt.rc('font', **{'family':'Times New Roman','size':16, 'weight':'bold', 'style':'normal'})
	# plt.rc('font', family='serif')
	plt.rc('font', serif='Times New Roman')
	plt.rc('font', family='serif')
	plt.rc('font', size='12')
	
	# Set figure size and dpi
	plt.figure(figsize=(7.0, 7.0), dpi=96, facecolor='w', edgecolor='w')

	# plt.subplot(121)
	plt.plot(X1, Z1, 'k', linewidth=1)
	plt.plot(X2, Z2, '--r', linewidth=2)
	# plt.gca().yaxis.set_ticklabels([])
	# plt.gca().axis([0,150,0,3])
	plt.locator_params(axis = 'x', nbins = 5)
	plt.locator_params(axis = 'y', nbins = 5)
	# plt.gca().set_aspect(50.0)
	# plt.grid()

	# plt.savefig(r"d:\Downloads\Calc\fig3.tiff", format='tiff', dpi=300)
	plt.savefig(r"d:\Downloads\Calc\fig_ref_axis.png", dpi=300)
	print(max(Z1))
	print(max(Z2))
	print(abs(100-100*max(Z2)/max(Z1)))
	
	plt.show()

def draw_Z_compare3(field1, field2, field3):
	X1 = field1.z * 1.0e03
	Z1 = numpy.absolute(field1.p[0, 0, :])
	X2 = field2.z * 1.0e03
	Z2 = numpy.absolute(field2.p[0, 0, :])
	X3 = field3.z * 1.0e03
	Z3 = numpy.absolute(field3.p[0, 0, :])
	#-----------------------------------

	# Set Font
	# plt.rc('font', **{'family':'Times New Roman','size':16, 'weight':'bold', 'style':'normal'})
	# plt.rc('font', family='serif')
	plt.rc('font', serif='Times New Roman')
	plt.rc('font', family='serif')
	plt.rc('font', size='24')
	
	# Set figure size and dpi
	plt.figure(figsize=(7.0, 7.0), dpi=100, facecolor='w', edgecolor='w')

	# plt.subplot(121)
	plt.plot(X1, Z1, color = '0.4', linewidth=7)
	plt.plot(X2, Z2, 'k', linewidth=2)
	plt.plot(X3, Z3, '--k', linewidth=3)

	
	plt.show()

def draw_Z_compare3_subplot(field1, field2, field3, field4):
	X1 = field1.z * 1.0e03
	Z1 = numpy.absolute(field1.p[0, 0, :])
	X2 = field2.z * 1.0e03
	Z2 = numpy.absolute(field2.p[0, 0, :])
	X3 = field3.z * 1.0e03
	Z3 = numpy.absolute(field3.p[0, 0, :])
	X4 = field4.z * 1.0e03
	Z4 = numpy.absolute(field4.p[0, 0, :])
	#-----------------------------------

	# Set Font
	# plt.rc('font', **{'family':'Times New Roman','size':16, 'weight':'bold', 'style':'normal'})
	# plt.rc('font', family='serif')
	plt.rc('font', serif='Times New Roman')
	plt.rc('font', family='serif')
	plt.rc('font', size='12')
	# plt.rc('text', usetex=True)
	
	# Set figure size and dpi
	plt.figure(dpi=96, facecolor='w', edgecolor='w')
	# plt.figure(figsize=(10.0, 5.0), dpi=96, facecolor='w', edgecolor='w')

	plt.subplot(121)
	plt.plot(X2, Z2, color = '0.4', linewidth=4)
	plt.plot(X1, Z1*1.12, '-.k', linewidth=2)
	# plt.xlabel(r'$\it{z}$, мм')
	# plt.ylabel(r'$\it{p/p}_0$')
	# plt.grid()
	plt.gca().set_aspect(0.9)
	# plt.title("(а)")
	

	plt.subplot(122)
	plt.plot(X2, Z2, color = '0.4', linewidth=4)
	plt.plot(X3, Z3, 'k', linewidth=1)
	plt.plot(X4, Z4, '--k', linewidth=2)
	plt.gca().yaxis.set_ticklabels([])
	# plt.xlabel(r'$\it{z}$, мм')
	# plt.title("(б)")
	# plt.grid()
	plt.gca().set_aspect(0.9)

	plt.savefig(r"D:\fig1.tiff", format='tiff', dpi=300)
	plt.show()


def draw_Z_compare_fig_ready(field1, field2):
	X1 = field1.z * 1.0e03
	Z1 = numpy.absolute(field1.p[0, 0, :])
	X2 = field2.z * 1.0e03
	Z2 = numpy.absolute(field2.p[0, 0, :])
	#-----------------------------------

	# Set Font
	# plt.rc('font', **{'family':'Times New Roman','size':16, 'weight':'bold', 'style':'normal'})
	# plt.rc('font', family='serif')
	plt.rc('font', serif='Times New Roman')
	plt.rc('font', family='serif')
	plt.rc('font', size='20')
	
	# Set figure size and dpi
	plt.figure(figsize=(6.0, 6.0), dpi=100, facecolor='w', edgecolor='w')

	# plt.subplot(121)
	plt.plot(X1, Z1, color = '0', linewidth=2)
	plt.plot(X2, Z2, '--k', linewidth=2)
	# plt.gca().yaxis.set_ticklabels([])
	plt.gca().axis([0,150,0,3])
	plt.locator_params(axis = 'x', nbins = 5)
	plt.locator_params(axis = 'y', nbins = 5)
	plt.grid()

	plt.savefig(r"D:\fig1.tiff", format='tiff')
	
	plt.show()


if __name__ == '__main__':
	# main_1_el()
	# main_compare()
	main()
