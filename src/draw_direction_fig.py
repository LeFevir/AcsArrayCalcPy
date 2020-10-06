# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt


def main():
	draw_dir()

def calc_dir(a,c,f):
	"""Диаграмма направленности одного элемента

	"""
	import math
	import numpy
	from scipy.special import j1

	
	wave_number = 2 * math.pi * f / c

	theta_min = -math.pi/2
	theta_max = math.pi/2

	thetas = numpy.linspace(theta_min, theta_max, 100, retstep=False)
	sines = numpy.sin(thetas)

	Dir = numpy.abs(2*j1(wave_number*a*sines) / (wave_number*a*sines))

	return numpy.degrees(thetas), Dir


def draw_dir():
	a = 3.5e-03   # element_radius
	c = 1500 # speed of sound
	f = 1.0e06 # freq

	deg1, d1 = calc_dir(a,c,f)

	deg2, d2 = calc_dir(a/2,c,f)

	deg3, d3 = calc_dir(a*2,c,f)
	#-----------------------------------

	# print numpy.max(Z)
	# print numpy.shape(X)

	# Set Font
	plt.rc('font', **{'family':'serif','serif':['Times New Roman'],'size':14})
	# Set figure size and dpi
	plt.figure(figsize=(5.0, 5.0), dpi=100, facecolor='w', edgecolor='k')

	plt.plot(deg1, d1, 'k')
	plt.plot(deg2, d2, '-.k')
	plt.plot(deg3, d3, '--k')
	# plt.xlabel('x, mm')
	plt.xlabel(r'$\it{\theta}$, $^{\circ}$')
	plt.legend([r"a = 3.5 мм = 2.33$\lambda$", r"a = 1.75 мм = 1.16$\lambda$", r"a = 7.0 мм"])
	plt.ylabel(r'$\it{F(\theta)}$')
	# plt.ylabel('y, mm')
	# plt.locator_params(nbins=6)
	# cb = plt.colorbar(cs, orientation='vertical', format='%.1f')
	# cb.set_label(r'$\it{p/p_{0}}$')
	# cb.locator = plt.MaxNLocator(nbins=6)
	# cb.update_ticks()
	# plt.gca().set_aspect('equal')

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

	# plt.savefig(r"fig_XY.png")
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
	plt.figure(figsize=(6.0, 8.0), dpi=100, facecolor='w', edgecolor='k')

	# plt.subplot(121)
	plt.contourf(X, Y, Z, 100, cmap=plt.cm.jet)
	plt.colorbar(orientation = 'vertical')
	plt.xlabel(r'$\it{x}$, mm')
	plt.ylabel(r'$\it{z}$, mm')

	plt.ylim(0, 50)
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

	plt.savefig(r"fig_XZ.png")
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
	plt.figure(figsize=(10.0, 8.0), dpi=100, facecolor='w', edgecolor='k')

	# plt.subplot(121)
	plt.contourf(Y, X, Z, 100, cmap=plt.cm.jet)
	plt.colorbar(orientation = 'vertical')
	plt.xlabel(r'$\it{z}$, mm')
	plt.ylabel(r'$\it{y}$, mm')
	# plt.ylim(0, 50)
	# plt.title('Amplitude')
	plt.gca().set_aspect('equal')

	# plt.subplot(122)
	# plt.contourf(X,Y,Z2,20,cmap=plt.cm.jet)
	# plt.title('Phase')
	# plt.xlabel(r'$\it{x}$, mm')
	# plt.gca().set_aspect('equal')

	# plt.suptitle("Distance from trans = 10 mm")

	# plt.savefig(ur"d:\Educ\АКУСТИКА\Philips R&D\AcField Cylindrical shell\Py\fig_YZ.png")
	plt.show()


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
	plt.figure(figsize=(8.0, 7.0), dpi=100, facecolor='w', edgecolor='k')

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
	plt.figure(figsize=(10.0, 4.0), dpi=100, facecolor='w', edgecolor='k')

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
	plt.figure(figsize=(10.0, 4.0), dpi=100, facecolor='w', edgecolor='k')

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

if __name__ == '__main__':
	main()
