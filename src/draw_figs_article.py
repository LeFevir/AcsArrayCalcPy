# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import output
import numpy


def main():
	field1 = output.restore_field_from_disk(r"d:\yandex_disk\!Расчеты НОВЫЕ\Решетки\Random 128 My\one_el\field_plane")
	field2 = output.restore_field_from_disk(r"d:\yandex_disk\!Расчеты НОВЫЕ\Решетки\Random 256 Gavrilov\one_el\field_plane")
	field3 = output.restore_field_from_disk(r"d:\yandex_disk\!Расчеты НОВЫЕ\Решетки\Random 1024 Gavrilov\one_el\field_plane")

	ribs_phantom = {
		# Ребра - плоские горизонтальные полоски одинаковой толщины

		# Расстояние до плоскости ребер от решетки
		'dist_to_ribs_plane': 45.0e-03,
		# 'dist_to_ribs_plane': 65.0e-03,

		# Количество ребер
		'ribs_count': 7,

		# Ширина ребра
		'rib_width': 18.0e-03,

		# Ширина щели между ребрами
		'gap_width': 14.0e-03,

		# Координата Y нижней грани нижнего ребра
		'bottom_coord': 14.0e-03
	}

	draw_figs([field1, field2, field3], ribs_phantom)
	

def draw_figs(fields, ribs_phantom):
	draw_coefs = 1.0e03
	Y0, X0 = numpy.meshgrid(fields[0].z, fields[0].y)
	X0 *= draw_coefs
	Y0 *= draw_coefs
	X = [X0 for i in range(3)]
	Y = [Y0 for i in range(3)]
	is1 = []
	is2 = []
	y1_s = []
	y2_s = []
	
	Z = [numpy.absolute(field.p[0, :, :]) for field in fields]
	for i in range(3):
		(is1, is2) = find_isoline(Z[i], fields[i].y)
		# is1*=draw_coefs
		# is2*=draw_coefs
		y1_s.append(draw_coefs*numpy.take(fields[i].y, is1))
		y2_s.append(draw_coefs*numpy.take(fields[i].y, is2))
		# y1_s.append(draw_coefs*is1)
		# y2_s.append(draw_coefs*is2)

	# print(X)

	#-----------------------------------

	# Set Font
	plt.rc('font', **{'family':'serif','serif':['Times New Roman'],'size':9})
	# Set figure size and dpi
	fig = plt.figure(figsize=(7.0, 2.0), dpi=96, facecolor='w', edgecolor='k')
	# fig.tight_layout(pad=0.0)
	# fig = plt.figure(figsize=(14.0, 4.0), dpi=96)

	grid = AxesGrid(fig, 111, nrows_ncols=(1, 3),
		axes_pad = 0.15,
		share_all=True,
		label_mode = "L",
		cbar_location = "right",
		cbar_mode="single",
		cbar_size=0.10
		)

	# Normalizing all distributions
	max_z = numpy.max(Z)

	# Ribs
	
	gap_bot_y = 0.0 - ribs_phantom['gap_width']/2.0
	gap_top_y = ribs_phantom['gap_width']/2.0
	ribs_x = ribs_phantom["dist_to_ribs_plane"]
	ribs_x *= draw_coefs
	gap_bot_y *= draw_coefs
	gap_top_y *= draw_coefs

	# pl = plt.plot(fields[2].z*draw_coefs, y2_s[2], '--w', lw=1.5)
	ax = plt.gca()

	pc = []
	for i in range(3):
		cs = grid[i].contourf(Y[i], X[i], Z[i], 100, cmap=plt.cm.jet, vmax=max_z)
		pc.append(cs)
		grid[i].plot(fields[i].z*draw_coefs, y1_s[i], '--w', lw=1.5)
		grid[i].plot(fields[i].z*draw_coefs, y2_s[i], '--w', lw=1.5)
		grid[i].plot(ribs_x*numpy.ones(100), numpy.linspace(fields[i].y[0]*draw_coefs, gap_bot_y, 100), 'w', lw=3)
		grid[i].plot(ribs_x*numpy.ones(100), numpy.linspace(gap_top_y, fields[i].y[-1]*draw_coefs, 100), 'w', lw=3)
		# grid[i].annotate('ребро', xy=(40, -15), xytext=(40, -15), color ='w')

	grid.cbar_axes[0].colorbar(pc[0])
	# This affects all axes as share_all = True.
	grid.axes_llc.set_aspect('equal')
	grid.axes_llc.locator_params(nbins=5)
	grid.axes_llc.set_xlim([0, 60])

	# grid.axes_llc.set_clim([0, 2.5])

	# plt.xlabel(r'$\it{z}$, mm')
	# plt.ylabel(r'$\it{y}$, mm')
	# plt.ylim(0, 50)
	# plt.title('Amplitude')
	# plt.gca().set_aspect('equal')
	
	
	# plane_y = numpy.arange(-20,20,1)
	# plane_z = 45.0
	# plt.plot(plane_z*numpy.ones(40), plane_y, '--w')

	# plt.savefig(r"d:\Downloads\Calc\fig_table_iso.png", dpi=300, pad_inches=0.0)
	# plt.savefig(r"d:\Downloads\Calc\fig3.eps", dpi=600, pad_inches=0.0)
	plt.savefig(r"d:\Downloads\Calc\fig3.png", dpi=300, pad_inches=0.0)
	plt.show()


def find_isoline(z, y):
	# На каждой линии по Х найти линию по -6дБ
	(x_len, y_len) = numpy.shape(z)
	isoline_1 = []
	isoline_2 = []
	# print(z.max())
	# print(numpy.where(z==z.max()))
	for i in range(y_len):
		ar = z[:, i]
		m = numpy.max(ar)
		l = m*0.5
		m_x = numpy.argmax(ar)
		if len(numpy.abs(ar[:m_x]-l)) != 0:
			l1 = numpy.abs(ar[:m_x]-l).argmin()
			isoline_1.append(l1)
		else:
			isoline_1.append(0)
		if len(numpy.abs(ar[m_x:]-l)) != 0:
			l2 = numpy.abs(ar[m_x:]-l).argmin()
			isoline_2.append(l2+m_x)
		else:
			isoline_2.append(x_len)

	return isoline_1, isoline_2

if __name__ == '__main__':
	main()
	# draw_figs()
