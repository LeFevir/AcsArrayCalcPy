# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy
import transducer
import rayleigh
import codecs
import output


def main():
	coef_fill = 0.52
	tr = []
	trans1 = transducer.transducer_from_file(r"d:\yandex_disk\!Расчеты НОВЫЕ\Решетки\Ideal\Ribs off\array.txt")	
	tris = rayleigh.restore_tris_from_disk(r'd:\yandex_disk\!Расчеты НОВЫЕ\Решетки\Ideal\Ribs off\tris.bin')
	tr.append({'trans': trans1, 'tris': tris})
	trans2 = transducer.transducer_from_file(r"d:\yandex_disk\!Расчеты НОВЫЕ\Решетки\Random 128 My\array.txt")
	trans2.add_elements_from_file(r"d:\yandex_disk\!Расчеты НОВЫЕ\Решетки\Random 128 My\array_elements.txt")
	trans2_off_file = r"d:\yandex_disk\!Расчеты НОВЫЕ\Решетки\Random 128 My\off\switched_off.txt"
	tr.append({'trans': trans2, 'off': trans2_off_file})
	trans3 = transducer.transducer_from_file(r"d:\yandex_disk\!Расчеты НОВЫЕ\Решетки\Random 256 Gavrilov\array.txt")
	trans3.add_elements_from_file(r"d:\yandex_disk\!Расчеты НОВЫЕ\Решетки\Random 256 Gavrilov\array_elements.txt")
	trans3_off_file = r"d:\yandex_disk\!Расчеты НОВЫЕ\Решетки\Random 256 Gavrilov\off_el\switched_off.txt"
	tr.append({'trans': trans3, 'off': trans3_off_file})
	trans4 = transducer.transducer_from_file(r"d:\yandex_disk\!Расчеты НОВЫЕ\Решетки\Random 1024 Gavrilov\array.txt")
	trans4.add_elements_from_file(r"d:\yandex_disk\!Расчеты НОВЫЕ\Решетки\Random 1024 Gavrilov\array_elements.txt")
	trans4_off_file = r"d:\yandex_disk\!Расчеты НОВЫЕ\Решетки\Random 1024 Gavrilov\off_el\switched_off.txt"
	tr.append({'trans': trans4, 'off': trans4_off_file})

	field = output.restore_field_from_disk(r"d:\yandex_disk\!Расчеты НОВЫЕ\Решетки\Ideal\Ribs off\field_plane_ribs")
	max_p = numpy.max(numpy.absolute(field.p))
	print(max_p**2)
	print(rayleigh.calc_power(field))
	field.p *= numpy.sqrt(coef_fill)
	field.vn *= numpy.sqrt(coef_fill)
	tr[0]['ribs'] = field
	field = output.restore_field_from_disk(r"d:\yandex_disk\!Расчеты НОВЫЕ\Решетки\Random 128 My\off_new\field_plane_ribs")
	max_p = numpy.max(numpy.absolute(field.p))
	print(max_p**2)
	print(rayleigh.calc_power(field))
	tr[1]['ribs'] = field
	field = output.restore_field_from_disk(r"d:\yandex_disk\!Расчеты НОВЫЕ\Решетки\Random 256 Gavrilov\off_el_new\field_plane_ribs")
	max_p = numpy.max(numpy.absolute(field.p))
	print(max_p**2)
	print(rayleigh.calc_power(field))
	tr[2]['ribs'] = field
	field = output.restore_field_from_disk(r"d:\yandex_disk\!Расчеты НОВЫЕ\Решетки\Random 1024 Gavrilov\off_el_new\field_plane_ribs")
	max_p = numpy.max(numpy.absolute(field.p))
	print(max_p**2)
	print(rayleigh.calc_power(field))
	tr[3]['ribs'] = field




	field = output.restore_field_from_disk(r"d:\yandex_disk\!Расчеты НОВЫЕ\Решетки\Ideal\Ribs off\focus_plane\field_plane")
	max_p = numpy.max(numpy.absolute(field.p))
	print(max_p**2)
	print(rayleigh.calc_power(field))
	field.p *= numpy.sqrt(coef_fill)
	field.vn *= numpy.sqrt(coef_fill)
	tr[0]['focal'] = field
	field = output.restore_field_from_disk(r"d:\yandex_disk\!Расчеты НОВЫЕ\Решетки\Random 128 My\off_new\focus_plane\field_plane")
	max_p = numpy.max(numpy.absolute(field.p))
	print(max_p**2)
	print(rayleigh.calc_power(field))
	tr[1]['focal'] = field
	field = output.restore_field_from_disk(r"d:\yandex_disk\!Расчеты НОВЫЕ\Решетки\Random 256 Gavrilov\off_el_new\focus_plane\field_plane")
	max_p = numpy.max(numpy.absolute(field.p))
	print(max_p**2)
	print(rayleigh.calc_power(field))
	tr[2]['focal'] = field
	field = output.restore_field_from_disk(r"d:\yandex_disk\!Расчеты НОВЫЕ\Решетки\Random 1024 Gavrilov\off_el_new\focus_plane\field_plane")
	max_p = numpy.max(numpy.absolute(field.p))
	print(max_p**2)
	print(rayleigh.calc_power(field))
	tr[3]['focal'] = field
	
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

	# draw_figs(trans1, tris, trans2, trans2_off_file, trans3, trans3_off_file, trans4, trans4_off_file)
	draw_figs(tr)


def draw_transducer_tris(trans, tris, ax):
	# Draw transducer trans
	#Draw in mm's
	draw_coef = 1000

	#On elements color    
	on_color = '0.6'
	#Off elements color
	off_color = 'w'

	# Set Font
	# plt.rc('font',**{'family':'serif','serif':['Times New Roman'],'size':12})
	# Set figure size and dpi
	# fig = plt.figure(figsize=(7.0, 7.0), dpi=96, facecolor='w')
		
	# # Draw central hole
	# circle = plt.Circle((0.0, 0.0), radius=trans.hole_radius*draw_coef, linewidth=1, facecolor='w', edgecolor='k')
	# plt.gca().add_patch(circle)

	x_s = [tri.xc*draw_coef for tri in tris]
	y_s = [tri.yc*draw_coef for tri in tris]
	z_s = [tri.zc*draw_coef for tri in tris]

	# Draw elements
	ax.plot(x_s, y_s,'.', color=on_color, lw=0)
	# plt.xlabel(r'$\it{x}$, mm')
	# plt.ylabel(r'$\it{y}$, mm')
	# plt.gca().set_aspect('equal')
	ax.axis([-100, 100, -100,100])
	# Draw aperture of trans as a circle
	circle = plt.Circle((0.0, 0.0), radius=trans.aperture*draw_coef/2, linewidth=2, facecolor='w', edgecolor='k')
	ax.add_patch(circle)


def draw_transducer(trans, off_file, ax):
	# Draw transducer trans
	#Draw in mm's
	draw_coef = 1000

	#On elements color    
	on_color = '0.6'
	#Off elements color
	off_color = 'w'

	# Draw aperture of trans as a circle
	circle = plt.Circle((0.0, 0.0), radius=trans.aperture*draw_coef/2, linewidth=2, facecolor='w', edgecolor='k')
	ax.add_patch(circle)
		
	# Draw central hole
	circle = plt.Circle((0.0, 0.0), radius=trans.hole_radius*draw_coef, linewidth=1, facecolor='w', edgecolor='k')
	ax.add_patch(circle)

	off_ids = []
	if off_file != None:
		with codecs.open(off_file, encoding='utf-8-sig') as f:
			off_ids = [int(line.strip()) for line in f]
	
	# Draw elements
	for el in trans.elements:
		if el['id'] in off_ids:
			circle = plt.Circle((el['center_x']*draw_coef, el['center_y']*draw_coef), radius=trans.element_radius*draw_coef, linewidth=trans.element_radius*draw_coef*0.10, facecolor=off_color, edgecolor='k')
		else:
			circle = plt.Circle((el['center_x']*draw_coef, el['center_y']*draw_coef), radius=trans.element_radius*draw_coef, linewidth=trans.element_radius*draw_coef*0.10, facecolor=on_color, edgecolor='k')
		ax.add_patch(circle)

	# ax.axis('scaled') # to scale circles properly

	# plt.xlabel(r'$\it{x}$, mm')
	# plt.ylabel(r'$\it{y}$, mm')


def draw_figs(tr):
	#-----------------------------------

	# Set Font
	plt.rc('font', **{'family':'serif','serif':['Times New Roman'],'size':9})
	# Set figure size and dpi
	fig = plt.figure(figsize=(7.0, 5.0), dpi=96, facecolor='w', edgecolor='k')
	# fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
	# fig.tight_layout(pad=0.0)
	# fig = plt.figure(figsize=(14.0, 4.0), dpi=96)

	grid = AxesGrid(fig, (0.1, 0.25, 0.8, 0.9), nrows_ncols=(2, 4),
		axes_pad = 0.35,
		share_all=True,
		label_mode = "L",
		cbar_location = "right",
		cbar_mode="edge",
		cbar_size="10%",
		cbar_pad="5%"
		)

	grid2 = AxesGrid(fig, (0.1, 0.1, 0.8, 0.25), nrows_ncols=(1, 4),
		axes_pad = 0.35,
		share_all=True,
		label_mode = "L",
		cbar_location = "right",
		cbar_mode="edge",
		cbar_size="10%",
		cbar_pad="5%"
		)
	# grid = AxesGrid(fig, 212, nrows_ncols=(1, 4),
	# 	axes_pad = 0.35,
	# 	share_all=True,
	# 	label_mode = "L",
	# 	cbar_location = "right",
	# 	cbar_mode="None",
	# 	cbar_size=0.10
	# 	)
	cs1 = []
	cs2 = []
	draw_transducer_tris(tr[0]['trans'], tr[0]['tris'], grid[0])
	for i in range(1,4):
		draw_transducer(tr[i]['trans'], tr[i]['off'], grid[i])
	for i in range(0,4):
		cs1.append(draw_XY_with_ribs(tr[i]['ribs'], grid[i+4]))
	for i in range(0,4):
		cs2.append(draw_XY(tr[i]['focal'], grid2[i]))

	# grid[9].axis([-25, 25, -25, 25])
	# draw_transducer(trans2, trans2_off_file, grid[1])
	# draw_transducer(trans3, trans3_off_file, grid[2])
	# draw_transducer(trans4, trans4_off_file, grid[3])
	# grid[i].plot(fields[i].z*1000, y2_s[i], '--w', lw=1.5)
	# grid[i].plot(ribs_x*numpy.ones(100), numpy.linspace(rib1_bot_y, rib1_top_y, 100), 'w', lw=3)
	# grid[i].annotate('ребро', xy=(40, -15), xytext=(40, -15), color ='w')

	fig.delaxes(grid.cbar_axes[0])
	grid.cbar_axes[1].colorbar(cs1[1])
	grid2.cbar_axes[0].colorbar(cs2[1])
	# This affects all axes as share_all = True.
	grid.axes_llc.set_aspect('equal')
	# grid.axes_llc.set_xlim([0, ])
	# grid2.axes_llc.set_clim([0, 1300])
	# grid.axes_llc.locator_params(nbins=5)
	# grid.axes_llc.set_xlim([0, 60])

	# grid.axes_llc.set_clim([0, 2.5])

	# plt.xlabel(r'$\it{z}$, mm')
	# plt.ylabel(r'$\it{y}$, mm')
	# plt.ylim(0, 50)
	# plt.title('Amplitude')
	# plt.gca().set_aspect('equal')
	
	
	# plane_y = numpy.arange(-20,20,1)
	# plane_z = 45.0
	# plt.plot(plane_z*numpy.ones(40), plane_y, '--w')

	plt.savefig(r"d:\Downloads\Calc\fig_table_surf.png", dpi=300, pad_inches=0.0)
	# plt.savefig(r"d:\Downloads\Calc\fig2.eps", dpi=600, pad_inches=0.0)
	plt.savefig(r"d:\Downloads\Calc\fig2.png", dpi=300, pad_inches=0.0)
	# plt.savefig(r"d:\Downloads\Calc\fig2.", dpi=600, pad_inches=0.0)
	plt.show()


def find_isoline(z, y):
	# На каждой линии по Х найти линию по -3дБ
	(x_len, y_len) = numpy.shape(z)
	isoline_1 = []
	isoline_2 = []
	# y1_s = []
	# y2_s = []
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

	# y1_s.append(numpy.take(y, isoline_1))
	# y2_s.append(numpy.take(y, isoline_2))

	return isoline_1, isoline_2
	# return y1_s, y2_s


def draw_XY_with_ribs(field, ax):
	Y, X = numpy.meshgrid(field.y, field.x)
	X *= 1.0e03
	Y *= 1.0e03
	Z = numpy.absolute(field.p[:, :, 0]) * numpy.absolute(field.p[:, :, 0])
	# Z =  0.5 * numpy.real(field.p[:,:,0] * numpy.conj(field.vn[:,:,0]*1000.0 * 1500.0))

	# Z *= 1000.0 * 1500.0
	# print("Full Energy = ", numpy.sum(Z)*0.1*0.1*1.0e-06)
	# print(numpy.max(Z))

	#-----------------------------------
	cs = ax.contourf(X, Y, Z, 100, cmap=plt.cm.jet, vmax=5.65)
	# plt.colorbar(orientation = 'vertical')
	# plt.xlabel('x, mm')
	# plt.xlabel(r'$\it{x}$, mm')
	# plt.ylabel(r'$\it{y}$, mm')
	# plt.ylabel('y, mm')
	# plt.locator_params(nbins=6)
	# cb = plt.colorbar(cs, orientation='vertical', format='%.1f')
	# cb.set_label(r'$\it{p/p_{0}}$')
	# cb.locator = plt.MaxNLocator(nbins=6)
	# cb.update_ticks()
	# plt.gca().set_aspect('equal')

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
	ax.plot(rib_x, rib1_bot_y*numpy.ones(200), '-.w', lw=0.6)
	ax.plot(rib_x, rib2_bot_y*numpy.ones(200), '-.w', lw=0.6)
	ax.plot(rib_x, rib3_top_y*numpy.ones(200), '-.w', lw=0.6)
	ax.plot(rib_x, rib2_top_y*numpy.ones(200), '-.w', lw=0.6)
	ax.plot(rib_x, rib1_top_y*numpy.ones(200), '-.w', lw=0.6)
	ax.plot(rib_x, rib3_bot_y*numpy.ones(200), '-.w', lw=0.6)
	ax.plot(rib_x, rib4_top_y*numpy.ones(200), '-.w', lw=0.6)
	ax.plot(rib_x, rib4_bot_y*numpy.ones(200), '-.w', lw=0.6)
	ax.plot(rib_x, rib5_top_y*numpy.ones(200), '-.w', lw=0.6)
	ax.plot(rib_x, rib5_bot_y*numpy.ones(200), '-.w', lw=0.6)

	return cs


def draw_XY(field, ax):
	Y, X = numpy.meshgrid(field.y, field.x)
	X *= 1.0e03
	Y *= 1.0e03
	Z = numpy.absolute(field.p[:, :, 0]) * numpy.absolute(field.p[:, :, 0])
	# Z2 = numpy.angle(field.p[:, :, 0])

	#-----------------------------------

	# print numpy.max(Z)
	# print numpy.shape(X)
	
	cs = ax.contourf(X, Y, Z, 100, cmap=plt.cm.jet, vmax=500)
	# plt.colorbar(orientation = 'vertical')
	# plt.xlabel('x, mm')
	# plt.xlabel(r'$\it{x}$, mm')
	# plt.ylabel(r'$\it{y}$, mm')
	# plt.ylabel('y, mm')
	# plt.locator_params(nbins=6)
	# cb = plt.colorbar(cs, orientation='vertical', format='%.1f')
	# cb.set_label(r'$\it{p/p_{0}}$')
	# cb.locator = plt.MaxNLocator(nbins=6)
	# cb.update_ticks()
	# plt.gca().set_aspect('equal')

	return cs

if __name__ == '__main__':
	main()
	# draw_figs()
