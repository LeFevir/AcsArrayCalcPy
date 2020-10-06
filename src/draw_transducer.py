# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import transducer
import codecs

import rayleigh

def draw_transducer_3_compare(trans1, trans2, trans3):
	off=None
	# Draw transducer trans
	#Draw in mm's
	draw_coef = 1000

	#On elements color    
	on_color = '0.8'
	#Off elements color
	off_color = 'w'

	# Set Font
	plt.rc('font',**{'family':'serif','serif':['Times New Roman'],'size':12})
	# Set figure size and dpi
	plt.figure(figsize=(6.5, 5.0), dpi=96, facecolor='w')

	plt.subplot(131)
	# Draw aperture of trans as a circle
	circle = plt.Circle((0.0, 0.0), radius=trans1.aperture*draw_coef/2, linewidth=1, facecolor='w', edgecolor='k')
	plt.gca().add_patch(circle)
		
	# Draw central hole
	circle = plt.Circle((0.0, 0.0), radius=trans1.hole_radius*draw_coef, linewidth=1, facecolor='w', edgecolor='k')
	plt.gca().add_patch(circle)

	off_ids = []
	if off != None:
		with codecs.open(off, encoding='utf-8-sig') as f:
			off_ids = [int(line.strip()) for line in f]
	
	# Draw elements
	for el in trans1.elements:
		if el['id'] in off_ids:
			circle = plt.Circle((el['center_x']*draw_coef, el['center_y']*draw_coef), radius=trans1.element_radius*draw_coef, linewidth=0.1, facecolor=off_color, edgecolor='k')
		else:
			circle = plt.Circle((el['center_x']*draw_coef, el['center_y']*draw_coef), radius=trans1.element_radius*draw_coef, linewidth=0.7, facecolor=on_color, edgecolor='k')
		plt.gca().add_patch(circle)

	plt.axis('scaled') # to scale circles properly

	# plt.xlabel(r'$\it{x}$, mm')
	# plt.ylabel(r'$\it{y}$, mm')

	plt.subplot(132)
	# Draw aperture of trans as a circle
	circle = plt.Circle((0.0, 0.0), radius=trans2.aperture*draw_coef/2, linewidth=1, facecolor='w', edgecolor='k')
	plt.gca().add_patch(circle)
		
	# Draw central hole
	circle = plt.Circle((0.0, 0.0), radius=trans2.hole_radius*draw_coef, linewidth=1, facecolor='w', edgecolor='k')
	plt.gca().add_patch(circle)

	off_ids = []
	if off != None:
		with codecs.open(off, encoding='utf-8-sig') as f:
			off_ids = [int(line.strip()) for line in f]
	
	# Draw elements
	for el in trans2.elements:
		if el['id'] in off_ids:
			circle = plt.Circle((el['center_x']*draw_coef, el['center_y']*draw_coef), radius=trans2.element_radius*draw_coef, linewidth=0.1, facecolor=off_color, edgecolor='k')
		else:
			circle = plt.Circle((el['center_x']*draw_coef, el['center_y']*draw_coef), radius=trans2.element_radius*draw_coef, linewidth=0.7, facecolor=on_color, edgecolor='k')
		plt.gca().add_patch(circle)

	plt.axis('scaled') # to scale circles properly
	# plt.xlabel(r'$\it{x}$, mm')
	# plt.ylabel(r'$\it{y}$, mm')
	plt.gca().yaxis.set_ticklabels([])

	plt.subplot(133)
	# Draw aperture of trans as a circle
	circle = plt.Circle((0.0, 0.0), radius=trans3.aperture*draw_coef/2, linewidth=1, facecolor='w', edgecolor='k')
	plt.gca().add_patch(circle)
		
	# Draw central hole
	circle = plt.Circle((0.0, 0.0), radius=trans3.hole_radius*draw_coef, linewidth=1, facecolor='w', edgecolor='k')
	plt.gca().add_patch(circle)

	off_ids = []
	if off != None:
		with codecs.open(off, encoding='utf-8-sig') as f:
			off_ids = [int(line.strip()) for line in f]
	
	# Draw elements
	for el in trans3.elements:
		if el['id'] in off_ids:
			circle = plt.Circle((el['center_x']*draw_coef, el['center_y']*draw_coef), radius=trans3.element_radius*draw_coef, linewidth=0.1, facecolor=off_color, edgecolor='k')
		else:
			circle = plt.Circle((el['center_x']*draw_coef, el['center_y']*draw_coef), radius=trans3.element_radius*draw_coef, linewidth=0.7, facecolor=on_color, edgecolor='k')
		plt.gca().add_patch(circle)

	plt.axis('scaled') # to scale circles properly
	# plt.xlabel(r'$\it{x}$, mm')
	# plt.ylabel(r'$\it{y}$, mm')
	plt.tight_layout()
	plt.gca().yaxis.set_ticklabels([])
	
	plt.savefig(r"C:\Downloads\fig_transes.tiff", format='tiff', dpi=300)
	plt.show()


def draw_transducer(trans, off=None):
	# Draw transducer trans
	#Draw in mm's
	draw_coef = 1000

	#On elements color    
	on_color = '0.6'
	#Off elements color
	off_color = 'w'

	# Set Font
	plt.rc('font',**{'family':'serif','serif':['Times New Roman'],'size':14})
	# Set figure size and dpi
	plt.figure(figsize=(5.0, 5.0), dpi=100, facecolor='w')
	
	# Draw aperture of trans as a circle
	circle = plt.Circle((0.0, 0.0), radius=trans.aperture*draw_coef/2, linewidth=2, facecolor='w', edgecolor='k')
	plt.gca().add_patch(circle)
		
	# Draw central hole
	circle = plt.Circle((0.0, 0.0), radius=trans.hole_radius*draw_coef, linewidth=1, facecolor='w', edgecolor='k')
	plt.gca().add_patch(circle)

	off_ids = []
	if off != None:
		with codecs.open(off, encoding='utf-8-sig') as f:
			off_ids = [int(line.strip()) for line in f]
	
	# Draw elements
	for el in trans.elements:
		if el['id'] in off_ids:
			circle = plt.Circle((el['center_x']*draw_coef, el['center_y']*draw_coef), radius=trans.element_radius*draw_coef, linewidth=1, facecolor=off_color, edgecolor='k')
		else:
			circle = plt.Circle((el['center_x']*draw_coef, el['center_y']*draw_coef), radius=trans.element_radius*draw_coef, linewidth=1, facecolor=on_color, edgecolor='k')
		plt.gca().add_patch(circle)

	plt.axis('scaled') # to scale circles properly

	plt.xlabel(r'$\it{x}$, mm')
	plt.ylabel(r'$\it{y}$, mm')

	
	plt.savefig(r"fig_trans.png")
	plt.show()


def draw_transducer_tris(trans, tris):
	# Draw transducer trans
	#Draw in mm's
	draw_coef = 1000

	#On elements color    
	on_color = '0.6'
	#Off elements color
	off_color = 'w'

	# Set Font
	plt.rc('font',**{'family':'serif','serif':['Times New Roman'],'size':12})
	# Set figure size and dpi
	plt.figure(figsize=(7.0, 7.0), dpi=96, facecolor='w')
		
	# # Draw central hole
	# circle = plt.Circle((0.0, 0.0), radius=trans.hole_radius*draw_coef, linewidth=1, facecolor='w', edgecolor='k')
	# plt.gca().add_patch(circle)

	x_s = [tri.xc*draw_coef for tri in tris]
	y_s = [tri.yc*draw_coef for tri in tris]
	z_s = [tri.zc*draw_coef for tri in tris]

	# Draw elements
	plt.figure(1)
	plt.plot(x_s, y_s,'.', color='0.4', lw=0)
	plt.xlabel(r'$\it{x}$, mm')
	plt.ylabel(r'$\it{y}$, mm')
	plt.gca().set_aspect('equal')
	plt.axis([-100, 100, -100,100])
	# Draw aperture of trans as a circle
	circle = plt.Circle((0.0, 0.0), radius=trans.aperture*draw_coef/2, linewidth=2, facecolor='w', edgecolor='k')
	plt.gca().add_patch(circle)
	
	# plt.figure(2)
	# plt.plot(z_s, y_s,'r.')
	# plt.xlabel(r'$\it{z}$, mm')
	# plt.ylabel(r'$\it{y}$, mm')
	# plt.gca().set_aspect('equal')

	plt.savefig(r"d:\yandex_disk\!Расчеты НОВЫЕ\Решетки\Ideal\Ribs off\fig_trans.png")
	plt.show()


if __name__ == '__main__':
	# test case
	# trans = transducer.transducer_from_file(r"..\test_files\array.txt")
	# trans.add_elements_from_file(r"..\test_files\array_elements.txt")
	# draw_transducer_3d(trans)

	
	# trans = transducer.transducer_from_file(r"h:\Downloads\New_calc\array.txt")
	# trans.add_elements_from_file(r"h:\Downloads\New_calc\array_elements.txt")
	# draw_transducer(trans, off = r"h:\Downloads\New_calc\switched_off.txt")
	draw_transducer(trans)

	# Draw 3 trans
	# trans1 = transducer.transducer_from_file(r"c:\YandexDisk\!Расчеты НОВЫЕ\Решетки\Random 256 Gavrilov\array.txt")
	# trans1.add_elements_from_file(r"c:\YandexDisk\!Расчеты НОВЫЕ\Решетки\Random 256 Gavrilov\array_elements.txt")
	# trans2 = transducer.transducer_from_file(r"c:\YandexDisk\!Расчеты НОВЫЕ\Решетки\Regular 256 Gavrilov\array.txt")
	# trans2.add_elements_from_file(r"c:\YandexDisk\!Расчеты НОВЫЕ\Решетки\Regular 256 Gavrilov\array_elements.txt")
	# trans3 = transducer.transducer_from_file(r"c:\YandexDisk\!Расчеты НОВЫЕ\Решетки\Random 1024 Gavrilov\array.txt")
	# trans3.add_elements_from_file(r"c:\YandexDisk\!Расчеты НОВЫЕ\Решетки\Random 1024 Gavrilov\array_elements.txt")
	# draw_transducer_3_compare(trans1, trans2, trans3)
	
	# trans = transducer.transducer_from_file(r"d:\yandex_disk\!Расчеты НОВЫЕ\Решетки\Ideal\Ribs off\array.txt")	
	# tris = rayleigh.restore_tris_from_disk(r'c:\Downloads\Calc\2014-07-02_14-47-10\tris.bin')
	# tris = rayleigh.restore_tris_from_disk(r'd:\yandex_disk\!Расчеты НОВЫЕ\Решетки\Ideal\Ribs off\tris.bin')
	# draw_transducer_tris(trans, tris)
