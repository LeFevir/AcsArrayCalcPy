# -*- coding: utf-8 -*-
import rayleigh


def test_slice_trans_plane_on_tris():
	tris = rayleigh.slice_trans_plane_on_tris(aperture=2.0, num_of_tris_on_line=4)
	assert len(tris) == 8


def test_compute_centroids():
	tri = rayleigh.Triangle(0.0, 0.0, 3.0, 3.0, 3.0, 0.0)
	rayleigh.project_tris_on_trans((tri,), 10.0)
	tri.compute_centroids()
	assert tri.xc == 2.0
	assert tri.yc == 1.0
	assert tri.zc != 0.0


def test_project_tris_on_trans():
	tri = rayleigh.Triangle(1.0, 2.0, 3.0, 3.0, 3.0, 3.0)
	tri.project_on_trans(10.0)
	assert tri.z1 != 0.0
	assert tri.z2 != 0.0
	assert tri.z3 != 0.0

if __name__ == '__main__':
	pass
