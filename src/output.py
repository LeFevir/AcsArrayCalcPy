# -*- coding: utf-8 -*-
"""

"""
import struct
import logging
import pickle


def save_field_on_disk(field, file_name="field"):
    logging.info("Saving Field on disk...")
    with open(file_name, 'wb') as f:
        pickle.dump(field, f)


def restore_field_from_disk(file_name="field"):
    logging.info("Restoring Field...")
    with open(file_name, 'rb') as f:
        field = pickle.load(f)
    return field


def PrintCylAmpPhasePhiZ_Binary(field, pointOnRo):
    logging.info("Binary FieldCylComplex output has been started")

    file_amp = open('amp_pres_phiz.bin', 'wb')
    file_phase = open('phase_pres_phiz.bin', 'wb')

    file_amp.write(struct.pack('<qq', field.n_phi, field.n_z))
    file_phase.write(struct.pack('<qq', field.n_phi, field.n_z))

    i = pointOnRo
    for j, phi in enumerate(field.phi):
        for k, z in enumerate(field.z):
            file_amp.write(struct.pack('<ddd', phi, z, field.amp(i, j, k)))
            file_phase.write(struct.pack('<ddd', phi, z, field.phase(i, j, k)))

    file_amp.close()
    file_phase.close()

    logging.info("Binary FieldCylComplex output has been finished")


def PrintCylAmpPhaseYZ_Binary(field, point_on_x):
    logging.info("Binary Field Plane output has been started")

    file_amp = open('amp_pres_yz.bin', 'wb')
    file_phase = open('phase_pres_yz.bin', 'wb')

    file_amp.write(struct.pack('<qq', field.n_y, field.n_z))
    file_phase.write(struct.pack('<qq', field.n_y, field.n_z))

    i = point_on_x
    for j, y in enumerate(field.y):
        for k, z in enumerate(field.z):
            file_amp.write(struct.pack('<ddd', y, z, field.amp(i, j, k)))
            file_phase.write(struct.pack('<ddd', y, z, field.phase(i, j, k)))

    file_amp.close()
    file_phase.close()

    logging.info("Binary Field Plane output has been finished")


def PrintCylAmpPhaseXZ_Binary(field, point_on_y):
    logging.info("Binary Field Plane output has been started")

    file_amp = open('amp_pres_xz.bin', 'wb')
    file_phase = open('phase_pres_xz.bin', 'wb')

    file_amp.write(struct.pack('<qq', field.n_x, field.n_z))
    file_phase.write(struct.pack('<qq', field.n_x, field.n_z))

    j = point_on_y
    for i, x in enumerate(field.x):
        for k, z in enumerate(field.z):
            file_amp.write(struct.pack('<ddd', x, z, field.amp(i, j, k)))
            file_phase.write(struct.pack('<ddd', x, z, field.phase(i, j, k)))

    file_amp.close()
    file_phase.close()

    logging.info("Binary Field Plane output has been finished")


def PrintCylAmpPhaseXY_Binary(field, point_on_z):
    logging.info("Binary Field Plane output has been started")

    file_amp = open('amp_pres_xy.bin', 'wb')
    file_phase = open('phase_pres_xy.bin', 'wb')

    file_amp.write(struct.pack('<qq', field.n_x, field.n_y))
    file_phase.write(struct.pack('<qq', field.n_x, field.n_y))

    k = point_on_z
    for i, x in enumerate(field.x):
        for j, y in enumerate(field.y):
            file_amp.write(struct.pack('<ddd', x, y, field.p_amp(i, j, k)))
            file_phase.write(struct.pack('<ddd', x, y, field.p_phase(i, j, k)))

    file_amp.close()
    file_phase.close()

    logging.info("Binary Field Plane output has been finished")


def PrintCylAmpPhaseX_Binary(field, point_on_y, point_on_z):
    logging.info("Binary Field Plane output has been started")

    file_amp = open('amp_pres_x.bin', 'wb')
    file_phase = open('phase_pres_x.bin', 'wb')

    file_amp.write(struct.pack('<q', field.n_x))
    file_phase.write(struct.pack('<q', field.n_x))

    k = point_on_z
    j = point_on_y
    for i, x in enumerate(field.x):
            file_amp.write(struct.pack('<dd', x, field.amp(i, j, k)))
            file_phase.write(struct.pack('<dd', x, field.phase(i, j, k)))

    file_amp.close()
    file_phase.close()

    logging.info("Binary Field Plane output has been finished")


def PrintCylAmpPhaseY_Binary(field, point_on_x, point_on_z):
    logging.info("Binary Field Plane output has been started")

    file_amp = open('amp_pres_y.bin', 'wb')
    file_phase = open('phase_pres_y.bin', 'wb')

    file_amp.write(struct.pack('<q', field.n_y))
    file_phase.write(struct.pack('<q', field.n_y))

    i = point_on_x
    k = point_on_z
    for j, y in enumerate(field.y):
            file_amp.write(struct.pack('<dd', y, field.amp(i, j, k)))
            file_phase.write(struct.pack('<dd', y, field.phase(i, j, k)))

    file_amp.close()
    file_phase.close()

    logging.info("Binary Field Plane output has been finished")



def PrintCylAmpPhaseZ_Binary(field, point_on_x, point_on_y):
    logging.info("Binary Field Plane output has been started")

    file_amp = open('amp_pres_z.bin', 'wb')
    file_phase = open('phase_pres_z.bin', 'wb')

    file_amp.write(struct.pack('<q', field.n_z))
    file_phase.write(struct.pack('<q', field.n_z))

    i = point_on_x
    j = point_on_y
    for k, z in enumerate(field.z):
            file_amp.write(struct.pack('<dd', z, field.amp(i, j, k)))
            file_phase.write(struct.pack('<dd', z, field.phase(i, j, k)))

    file_amp.close()
    file_phase.close()

    logging.info("Binary Field Plane output has been finished")
