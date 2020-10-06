# -*- coding: utf-8 -*-

"""
Copyright 2012. Sergey Ilin.
Lab 366, Acoustic Department, Faculty of Physics
Lomonosov Moscow State University

"""


import logging
import codecs


class Medium:
    def __init__(self, name, density, speed_of_sound, attenuation=0.0, attenuation_dB_m=0.0):
        self.name = name
        self.density = density
        self.speed_of_sound = speed_of_sound
        if attenuation_dB_m != 0.0:
            self.attenuation = 0.115129254 * attenuation_dB_m
        else:
            self.attenuation = attenuation
        logging.info('Medium has been set: ' + str(self))

    def __str__(self):
        return 'Medium "{}" has density = {} kg/m3, speed of sound = {} m/s, attenuation coef = {} Np/m'.format(
            str(self.name),
            str(self.density),
            str(self.speed_of_sound),
            str(self.attenuation)
            )


def medium_from_file(file_path):
    # Adding new medium from file
    # 'utf-8-sig' means that it skips BOM in UTF files created in Notepad
    with codecs.open(file_path, encoding='utf-8-sig') as f:
        lines = [line.strip() for line in f]

    if len(lines) < 3:
        logging.info("Couldn't read array parameters from file " + file_path)
        return

    medium = Medium(
        name = lines[0],
        density = float(lines[1]),
        speed_of_sound = float(lines[2])
        )
    return medium


if __name__ == '__main__':
    # test case
    # medium = Medium("Water", 1000.0, 1488.0, attenuation_dB_m = 10.0)
    medium = medium_from_file("water_medium.txt")
    print(medium)
