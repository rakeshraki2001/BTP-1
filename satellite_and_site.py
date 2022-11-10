import csv
import math
from math import cos
from math import sin
from math import radians as rad
import numpy as np

import matplotlib.pyplot as plt


class Point:
    def __init__(self, la, lo):
        assert ((la <= 90) and (la >= -90)), f"Invalid Latitude"
        assert ((lo <= 180) and (lo > -180)), f"Invalid Longitude"
        self.lat = la
        self.lon = lo


class Satellite(Point):
    all = []

    def __init__(self, la: float, lo: float, height=500.0):
        assert height > 0, f"Invalid height"
        self.height = height
        Point.__init__(self,la,lo)
        self.connected = False
        Satellite.all.append(self)

    @staticmethod
    def create_satellites(num: int, noflong: int, height=500.0):
        for la in range(0, 90, int(90/((num+1)//2))):
            for lo in range(180, -180, -360//noflong):
                Satellite(la, lo, height)

        if num%2 == 0:
            num = num + 2

        for la in range(0, -90, -int(90/((num+1)//2))):
            if la == 0:
                continue
            for lo in range(180, -180, -360//noflong):
                Satellite(la, lo, height)

    @staticmethod
    def create_satellites_float(num: int, noflong: int, height=500.0):
        la = 0.0
        while la <= 90.0:
            lon = 180.0
            while lon > -180.0:
                # print(la,lon)
                Satellite(la, lon, height)
                lon -= 360.0/noflong
            la += 90.0/((num+1)//2)

        if num%2 == 0:
            num = num + 2

        la = -90.0
        while la < 0.0:
            lon = 180.0
            while lon > -180.0:
                Satellite(la, lon, height)
                lon -= 360.0 / noflong
            la += 90.0 / ((num + 1) // 2)

    def __repr__(self):
        return f"Satellite({self.lat}, {self.lon}, {self.height})"


class Site(Point):
    all = []
    recievers_per_site = 5

    def __init__(self, lat: float, lon: float):
        super().__init__(
            lat, lon
        )
        self.recievers = Site.recievers_per_site
        self.connected = False
        Site.all.append(self)

    @classmethod
    def create_sites_from_csv(cls, filename='sites.csv', change_in_lon=0):
        with open(filename, 'r') as f:
            reader = csv.DictReader(f)
            sites = list(reader)
            for site in sites:
                lat = float(site.get('lat'))
                lon = float(site.get('lon'))
                if lon+change_in_lon > 180 and change_in_lon != 0:
                    lon = -180 + change_in_lon - (180-lon)
                else:
                    lon = lon + change_in_lon
                Site(lat, lon)

    def __repr__(self):
        return f"Site({self.lat}, {self.lon}, {self.recievers})"


if __name__ == '__main__':
    Site.create_sites_from_csv(change_in_lon = 0)
    print(Site.all)
    print(len(Site.all))




