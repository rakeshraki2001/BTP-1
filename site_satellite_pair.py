
import csv
import math
from math import cos
from math import sin
from math import radians as rad
import numpy as np

import matplotlib.pyplot as plt

from satellite_and_site import Satellite
from satellite_and_site import Site


class Pairs(Site, Satellite):
    all = []
    valid = []
    pi = 3.1415926535
    re = 6371.0
    atm_boundary_height = 2.0
    ratm = re + atm_boundary_height
    diT = 0.2
    dgR = 2
    wavelength = 737*pow(10, -9)
    alpha = 0.028125

    def __init__(self, site: Site, satellite: Satellite):
        self.silat = site.lat
        self.silon = site.lon
        self.salat = satellite.lat
        self.salon = satellite.lon
        self.sahei = satellite.height
        Pairs.all.append(self)

    def central_angle(self):
        return math.acos(cos(rad(self.salat))*cos(rad(self.silat))*cos(rad(self.salon-self.silon)) + sin(rad(self.salat))*sin(rad(self.silat)))

    def central_angle_indegrees(self):
        return (180/self.pi)*self.central_angle()

    def visibility_test(self):
        re = Pairs.re
        rs = Pairs.re + self.sahei
        return self.central_angle() <= math.acos(re/rs)

    def elevation_angle(self):
        rs = Pairs.re + self.sahei
        if self.visibility_test():
            return math.acos(rs*sin(self.central_angle())/self.dis_sat2site())
        else:
            return -1

    def elevation_angle_indegrees(self):
        return (180/self.pi)*self.elevation_angle()

    def dis_sat2site(self):
        re = Pairs.re
        rs = Pairs.re+self.sahei
        if self.visibility_test():
            return rs*(math.sqrt(1+math.pow(re/rs, 2)-2*(re/rs)*cos(self.central_angle())))
        else:
            return -1

    def dis_site_atmboundary(self):
        re = Pairs.re
        ratm = Pairs.ratm
        a = 90 + self.elevation_angle_indegrees()
        b = (180/self.pi)*math.asin(re*sin(rad(a))/ratm)
        c = 180 - a - b
        return math.sqrt(pow(re,2) + pow(ratm,2) - 2*re*ratm*cos(rad(c)))

    def nig_freespace(self):
        return (pow(self.pi,2)*pow(self.diT,2)*pow(self.dgR,2)/16)/pow(self.wavelength*self.dis_sat2site()*1000,2)

    def nig_atmboundary(self):
        return pow(math.e, -self.alpha*self.dis_site_atmboundary())

    def nig(self):
        return self.nig_freespace()*self.nig_atmboundary()

    def possible(self, threshold = 0.1):
        return self.visibility_test() and (self.nig() >= threshold)

    @classmethod
    def create_pairs(cls,height=500.0,threshold = 0.1):
        Site.create_sites_from_csv()
        Satellite.create_satellites_float(10, 10, height)
        for site in Site.all:
            for satellite in Satellite.all:
                pai = Pairs(site, satellite)
                if pai.possible(threshold):
                    cls.valid.append(pai)

    def __repr__(self):
        return f'Pair: site[{self.silat}, {self.silon}] satellite[{self.salat}, {self.salon}, {self.sahei}]'


#pa = Pairs(Site(0,0),Satellite(0,30,1000))
#print(pa.visibility_test())
#print(pa.central_angle())
#print(pa.central_angle_indegrees())
#print(pa.dis_sat2site())
#print(pa.elevation_angle())
#print(pa.elevation_angle_indegrees())
#print(pa.dis_site_atmboundary())
#print(pa.nig_atmboundary())
#print(pa.nig_freespace())
#print(pa.nig())


def plot_nig_variation(height, unvisible_nig=-1):
    lat = np.outer(np.linspace(0, 90, 90, endpoint=False), np.ones(90))
    lon = lat.copy().T

    z = np.outer(np.linspace(0, 90, 90, endpoint=False), np.ones(90))

    for h in height:
        for i in range(0,90):
            for j in range(0,90):
                pai = Pairs(Site(0,0),Satellite(i,j,h))
                if pai.visibility_test():
                    z[i][j] = pai.nig()
                else:
                    z[i][j] = unvisible_nig

        fig = plt.figure(figsize =(14, 9))
        ax = plt.axes(projection ='3d')
        # print(len(pai.all))
        ax.plot_surface(lat, lon, z)

        plt.show()


hei = [500,600,700,800]
bighei = [700,1000,5000,10000,20000,25000,30000,35000]

if __name__ == '__main_':
    Pairs.create_pairs(25000, 0)
    print(Site.all)
    print(Satellite.all)
    print(len(Satellite.all))
    print(len(Pairs.all))
    print(len(Pairs.valid))
    #plot_nig_variation(hei)
    plot_nig_variation(bighei,0)


def percentage_of_valid_pairs(num: int, noflong: int, height=800.0, threshold=0.1):
    Site.create_sites_from_csv()
    Satellite.create_satellites_float(num, noflong, height)
    #print('nofsatellites', len(Satellite.all))
    #print('nofsites', len(Site.all))
    total, count = 0, 0
    for sat in Satellite.all:
        for site in Site.all:
            # print(sat,site,Pairs(site,sat).possible(0))
            total += 1
            if Pairs(site, sat).possible(threshold):
                count += 1
    #print(total//2, count//2, (count / total) * 100/2)
    Site.all.clear()
    Satellite.all.clear()
    return total//2, count//2, (count / total) * 100/2


def percentage_of_valid_conns(num: int, noflong: int, height=800.0, threshold=0.1):
    Site.create_sites_from_csv()
    Satellite.create_satellites_float(num, noflong, height)
    #print('nofsatellites', len(Satellite.all))
    #print('nofsites', len(Site.all))
    total, count = 0, 0
    for sat in Satellite.all:
        for s1 in Site.all:
            for s2 in Site.all:
                if s1 != s2:
                    # print(sat,site,Pairs(site,sat).possible(0))
                    total += 1
                    if Pairs(s1, sat).possible(threshold) and Pairs(s2, sat).possible(threshold):
                        count += 1
    #print(total//2, count//2, (count / total) * 100/2)
    Site.all.clear()
    Satellite.all.clear()
    return total//2, count//2, (count / total) * 100/2


def percentage_of_valid_conns_per_satellite(num: int, noflong: int, height=800.0, threshold=0.1):
    conns_per_sat = []
    Site.create_sites_from_csv()
    Satellite.create_satellites_float(num, noflong, height)
    #print('nofsatellites', len(Satellite.all))
    #print('nofsites', len(Site.all))
    for sat in Satellite.all:
        total, count = 0, 0
        for s1 in Site.all:
            for s2 in Site.all:
                if s1 != s2:
                    # print(sat,site,Pairs(site,sat).possible(0))
                    total += 1
                    if Pairs(s1, sat).possible(threshold) and Pairs(s2, sat).possible(threshold):
                        count += 1
        #print(total, count, (count / total) * 100)
        conns_per_sat.append([total//2, count//2, (count / total) * 100/2])
        Site.all.clear()
        Satellite.all.clear()
    return conns_per_sat


if __name__ == '__main_':
    #print(percentage_of_valid_pairs(num=100, noflong=100, height=800, threshold=0))
    print(percentage_of_valid_conns(num=10, noflong=10, height=35000, threshold=0.00001))
    #lis = percentage_of_valid_conns_per_satellite(num=10, noflong=10, height=800, threshold=0)
    #print(lis)
    #print(len(lis))
    # lis = percentage_of_valid_conns_per_satellite(num=100, noflong=100, height=800, threshold=0)
    # print(len(lis))

if __name__ == '__main__':
    height = []
    valid_conn = []
    for h in range(500, 35000, 100):
        height += [h]
        valid_conn += [percentage_of_valid_conns(num=10, noflong=10, height=h, threshold=0)[1]]
        Site.all.clear()
        Satellite.all.clear()
    print(valid_conn)
    plt.plot(height, valid_conn)
    # naming the x axis
    plt.xlabel('height')
    # naming the y axis
    plt.ylabel('no_of_valid_connections')
    # giving a title to my graph
    plt.title('NofConnections Vs HEIGHT')
    # function to show the plot
    plt.show()