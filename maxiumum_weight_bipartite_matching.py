import csv
import math
from math import cos
from math import sin
from math import radians as rad
import numpy as np
import scipy
import matplotlib.pyplot as plt

from satellite_and_site import Satellite
from satellite_and_site import Site
from site_satellite_pair import Pairs
from site_satellite_pair import percentage_of_valid_conns
import two_sites_one_satellite
from two_sites_one_satellite import *
from scipy.optimize import linear_sum_assignment


class Connection:
    ebits_rate = pow(10, 9)

    @staticmethod
    def hungarian_sat_matching(num: int, noflong: int, height=500.0, threshold=0.01):
        Site.create_sites_from_csv()
        Satellite.create_satellites_float(num, noflong, height)
        nofsites = len(Site.all)
        nofsat = len(Satellite.all)
        # print(nofsat, nofsites)
        nig_satellites_two_sites = []
        total = 0
        for sat in Satellite.all:
            nig_satellite = []
            for i, s1 in enumerate(Site.all):
                for j, s2 in enumerate(Site.all[i + 1:], i + 1):
                    if Pairs(s1, sat).possible(threshold) and Pairs(s2, sat).possible(threshold):
                        # print(Pairs(s1, sat).possible(threshold),Pairs(s2, sat).possible(threshold))
                        nig_satellite += [Pairs(s1, sat).nig()+Pairs(s2, sat).nig()]
                        total += 1
                    else:
                        nig_satellite += [0]
            nig_satellites_two_sites += [nig_satellite]
        nig_satellites_two_sites = np.array(nig_satellites_two_sites)
        row_ind, col_ind = linear_sum_assignment(nig_satellites_two_sites, maximize=True)
        return Connection.ebits_rate * nig_satellites_two_sites[row_ind, col_ind].sum()



if __name__ == '__main__':
    height = [i for i in range(500,35000,100)]
    total_nig_mul = []
    nof_valid_conn_inf_rec = []
    for h in height:
        mul = Connection.greedy_sat_infinite_recievers_maxnigmul(10, 10, h, 0.00001)
        Site.all.clear()
        Satellite.all.clear()
        Connection.greedy_sat_inf_recievers_maxnigmul_all.clear()
        Connection.greedy_sat_infinite_recievers(10, 10, h, 0.00001)
        nof_valid_conn_inf_rec += [len(Connection.greedy_sat_inf_recievers_all)]
        Site.all.clear()
        Satellite.all.clear()
        Connection.greedy_sat_inf_recievers_all.clear()
        total_nig_mul += [mul]
        #print(mul-add, end='  ')
    fig, ax = plt.subplots()

    # Plot linear sequence, and set tick labels to the same color
    ax.plot(height, total_nig_mul, color='red', label='total_nig_mul')
    ax.tick_params(axis='y', labelcolor='red')

    # Generate a new Axes instance, on the twin-X axes (same position)
    ax2 = ax.twinx()

    # Plot exponential sequence, set scale to logarithmic and change tick color
    ax2.plot(height, nof_valid_conn_inf_rec, color='green', label='total_nig_valid_conn_add')
    ax2.tick_params(axis='y', labelcolor='green')
    ax.set_xlabel("height")
    ax.set_ylabel("total_nig_mul")
    ax2.set_ylabel("nof_valid_conn_inf_rec")
    plt.title(f'total_nig_multiplication and nof_valid_conn_inf_rec with infinite receivers Vs HEIGHT)')
    plt.show()

