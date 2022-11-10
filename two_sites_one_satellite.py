import csv
import math
from math import cos
from math import sin
from math import radians as rad
import numpy as np

import matplotlib.pyplot as plt
from scipy.optimize import linear_sum_assignment
from satellite_and_site import Satellite
from satellite_and_site import Site
from site_satellite_pair import Pairs
import networkx as nx
import random
import WBbM as WBbM


class Connection(Pairs, Site, Satellite):
    valid = []
    greedy_sat_all = []
    greedy_sat_n_recievers_all = []
    greedy_sat_inf_recievers_all = []
    greedy_sat_n_recievers_maxnigmul_all = []
    greedy_sat_inf_recievers_maxnigmul_all = []
    greedy_sat_n_recievers_maxnigadd_all = []
    greedy_sat_inf_recievers_maxnigadd_all = []
    greedy_sat_nig_all = []
    threshold = 0.1
    ebits_rate = pow(10, 9)

    def __init__(self, s1: Site, s2: Site, sat: Satellite):
        self.site1 = s1
        self.site2 = s2
        self.satellite = sat
        #if self.possible_conn():
            # Connection.valid.append(self)

    def possible_conn(self):
        return (Pairs(self.site1, self.satellite)).possible(Connection.threshold) and (Pairs(self.site2, self.satellite)).possible(Connection.threshold)

    @staticmethod
    def greedy_sat(num: int, noflong: int, height=500.0, threshold=0.01):
        Site.create_sites_from_csv()
        Satellite.create_satellites_float(num, noflong, height)
        nofsites = len(Site.all)
        nofsat = len(Satellite.all)
        # print(nofsat, nofsites)
        for sat in Satellite.all:
            for i, s1 in enumerate(Site.all):
                # print(Pairs(s1, sat).possible(Connection.threshold))
                if Pairs(s1, sat).possible(threshold) and not s1.connected and not sat.connected:
                    for j, s2 in enumerate(Site.all[i+1:], i+1):
                        if Pairs(s2, sat).possible(threshold) and not s2.connected and not sat.connected:
                            #print(Pairs(s1, sat).possible(threshold), Pairs(s2, sat).possible(threshold))
                            s1.connected = True
                            s2.connected = True
                            sat.connected = True
                            Connection.greedy_sat_all.append(Connection(s1,s2,sat))

    @staticmethod
    def greedy_sat_nig(num: int, noflong: int, height=500.0, threshold=0.01):
        Site.create_sites_from_csv()
        Satellite.create_satellites_float(num, noflong, height)
        nofsites = len(Site.all)
        nofsat = len(Satellite.all)
        # print(nofsat, nofsites)
        total = 0
        for sat in Satellite.all:
            site1 = site2 = Site.all[0]
            satellite = Satellite.all[0]
            found = False
            maxnig = 0
            for i, s1 in enumerate(Site.all):
                # print(Pairs(s1, sat).possible(Connection.threshold))
                if Pairs(s1, sat).possible(threshold) and not s1.connected and not sat.connected:
                    for j, s2 in enumerate(Site.all[i + 1:], i + 1):
                        if Pairs(s2, sat).possible(threshold) and not s2.connected and not sat.connected:
                            if Pairs(s1, sat).nig() + Pairs(s2, sat).nig() > maxnig:
                                maxnig = Pairs(s1, sat).nig() + Pairs(s2, sat).nig()
                                site1 = s1
                                site2 = s2
                                satellite = sat
                                found = True

            if found:
                Connection.greedy_sat_nig_all.append(Connection(site1, site2, satellite))
                total += maxnig * Connection.ebits_rate
                site1.connected = True
                site2.connected = True
                satellite.connected = True
        return total
    @staticmethod
    def greedy_sat_n_recievers(num: int, noflong: int, height=500.0, threshold=0.01, nofreci=5):
        Site.recievers_per_site = nofreci
        Site.create_sites_from_csv()
        # print(Site.all)
        Satellite.create_satellites_float(num, noflong, height)
        nofsites = len(Site.all)
        nofsat = len(Satellite.all)
        # print(nofsat, nofsites)
        for sat in Satellite.all:
            found = False
            for i, s1 in enumerate(Site.all):
                if found:
                    break
                # print(Pairs(s1, sat).possible(Connection.threshold))
                if Pairs(s1, sat).possible(threshold) and (s1.recievers != 0):
                    for j, s2 in enumerate(Site.all[i + 1:], i + 1):
                        if Pairs(s2, sat).possible(threshold) and (s2.recievers != 0):
                            s1.recievers -= 1
                            s2.recievers -= 1
                            sat.connected = True
                            #print(Pairs(s1, sat).possible(threshold),Pairs(s2, sat).possible(threshold))
                            Connection.greedy_sat_n_recievers_all.append(Connection(s1, s2, sat))
                            found = True
                            break

    @staticmethod
    def greedy_sat_infinite_recievers(num: int, noflong: int, height=500.0, threshold=0.01):
        Site.create_sites_from_csv()
        Satellite.create_satellites_float(num, noflong, height)
        nofsites = len(Site.all)
        nofsat = len(Satellite.all)
        # print(nofsat, nofsites)
        for sat in Satellite.all:
            found = False
            for i, s1 in enumerate(Site.all):
                if found:
                    break
                # print(Pairs(s1, sat).possible(Connection.threshold))
                if Pairs(s1, sat).possible(threshold) and not sat.connected:
                    for j, s2 in enumerate(Site.all[i + 1:], i + 1):
                        if Pairs(s2, sat).possible(threshold) and not sat.connected:
                            # print(Pairs(s1, sat).possible(threshold),Pairs(s2, sat).possible(threshold))
                            Connection.greedy_sat_inf_recievers_all.append(Connection(s1, s2, sat))
                            sat.connected = True
                            found = True
                            break

    @staticmethod
    def greedy_sat_n_recievers_maxnigmul(num: int, noflong: int, height=500.0, threshold=0.01, nofreci=5):
        Site.recievers_per_site = nofreci
        Site.create_sites_from_csv()
        # print(Site.all)
        Satellite.create_satellites_float(num, noflong, height)
        nofsites = len(Site.all)
        nofsat = len(Satellite.all)
        # print(nofsat, nofsites)
        total = 0
        for sat in Satellite.all:
            maxnig = 0
            site1 = site2 = Site.all[0]
            satellite = Satellite.all[0]
            found = False
            for i, s1 in enumerate(Site.all):
                # print(Pairs(s1, sat).possible(Connection.threshold))
                if Pairs(s1, sat).possible(threshold) and (s1.recievers != 0):
                    for j, s2 in enumerate(Site.all[i + 1:], i + 1):
                        if Pairs(s2, sat).possible(threshold) and (s2.recievers != 0):
                            if not found:
                                s1.recievers -= 1
                                s2.recievers -= 1
                            if Pairs(s1,sat).nig()*Pairs(s2,sat).nig()>maxnig:
                                maxnig = Pairs(s1,sat).nig()*Pairs(s2,sat).nig()
                                site1 = s1
                                site2 = s2
                                satellite = sat
                            # print(Pairs(s1, sat).possible(threshold),Pairs(s2, sat).possible(threshold))
                            found = True
            if found:
                satellite.connected = True
                total += (Pairs(site1, satellite).nig() + Pairs(site2, satellite).nig())*Connection.ebits_rate
                Connection.greedy_sat_n_recievers_maxnigmul_all.append(Connection(site1, site2, satellite))
        return total

    @staticmethod
    def greedy_sat_infinite_recievers_maxnigmul(num: int, noflong: int, height=500.0, threshold=0.01):
        Site.create_sites_from_csv()
        # print(Site.all)
        Satellite.create_satellites_float(num, noflong, height)
        nofsites = len(Site.all)
        nofsat = len(Satellite.all)
        # print(nofsat, nofsites)
        total = 0
        for sat in Satellite.all:
            maxnig = 0
            site1 = site2 = Site.all[0]
            satellite = Satellite.all[0]
            found = False
            for i, s1 in enumerate(Site.all):
                # print(Pairs(s1, sat).possible(Connection.threshold))
                if Pairs(s1, sat).possible(threshold):
                    for j, s2 in enumerate(Site.all[i + 1:], i + 1):
                        if Pairs(s2, sat).possible(threshold):
                            if Pairs(s1, sat).nig() * Pairs(s2, sat).nig() > maxnig:
                                maxnig = Pairs(s1, sat).nig() * Pairs(s2, sat).nig()
                                site1 = s1
                                site2 = s2
                                satellite = sat
                            # print(Pairs(s1, sat).possible(threshold),Pairs(s2, sat).possible(threshold))
                            found = True
            if found:
                satellite.connected = True
                total += (Pairs(site1, satellite).nig() + Pairs(site2, satellite).nig()) * Connection.ebits_rate
                Connection.greedy_sat_inf_recievers_maxnigmul_all.append(Connection(site1, site2, satellite))
        return total

    @staticmethod
    def greedy_sat_n_recievers_maxnigadd(num: int, noflong: int, height=500.0, threshold=0.01, nofreci=5):
        Site.recievers_per_site = nofreci
        Site.create_sites_from_csv()
        # print(Site.all)
        Satellite.create_satellites_float(num, noflong, height)
        nofsites = len(Site.all)
        nofsat = len(Satellite.all)
        # print(nofsat, nofsites)
        total = 0
        for sat in Satellite.all:
            maxnig = 0
            site1 = site2 = Site.all[0]
            satellite = Satellite.all[0]
            found = False
            for i, s1 in enumerate(Site.all):
                # print(Pairs(s1, sat).possible(Connection.threshold))
                if Pairs(s1, sat).possible(threshold) and (s1.recievers != 0):
                    for j, s2 in enumerate(Site.all[i + 1:], i + 1):
                        if Pairs(s2, sat).possible(threshold) and (s2.recievers != 0):
                            if not found:
                                s1.recievers -= 1
                                s2.recievers -= 1
                            if Pairs(s1, sat).nig() + Pairs(s2, sat).nig() > maxnig:
                                maxnig = Pairs(s1, sat).nig() + Pairs(s2, sat).nig()
                                site1 = s1
                                site2 = s2
                                satellite = sat
                            # print(Pairs(s1, sat).possible(threshold),Pairs(s2, sat).possible(threshold))
                            found = True
            if found:
                satellite.connected = True
                total += maxnig * Connection.ebits_rate
                Connection.greedy_sat_n_recievers_maxnigadd_all.append(Connection(site1, site2, satellite))
        return total

    @staticmethod
    def greedy_sat_infinite_recievers_maxnigadd(num: int, noflong: int, height=500.0, threshold=0.01):
        Site.create_sites_from_csv()
        # print(Site.all)
        Satellite.create_satellites_float(num, noflong, height)
        nofsites = len(Site.all)
        nofsat = len(Satellite.all)
        # print(nofsat, nofsites)
        total = 0
        for sat in Satellite.all:
            maxnig = 0
            site1 = site2 = Site.all[0]
            satellite = Satellite.all[0]
            found = False
            for i, s1 in enumerate(Site.all):
                # print(Pairs(s1, sat).possible(Connection.threshold))
                if Pairs(s1, sat).possible(threshold):
                    for j, s2 in enumerate(Site.all[i + 1:], i + 1):
                        if Pairs(s2, sat).possible(threshold):
                            if Pairs(s1, sat).nig() + Pairs(s2, sat).nig() > maxnig:
                                maxnig = Pairs(s1, sat).nig() + Pairs(s2, sat).nig()
                                site1 = s1
                                site2 = s2
                                satellite = sat
                            # print(Pairs(s1, sat).possible(threshold),Pairs(s2, sat).possible(threshold))
                            found = True
            if found:
                satellite.connected = True
                total += maxnig * Connection.ebits_rate
                Connection.greedy_sat_inf_recievers_maxnigadd_all.append(Connection(site1, site2, satellite))
        return total

    @staticmethod
    def hungarian_sat_matching(num: int, noflong: int, height=500.0, threshold=0.01):
        Site.create_sites_from_csv()
        Satellite.create_satellites_float(num, noflong, height)
        nofsites = len(Site.all)
        nofsat = len(Satellite.all)
        # print(nofsat, nofsites)
        nig_satellites_two_sites = []
        total = 0

        for i, s1 in enumerate(Site.all):
            for j, s2 in enumerate(Site.all[i + 1:], i + 1):
                nig_satellite = []
                for sat in Satellite.all:
                    if Pairs(s1, sat).possible(threshold) and Pairs(s2, sat).possible(threshold):
                        # print(Pairs(s1, sat).possible(threshold),Pairs(s2, sat).possible(threshold))
                        nig_satellite += [Pairs(s1, sat).nig() + Pairs(s2, sat).nig()]
                        total += 1
                    else:
                        nig_satellite += [0]
                nig_satellites_two_sites += [nig_satellite]
        nig_satellites_two_sites = np.array(nig_satellites_two_sites)
        nofreps = int(nofsat//(nofsites*(nofsites-1)/2))+1
        #print(nofreps)
        #print(nofsat)
        #print(np.shape(np.tile(nig_satellites_two_sites, (nofreps,1))))
        #nig_satellites_two_sites = np.tile(nig_satellites_two_sites, (nofreps,1))
        row_ind, col_ind = linear_sum_assignment(nig_satellites_two_sites, maximize=True)
        Site.all.clear()
        Satellite.all.clear()
        return Connection.ebits_rate * nig_satellites_two_sites[row_ind, col_ind].sum()

    @staticmethod
    def hungarian_sat_full_matching(num: int, noflong: int, hei=500.0, thr=0.01):
        selected_edges = list()
        # ----------------------------------------------------------------------
        # Create a bipartite complete graph with random weights between 0 and 9.
        #num, noflong, height, th = 8, 9, 15000, 0.00001
        g = nx.Graph()
        Site.create_sites_from_csv()
        Satellite.create_satellites_float(num, noflong, hei)

        site_pair_all = []
        sat_all = []

        nsite = 0
        for i, s1 in enumerate(Site.all):
            for j, s2 in enumerate(Site.all[i + 1:], i + 1):
                g.add_node(nsite, data=[s1, s2])
                site_pair_all += [nsite]
                nsite += 1
                # print("Sites-"+str(i)+str(j), end=' ')

        tnsite = nsite

        nsat = 0 + nsite
        for i, sat in enumerate(Satellite.all):
            g.add_node(nsat, data=sat)
            sat_all += [nsat]
            nsat += 1
            # print("Satellite-"+str(i), end=' ')

        # print(site_pair_all)
        # print(sat_all)

        nsite = 0
        for i, s1 in enumerate(Site.all):
            for j, s2 in enumerate(Site.all[i + 1:], i + 1):
                nsat = 0 + tnsite
                for k, sat in enumerate(Satellite.all, 100):
                    if Pairs(s1, sat).possible(thr) and Pairs(s2, sat).possible(thr):
                        wei = (Pairs(s1, sat).nig() + Pairs(s2, sat).nig()) * pow(10, 9)
                        g.add_edge(nsite, nsat, weight=wei)
                    else:
                        g.add_edge(nsite, nsat, weight=0)
                    nsat += 1
                nsite += 1

        # print(g.number_of_edges())
        # """
        W = list()

        right = set(site_pair_all)
        left = set(sat_all)

        num_left, num_right = len(left), len(right)

        for node in left:
            node_weights = list()
            node_edges = g.edges(node, data=True)
            # print(node_edges)
            for node_edge in node_edges:
                node_weights.append(node_edge[2]['weight'])
            W.append(node_weights)

        # print(W)
        row_capacity, column_capacity = [1] * len(W), [1] * len(list(zip(*W)))
        # minimum paper cardinality
        ldp = 0
        # maximum paper cardinality with the minimum value of 2
        # can be constant or a list of all capacities.
        udp = len(Satellite.all)
        # udp = row_capacity
        # maximum papers one reviewer will review with the minimum value of 2
        # can be constant or a list of all capacities.
        uda = 1
        # uda = column_capacity
        # minimum papers every reviewer has to review
        lda = 0
        # ----------------------------------------------------------------------
        # Solve the bipartite b-matching problem with the WBbM algorithm.
        b_matching = WBbM.WBbM(num_left, num_right, [j for j in list(np.concatenate(W))], lda, uda, ldp, udp,
                               LogToConsole=0)
        results, total_weight = b_matching.Bb_matching(optimization_mode="max")

        for row_index in range(len(results)):
            for column_index in range(len(results[row_index])):
                if results[row_index][column_index] == 1:
                    selected_edges.append(
                        (list(right)[column_index], list(left)[row_index]))  # the order based on the gold-standard

        #print("Selected edges are:", selected_edges, "\nTotal weight:", total_weight)
        # ----------------------------------------------------------------------
        # Illustrate the matching output.
        # elarge = [(u, v) for (u, v, d) in g.edges(data=True) if d['weight'] > threshold]
        # esmall = [(u, v) for (u, v, d) in g.edges(data=True) if d['weight'] <= threshold]

        # plot_graph(g, elarge, esmall, "red")  # edges with a weight over the threshold
        # plot_graph(g, selected_edges, g.edges - selected_edges, "green")  # selected edges based on the b-matching algorithm
        Site.all.clear()
        Satellite.all.clear()
        return total_weight

    def __repr__(self):
        return f'Connection({self.site1}, {self.site2}, {self.satellite})'

if __name__ == '__main_':
    h,th,nr = 35000,0.00001,50
    Connection.greedy_sat(10, 10, h, th)
    #print(Connection.greedy_sat_all)
    print(len(Connection.greedy_sat_all))
    Site.all.clear()
    Satellite.all.clear()
    Connection.greedy_sat_all.clear()
    Connection.greedy_sat_n_recievers(10, 10, h, th, nr)
    #print(Connection.greedy_sat_n_recievers_all)
    print(len(Connection.greedy_sat_n_recievers_all))
    Site.all.clear()
    Satellite.all.clear()
    Connection.greedy_sat_n_recievers_all.clear()
    Connection.greedy_sat_infinite_recievers(10, 10, h, th)
    #print(Connection.greedy_sat_inf_recievers_all)
    print(len(Connection.greedy_sat_inf_recievers_all))
    Site.all.clear()
    Satellite.all.clear()
    Connection.greedy_sat_inf_recievers_all.clear()
    print(Connection.greedy_sat_n_recievers_maxnigmul(10,10,h,th,1))
    Site.all.clear()
    Satellite.all.clear()
    Connection.greedy_sat_n_recievers_maxnigmul_all.clear()
    print(Connection.greedy_sat_n_recievers_maxnigmul(10, 10, h, th, nr))
    Site.all.clear()
    Satellite.all.clear()
    Connection.greedy_sat_n_recievers_maxnigmul_all.clear()
    print(Connection.greedy_sat_infinite_recievers_maxnigmul(10,10,h,th))

if __name__ == '__main_':
    height = [i for i in range(500,35000,100)]
    valid_conn_inf_rec = []
    for h in height:
        Connection.greedy_sat_infinite_recievers(10, 10, h, 0.00001)
        valid_conn_inf_rec += [len(Connection.greedy_sat_inf_recievers_all)]
        Site.all.clear()
        Satellite.all.clear()
        Connection.greedy_sat_inf_recievers_all.clear()
    plt.plot(height, valid_conn_inf_rec)
    # naming the x axis
    plt.xlabel('height')
    # naming the y axis
    plt.ylabel('no_of_valid_connections with infinite recievers')
    # giving a title to my graph
    plt.title('no_of_valid_connections with infinite recievers Vs HEIGHT)')
    # function to show the plot
    plt.show()

if __name__ == '__main_':
    height = [500,1000,5000,10000,25000,35000]
    #nofrec = [1,2,3,5,10,20,30,40,50]
    nofrec = [i for i in range(1,60)]
    for h in height:
        valid_conn = []
        for n in nofrec:
            Connection.greedy_sat_n_recievers(10, 10, h, 0.00001, n)
            #print(Connection.greedy_sat_n_recievers_all)
            #print(len(Connection.greedy_sat_n_recievers_all), n)
            valid_conn += [len(Connection.greedy_sat_n_recievers_all)]
            Site.all.clear()
            Satellite.all.clear()
            Connection.greedy_sat_n_recievers_all.clear()
        #print(valid_conn)
        plt.plot(nofrec, valid_conn)
        # naming the x axis
        plt.xlabel('nofrecievers')
        # naming the y axis
        plt.ylabel('no_of_valid_connections')
        # giving a title to my graph
        plt.title(f'NofConnections Vs nofreceivers at HEIGHT(h={h})')
        # function to show the plot
        plt.show()


if __name__ == '__main_':
    height = [500,1000,5000,10000,25000,35000]
    #threshold = [0.01,0.01,0.001,0.0001,0.00001,0.00001]
    threshold = [0.00001,0.00001,0.00001,0.00001,0.00001,0.00001]
    #nofrec = [1,2,3,5,10,20,30,40,50]
    nofrec = [i for i in range(1,60)]
    for h, th in zip(height, threshold):
        nof_valid_conn = []
        total_nig_valid_conn = []
        for n in nofrec:
            total_nig_valid_conn += [Connection.greedy_sat_n_recievers_maxnigmul(10, 10, h, th, n)]
            #print(Connection.greedy_sat_n_recievers_all)
            #print(len(Connection.greedy_sat_n_recievers_all), n)
            nof_valid_conn += [len(Connection.greedy_sat_n_recievers_maxnigmul_all)]
            Site.all.clear()
            Satellite.all.clear()
            Connection.greedy_sat_n_recievers_maxnigmul_all.clear()
        #print(valid_conn)
        fig, ax = plt.subplots()

        # Plot linear sequence, and set tick labels to the same color
        ax.plot(nofrec, total_nig_valid_conn, color='red', label='total_nig_valid_conn')
        ax.tick_params(axis='y', labelcolor='red')

        # Generate a new Axes instance, on the twin-X axes (same position)
        ax2 = ax.twinx()

        # Plot exponential sequence, set scale to logarithmic and change tick color
        ax2.plot(nofrec, nof_valid_conn, color='green', label='nof_valid_conn')
        ax2.tick_params(axis='y', labelcolor='green')
        ax.set_xlabel("nofrecievers")
        ax.set_ylabel("total_nig_valid_conn")
        ax2.set_ylabel("nofvalid_conn")
        plt.title(f'total_max_nig_multiplication and Nofvalidconns Vs nofreceivers at HEIGHT(h={h})')
        plt.show()

if __name__ == '__main_':
    height = [500,1000,5000,10000,25000,35000]
    #threshold = [0.01,0.01,0.001,0.0001,0.00001,0.00001]
    threshold = [0.00001,0.00001,0.00001,0.00001,0.00001,0.00001]
    #nofrec = [1,2,3,5,10,20,30,40,50]
    nofrec = [i for i in range(1,60)]
    for h, th in zip(height, threshold):
        nof_valid_conn = []
        total_nig_valid_conn = []
        for n in nofrec:
            total_nig_valid_conn += [Connection.greedy_sat_n_recievers_maxnigadd(10, 10, h, th, n)]
            #print(Connection.greedy_sat_n_recievers_all)
            #print(len(Connection.greedy_sat_n_recievers_all), n)
            nof_valid_conn += [len(Connection.greedy_sat_n_recievers_maxnigadd_all)]
            Site.all.clear()
            Satellite.all.clear()
            Connection.greedy_sat_n_recievers_maxnigadd_all.clear()
        #print(valid_conn)
        fig, ax = plt.subplots()

        # Plot linear sequence, and set tick labels to the same color
        ax.plot(nofrec, total_nig_valid_conn, color='red', label='total_nig_valid_conn')
        ax.tick_params(axis='y', labelcolor='red')

        # Generate a new Axes instance, on the twin-X axes (same position)
        ax2 = ax.twinx()

        # Plot exponential sequence, set scale to logarithmic and change tick color
        ax2.plot(nofrec, nof_valid_conn, color='green', label='nof_valid_conn')
        ax2.tick_params(axis='y', labelcolor='green')
        ax.set_xlabel("nofrecievers")
        ax.set_ylabel("total_nig_valid_conn")
        ax2.set_ylabel("nofvalid_conn")
        plt.title(f'total_max_nig_addition and Nofvalidconns Vs nofreceivers at HEIGHT(h={h})')
        plt.show()

if __name__ == '__main_':
    height = [10000]#,1000,5000,10000,25000,35000]
    #threshold = [0.01,0.01,0.001,0.0001,0.00001,0.00001]
    threshold = [0.00001]#,0.00001,0.00001,0.00001,0.00001,0.00001]
    #nofrec = [1,2,3,5,10,20,30,40,50]
    nofrec = [i for i in range(1,60)]
    for h, th in zip(height, threshold):
        total_nig_valid_conn_mul = []
        total_nig_valid_conn_add = []
        for n in nofrec:
            total_nig_valid_conn_mul += [Connection.greedy_sat_n_recievers_maxnigmul(10, 10, h, th, n)]
            #print(Connection.greedy_sat_n_recievers_all)
            #print(len(Connection.greedy_sat_n_recievers_all), n)
            Site.all.clear()
            Satellite.all.clear()
            Connection.greedy_sat_n_recievers_maxnigmul_all.clear()
            total_nig_valid_conn_add += [Connection.greedy_sat_n_recievers_maxnigadd(10, 10, h, th, n)]
            # print(Connection.greedy_sat_n_recievers_all)
            # print(len(Connection.greedy_sat_n_recievers_all), n)
            Site.all.clear()
            Satellite.all.clear()
            Connection.greedy_sat_n_recievers_maxnigadd_all.clear()

        #print(valid_conn)
        fig, ax = plt.subplots()

        # Plot linear sequence, and set tick labels to the same color
        ax.plot(nofrec, total_nig_valid_conn_mul, color='red', label='total_nig_valid_conn_mul')
        ax.tick_params(axis='y', labelcolor='red')

        # Generate a new Axes instance, on the twin-X axes (same position)
        ax2 = ax.twinx()

        # Plot exponential sequence, set scale to logarithmic and change tick color
        ax2.plot(nofrec, total_nig_valid_conn_add, color='green', label='total_nig_valid_conn_add')
        ax2.tick_params(axis='y', labelcolor='green')
        ax.set_xlabel("nofrecievers")
        ax.set_ylabel("total_nig_valid_conn_mul")
        ax2.set_ylabel("total_nig_valid_conn_add")
        plt.title(f'total_max_nig_multiplication and total_max_nig_addition Vs nofreceivers at HEIGHT(h={h})')
        plt.show()

if __name__ == '__main__':
    height = [i for i in range(500,35000,100)]
    total_nig_add = []
    total_nig_mul = []
    diff = []
    for h in height:
        mul = Connection.greedy_sat_infinite_recievers_maxnigmul(10, 10, h, 0.00001)
        Site.all.clear()
        Satellite.all.clear()
        Connection.greedy_sat_inf_recievers_maxnigmul_all.clear()
        add = Connection.greedy_sat_infinite_recievers_maxnigadd(10, 10, h, 0.00001)
        Site.all.clear()
        Satellite.all.clear()
        Connection.greedy_sat_inf_recievers_maxnigadd_all.clear()
        total_nig_add += [add]
        total_nig_mul += [mul]
        diff += [mul - add]
        #print(mul-add, end='  ')
    fig, ax = plt.subplots()

    # Plot linear sequence, and set tick labels to the same color
    ax.plot(height, total_nig_mul, color='red', label='total_nig_valid_conn_mul')
    ax.tick_params(axis='y', labelcolor='red')

    # Generate a new Axes instance, on the twin-X axes (same position)
    ax2 = ax.twinx()

    # Plot exponential sequence, set scale to logarithmic and change tick color
    ax2.plot(height, total_nig_add, color='green', label='total_nig_valid_conn_add')
    ax2.tick_params(axis='y', labelcolor='green')
    ax.set_xlabel("height")
    ax.set_ylabel("total_nig_mul")
    ax2.set_ylabel("total_nig_add")

    ax3 = ax.twinx()
    ax3.spines.right.set_position(("axes", 1.07))
    # Plot exponential sequence, set scale to logarithmic and change tick color
    ax3.plot(height, diff, color='orange', label='Difference')
    ax3.tick_params(axis='y', labelcolor='orange')
    # ax.set_xlabel("height")
    # ax.set_ylabel("total_nig_mul")
    ax3.set_ylabel("Difference between these two")

    plt.title(f'total_nig_multiplication and total_nig_addition with infinite receivers Vs HEIGHT)')
    plt.show()

if __name__ == '__main_':
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
    plt.show()

if __name__ == '__main_':
    height = [i for i in range(500,35000,100)]
    total_nig_sat = []
    total_nig_hungarian = []
    for h in height:
        res = Connection.greedy_sat_nig(10, 10, h, 0.00001)
        Site.all.clear()
        Satellite.all.clear()
        Connection.greedy_sat_nig_all.clear()
        total_nig_sat += [res]
        total_nig_hungarian += [Connection.hungarian_sat_matching(10, 10, h, 0.00001)]
    fig, ax = plt.subplots()

    # Plot linear sequence, and set tick labels to the same color
    ax.plot(height, total_nig_sat, color='red', label='total_nig_sat')
    ax.tick_params(axis='y', labelcolor='red')

    # Generate a new Axes instance, on the twin-X axes (same position)
    ax2 = ax.twinx()

    # Plot exponential sequence, set scale to logarithmic and change tick color
    ax2.plot(height, total_nig_hungarian, color='green', label='total_nig_hungarian')
    ax2.tick_params(axis='y', labelcolor='green')
    ax.set_xlabel("height")
    ax.set_ylabel("total_nig_sat")
    ax2.set_ylabel("total_nig_hungarian")
    plt.title(f'total_nig_sat and total_nig_hungarian_method with infinite receivers Vs HEIGHT)')
    plt.show()

if __name__ == '__main_':
    height = [i for i in range(500, 35000, 100)]
    total_nig_sat_inf = []
    total_nig_hungarian = []
    diff = []
    for h in height:
        res = Connection.greedy_sat_infinite_recievers_maxnigadd(8, 9, h, 0.00001)
        Site.all.clear()
        Satellite.all.clear()
        Connection.greedy_sat_inf_recievers_maxnigadd_all.clear()
        total_nig_sat_inf += [res]
        total_nig_hungarian += [Connection.hungarian_sat_full_matching(8, 9, h, 0.00001)]
        diff += [total_nig_hungarian[-1] - total_nig_sat_inf[-1]]
    fig, ax = plt.subplots()

    # Plot linear sequence, and set tick labels to the same color
    ax.plot(height, total_nig_sat_inf, color='red', label='total_nig_sat_inf')
    ax.tick_params(axis='y', labelcolor='red')

    # Generate a new Axes instance, on the twin-X axes (same position)
    ax2 = ax.twinx()

    # Plot exponential sequence, set scale to logarithmic and change tick color
    ax2.plot(height, total_nig_hungarian, color='green', label='total_nig_hungarian')
    ax2.tick_params(axis='y', labelcolor='green')
    ax.set_xlabel("height")
    ax.set_ylabel("total_nig_sat_with_inf_recievers")
    ax2.set_ylabel("total_nig_full_hungarian")

    ax3 = ax.twinx()
    ax3.spines.right.set_position(("axes", 1.07))
    # Plot exponential sequence, set scale to logarithmic and change tick color
    ax3.plot(height, diff, color='orange', label='Difference')
    ax3.tick_params(axis='y', labelcolor='orange')
    # ax.set_xlabel("height")
    # ax.set_ylabel("total_nig_mul")
    ax3.set_ylabel("Difference between these two")

    plt.title(f'total_nig_sat_inf and total_nig_full_hungarian_method with infinite receivers Vs HEIGHT)')
    plt.show()