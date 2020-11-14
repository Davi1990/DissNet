"""
miscelallaneous functions and classes to extract connectivity metrics
Author: Davide Momi, PhD [momi.davide89@gmail.com], https://twitter.com/davemomi
"""

import numpy as np
import pandas as pd
from math import pi
import glob
import seaborn as sns
import matplotlib.pyplot as plt


class Connectivity_metrics(object):

    def __init__(self, matrices_files, net_label_txt, labels_dic):

        self.matrices_files = matrices_files
        self.net_label_txt = net_label_txt
        self.labels_dic = labels_dic

    def nodes_overall_conn(self):
        self.nodes_conn = []
        for subj in range(len(self.matrices_files)):
            self.matrix = pd.read_csv(self.matrices_files[subj], sep= ' ', header=None)
            self.matrix = np.array(self.matrix)
            np.fill_diagonal(self.matrix,0)
            for nodes in range(self.matrix.shape[0]):
                self._node_conn =  np.sum(self.matrix[nodes])
                self.nodes_conn.append(self._node_conn)
        self.nodes_conn = np.array(self.nodes_conn)
        self.nodes_conn = self.nodes_conn.reshape(len(self.matrices_files), self.matrix.shape[0])
        return self.nodes_conn


    def node_inner_net_conn(self, sbj_number, nodes_number):
        with open(self.net_label_txt) as f:
            net=f.read().splitlines()
        self.all_conn = np.zeros([sbj_number, nodes_number])
        for subj in range(len(self.matrices_files)):
            self.matrix = pd.read_csv(self.matrices_files[subj], sep= ' ', header=None)
            self.matrix = np.array(self.matrix)
            np.fill_diagonal(self.matrix,0)
            for network in net:
                for nodes in self.labels_dic[network]:
                    self.sub_matrix =self.matrix[nodes]
                    self.all_conn[subj, nodes] = np.sum(self.sub_matrix[self.labels_dic[network]])
        return self.all_conn


    def node_outer_net_conn(self, sbj_number, nodes_number):
        with open(self.net_label_txt) as f:
            net=f.read().splitlines()
        self.all_conn = np.zeros([sbj_number, nodes_number])
        for subj in range(len(self.matrices_files)):
            self.matrix = pd.read_csv(self.matrices_files[subj], sep= ' ', header=None)
            self.matrix = np.array(self.matrix)
            np.fill_diagonal(self.matrix,0)
            self.nodes_ranges = np.arange(len(self.labels_dic['nodes']))
            for network in net:
                self.outer_idx = np.setdiff1d(self.nodes_ranges, self.labels_dic[network])
                for nodes in self.outer_idx:
                    self.sub_matrix =self.matrix[nodes]
                    self.all_conn[subj, nodes] = np.sum(self.sub_matrix[self.outer_idx])
        return self.all_conn


    def net_inner_conn(self):
        with open(self.net_label_txt) as f:
            net=f.read().splitlines()
        self.all_conn = []
        for subj in range(len(self.matrices_files)):
            self.matrix = pd.read_csv(self.matrices_files[subj], sep= ' ', header=None)
            self.matrix = np.array(self.matrix)
            np.fill_diagonal(self.matrix,0)
            for network in net:
                self.subj_matrix = self.matrix[self.labels_dic[network]]
                self.subj_matrix = self.subj_matrix[:,self.labels_dic[network]]
                self.streamlines_sum = np.sum(np.sum(self.subj_matrix))
                self.conn_measure = self.streamlines_sum/len(self.labels_dic[network])
                self.all_conn.append(self.conn_measure)
        self.all_conn = np.array(self.all_conn)
        self.all_conn = self.all_conn.reshape(len(self.matrices_files), len(net))
        return self.all_conn
