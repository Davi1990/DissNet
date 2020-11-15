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
        '''

        computing the overall connectivity of each node
        regardless of network affiliation


        Returns
        -------

        float data : numpy array |
        numpy array (dim number of subject X number of node)
        representing the connectivity of each node regardless
        of network affiliation

        '''

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


    def node_inner_conn(self, sbj_number, nodes_number):
        '''
        computing the connectivity of each node with its own network

        Parameters
        ----------
        sbj_number: int |
        number of subjects
        nodes_number: int|
        number of nodes


        Returns
        -------

        float data : numpy array |
        numpy array (dim number of subject X number of node)
        representing the connectivity of each node with its own
        network

        '''

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
                    self.streamlines_sum = np.sum(self.sub_matrix[self.labels_dic[network]])
                    self.all_conn[subj, nodes] = self.streamlines_sum/self.labels_dic[network].shape[0]
        return self.all_conn


    def node_outer_conn(self, sbj_number, nodes_number):
        '''
        computing the connectivity of each node with the other nodes
        which don't belong to the same network

        Parameters
        ----------
        sbj_number: int |
        number of subjects
        nodes_number: int|
        number of nodes


        Returns
        -------

        float data : numpy array |
        numpy array (dim number of subject X number of node)
        representing the connectivity of each node with regions that
        are outsite the node's network

        '''

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
                    self.streamlines_sum = np.sum(self.sub_matrix[self.outer_idx])
                    self.all_conn[subj, nodes] = self.streamlines_sum/self.outer_idx.shape[0]
        return self.all_conn


    def net_inner_conn(self):
        '''
        computing the how much each network is connected with itself

        Returns
        -------

        float data : numpy array |
        numpy array (dim number of subject X number of network)
        representing the connectivity of each network with itself

        '''

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


    def net_outer_conn(self):
        '''
        computing the how much each network is connected with the other
        networks

        Returns
        -------

        float data : numpy array |
        numpy array (dim number of subject X number of network)
        representing the connectivity of each network with other networks

        '''


        with open(self.net_label_txt) as f:
            net=f.read().splitlines()
        self.all_conn = []
        for subj in range(len(self.matrices_files)):
            self.matrix = pd.read_csv(self.matrices_files[subj], sep= ' ', header=None)
            self.matrix = np.array(self.matrix)
            np.fill_diagonal(self.matrix,0)
            self.nodes_ranges = np.arange(len(self.labels_dic['nodes']))
            for network in net:
                self.outer_idx = np.setdiff1d(self.nodes_ranges, self.labels_dic[network])
                self.subj_matrix = self.matrix[self.labels_dic[network]]
                self.subj_matrix = self.subj_matrix[:,self.outer_idx]
                self.streamlines_sum = np.sum(np.sum(self.subj_matrix))
                self.conn_measure = self.streamlines_sum/self.outer_idx.shape[0]
                self.all_conn.append(self.conn_measure)
        self.all_conn = np.array(self.all_conn)
        self.all_conn = self.all_conn.reshape(len(self.matrices_files), len(net))
        return self.all_conn
