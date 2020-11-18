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


    def nodes_overall_conn(self, make_symmetric=True):
        '''

        computing the overall connectivity of each node
        regardless of network affiliation

        Parameters
        ----------

        make_symmetric: Boolean|
                        True indicate that the matrix is either upper
                        or lower triangular and need to be symmetrize
                        False indicate that the matrix is a full matrix already

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
            if make_symmetric==True:
                self.matrix = self.matrix + self.matrix.T - np.diag(self.matrix.diagonal())
            else:
                self.matrix = self.matrix
            np.fill_diagonal(self.matrix,0)
            for nodes in range(self.matrix.shape[0]):
                self._node_conn =  np.sum(self.matrix[nodes])
                self.nodes_conn.append(self._node_conn)
        self.nodes_conn = np.array(self.nodes_conn)
        self.nodes_conn = self.nodes_conn.reshape(len(self.matrices_files), self.matrix.shape[0])

        return self.nodes_conn



    def node_inner_conn(self, sbj_number, nodes_number, make_symmetric=True):
        '''
        computing the connectivity of each node with its own network

        Parameters
        ----------
        sbj_number: int |
                    number of subjects
        nodes_number: int|
                    number of nodes
        make_symmetric: Boolean|
                        True indicate that the matrix is either upper
                        or lower triangular and need to be symmetrize
                        False indicate that the matrix is a full matrix already

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
            if make_symmetric==True:
                self.matrix = self.matrix + self.matrix.T - np.diag(self.matrix.diagonal())
            else:
                self.matrix = self.matrix
            np.fill_diagonal(self.matrix,0)
            for network in net:
                for nodes in self.labels_dic[network]:
                    self.sub_matrix =self.matrix[nodes]
                    self.streamlines_sum = np.sum(self.sub_matrix[self.labels_dic[network]])
                    self.all_conn[subj, nodes] = self.streamlines_sum/self.labels_dic[network].shape[0]

        return self.all_conn



    def node_outer_conn(self, sbj_number, nodes_number, make_symmetric=True):
        '''
        computing the connectivity of each node with the other nodes
        which don't belong to the same network

        Parameters
        ----------
        sbj_number: int |
                    number of subjects
        nodes_number: int|
                    number of nodes
        make_symmetric: Boolean|
                        True indicate that the matrix is either upper
                        or lower triangular and need to be symmetrize
                        False indicate that the matrix is a full matrix already

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
            if make_symmetric==True:
                self.matrix = self.matrix + self.matrix.T - np.diag(self.matrix.diagonal())
            else:
                self.matrix = self.matrix
            np.fill_diagonal(self.matrix,0)
            self.nodes_ranges = np.arange(len(self.labels_dic['nodes']))
            for network in net:
                self.outer_idx = np.setdiff1d(self.nodes_ranges, self.labels_dic[network])
                for nodes in self.outer_idx:
                    self.sub_matrix =self.matrix[nodes]
                    self.streamlines_sum = np.sum(self.sub_matrix[self.outer_idx])
                    self.all_conn[subj, nodes] = self.streamlines_sum/self.outer_idx.shape[0]

        return self.all_conn



    def node_ranking(self, sbj_number, nodes_number, networks_number, make_symmetric=True):
        '''
        computing how much each node is connected with the each network

        Parameters
        ----------

        sbj_number: int |
                    number of subjects
        nodes_number: int|
                    number of nodes
        networks_number: int|
                    number of networks
        make_symmetric: Boolean|
                        True indicate that the matrix is either upper
                        or lower triangular and need to be symmetrize
                        False indicate that the matrix is a full matrix already
        Returns
        -------

        float data : numpy array |
        numpy a 3D array (dim number of subject X  number of network X number of network)
        representing the connectivity of each node with all the networks

        '''
        with open(self.net_label_txt) as f:
            net=f.read().splitlines()
        self.all_conn = np.zeros([sbj_number, nodes_number, networks_number])
        for subj in range(len(self.matrices_files)):
            self.matrix = pd.read_csv(self.matrices_files[subj], sep= ' ', header=None)
            self.matrix = np.array(self.matrix)
            if make_symmetric==True:
                self.matrix = self.matrix + self.matrix.T - np.diag(self.matrix.diagonal())
            else:
                self.matrix = self.matrix
            np.fill_diagonal(self.matrix,0)
            for nodes in range(self.matrix.shape[0]):
                self.node_conn = self.matrix[nodes]
                for network in net:
                    self.streamlines_sum =  np.sum(self.node_conn[self.labels_dic[network]])
                    self.all_conn[subj, nodes, net.index(network)] = self.streamlines_sum/self.labels_dic[network].shape[0]

        return self.all_conn



    def net_inner_conn(self, make_symmetric=True):
        '''
        computing the how much each network is connected with itself

        Parameters
        ----------

        make_symmetric: Boolean|
                        True indicate that the matrix is either upper
                        or lower triangular and need to be symmetrize
                        False indicate that the matrix is a full matrix already
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
            if make_symmetric==True:
                self.matrix = self.matrix + self.matrix.T - np.diag(self.matrix.diagonal())
            else:
                self.matrix = self.matrix
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



    def net_outer_conn(self, make_symmetric=True):
        '''
        computing how much each network is connected with the other
        networks

        Parameters
        ----------

        make_symmetric: Boolean|
                        True indicate that the matrix is either upper
                        or lower triangular and need to be symmetrize
                        False indicate that the matrix is a full matrix already

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
            if make_symmetric==True:
                self.matrix = self.matrix + self.matrix.T - np.diag(self.matrix.diagonal())
            else:
                self.matrix = self.matrix
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



    def net_ranking(self, sbj_number, nodes_number, make_symmetric=True, percentage_value=False):
        '''
        computing how much each node is connected with the each network

        Parameters
        ----------

        sbj_number: int |
                    number of subjects
        nodes_number: int|
                    number of nodes
        make_symmetric: Boolean|
                        True indicate that the matrix is either upper
                        or lower triangular and need to be symmetrize
                        False indicate that the matrix is a full matrix already
        percentage_value: Boolean|
                        True return values express in percentage_value
                        False return raw values

        
        Returns
        -------

        float data : numpy array |
        numpy a 3D array (dim number of subject X  number of network X number of network)
        representing the connectivity of each node with all the networks

        '''
        with open(self.net_label_txt) as f:
            net=f.read().splitlines()
        self.all_conn = self.node_ranking(sbj_number, nodes_number, len(net), make_symmetric=make_symmetric)
        self.all_conn_rank = np.zeros([sbj_number, len(net), len(net)])
        for subj in range(len(self.matrices_files)):
            self.subj2use = self.all_conn[subj,:,:]
            for network in net:
                self.net2use = self.subj2use[self.labels_dic[network],:]
                if percentage_value==False:
                    self.all_conn_rank[subj, net.index(network), :] = np.mean(self.net2use, axis=0)
                else:
                    self.all_conn_rank[subj, net.index(network), :] = 100* np.mean(self.net2use, axis=0)/np.sum(np.mean(self.net2use, axis=0))

        return self.all_conn_rank



    def all_standard_metrics(self, sbj_number, nodes_number, networks_number, make_symmetric=True, percentage_value=False):
        self.metrics_dict = {
        "nodes_overall_conn": self.nodes_overall_conn(make_symmetric=make_symmetric),
        "node_inner_conn": self.node_inner_conn(sbj_number, nodes_number, make_symmetric=make_symmetric),
        "node_outer_conn": self.node_outer_conn(sbj_number, nodes_number, make_symmetric=make_symmetric),
        "node_ranking": self.node_ranking(sbj_number, nodes_number, networks_number, make_symmetric=make_symmetric),
        "net_inner_conn": self.net_inner_conn(make_symmetric=make_symmetric),
        "net_outer_conn": self.net_outer_conn(make_symmetric=make_symmetric),
        "net_ranking": self.net_ranking(sbj_number, nodes_number, make_symmetric=make_symmetric, percentage_value=percentage_value)
        }

        return self.metrics_dict




class Graph_Theory(object):

    def __init__(self, matrices_files, net_label_txt, labels_dic):

        self.matrices_files = matrices_files
        self.net_label_txt = net_label_txt
        self.labels_dic = labels_dic

    def nodal_degree(self, sbj_number, nodes_number, make_symmetric=True):
        '''
        computing how much each node is connected with the each network


        Parameters
        ----------
        sbj_number: int |
                    number of subjects
        nodes_number: int|
                      number of nodes
        make_symmetric: Boolean|
                        True indicate that the matrix is either upper
                        or lower triangular and need to be symmetrize
                        False indicate that the matrix is a full matrix already

        Returns
        -------

        dict: : dictonary with the following keys |

        degree: int | Number of links connected to the node
                    The indegree
        in_degree: int | Number of inward links
        out_degree: int | Number of outward links
        joint_in_degree : int | number of vertices with in_degree>out_degree
        joint_out_degree : int | number of vertices with out_degree>in_degree
        joint_bilateral : int | number of vertices with in_degree==out_degree
        '''


        self.all_nodal_degree = {
        "degree": np.zeros([sbj_number, nodes_number]),
        "in_degree" : np.zeros([sbj_number, nodes_number]),
        "out_degree" : np.zeros([sbj_number, nodes_number]),
        "joint_in_degree" : np.zeros([sbj_number, nodes_number]),
        "joint_out_degree" : np.zeros([sbj_number, nodes_number]),
        "joint_bilateral" : np.zeros([sbj_number, nodes_number])
        }
        for subj in range(len(self.matrices_files)):
            self.matrix = pd.read_csv(self.matrices_files[subj], sep= ' ', header=None)
            self.matrix = np.array(self.matrix)
            if make_symmetric==True:
                self.matrix = self.matrix + self.matrix.T - np.diag(self.matrix.diagonal())
            else:
                self.matrix = self.matrix
            np.fill_diagonal(self.matrix,0)
            self.matrix_bin = self.matrix
            self.matrix_bin[self.matrix_bin != 0] = 1 #bin matrix
            self.indegree = np.sum(self.matrix_bin, axis=0)  # indegree = column sum of CIJ
            self.all_nodal_degree['in_degree'][subj] = self.indegree
            self.outdegree = np.sum(self.matrix_bin, axis=1)  # outdegree = row sum of CIJ
            self.all_nodal_degree['out_degree'][subj] = self.outdegree
            self.deg = self.indegree + self.outdegree  # degree = indegree+outdegree
            self.all_nodal_degree['degree'][subj] = self.deg

            self.szJ = np.max((self.indegree, self.outdegree)) + 1
            self.joint = np.zeros((self.szJ, self.szJ))

            for i in range(len(self.matrix_bin)):
                self.joint[self.indegree[i], self.outdegree[i]] += 1

            self.joint_od = np.sum(np.triu(self.joint, 1))
            self.joint_id = np.sum(np.tril(self.joint, -1))
            self.joint_bl = np.sum(np.diag(self.joint))
            self.all_nodal_degree['joint_in_degree'][subj] = self.joint_id
            self.all_nodal_degree['joint_out_degree'][subj] = self.joint_od
            self.all_nodal_degree['joint_bilateral'][subj] = self.joint_bl

        return self.all_nodal_degree
