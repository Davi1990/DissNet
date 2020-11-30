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
import bct as bct


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
        computing how much each node is connected with the other network

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

    def nodal_degree(self, sbj_number, nodes_number, make_symmetric=True, binarize=False):
        '''
        computing graph theory node measures regardless of network affiliation


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
        binarize: Boolean|
                        True will make the connectivity matrix binary
                        Default is False

        Returns
        -------

        dict: : dictonary with the following keys |

        degree: int | Number of links connected to the node
        in_degree: int | Number of inward links
        out_degree: int | Number of outward links
        joint_in_degree: int | number of vertices with in_degree>out_degree
        joint_out_degree: int | number of vertices with out_degree>in_degree
        joint_bilateral: int | number of vertices with in_degree==out_degree
        node_strength_dir: int | node strength (in-strength + out-strength)
        node_strength_undir: int | sum of weights of links connected to the node
        '''


        self.all_nodal_degree = {
        "degree": np.zeros([sbj_number, nodes_number]),
#        "in_degree" : np.zeros([sbj_number, nodes_number]),
#        "out_degree" : np.zeros([sbj_number, nodes_number]),
#        "joint_in_degree" : np.zeros([sbj_number, nodes_number]),
#        "joint_out_degree" : np.zeros([sbj_number, nodes_number]),
#        "joint_bilateral" : np.zeros([sbj_number, nodes_number]),
#        "node_strength_dir": np.zeros([sbj_number, nodes_number]),
        "node_strength_undir":np.zeros([sbj_number, nodes_number])
        }
        for subj in range(len(self.matrices_files)):
            self.matrix = pd.read_csv(self.matrices_files[subj], sep= ' ', header=None)
            self.matrix = np.array(self.matrix)
            if make_symmetric==True:
                self.matrix = self.matrix + self.matrix.T - np.diag(self.matrix.diagonal())
            else:
                self.matrix = self.matrix

            if binarize==True:
                self.matrix = bct.algorithms.binarize(self.matrix)
            else:
                self.matrix = self.matrix

            np.fill_diagonal(self.matrix,0)

            self.deg = bct.algorithms.degrees_und(self.matrix)
            # self.all_nodal_degree['in_degree'][subj] = self.inp
            # self.all_nodal_degree['out_degree'][subj] = self.od
            self.all_nodal_degree['degree'][subj] = self.deg

            # self.J, self.J_od, self.J_id, self.J_bl = bct.algorithms.jdegree(self.matrix)
            # self.all_nodal_degree['joint_in_degree'][subj] = self.J_id
            # self.all_nodal_degree['joint_out_degree'][subj] = self.J_od
            # self.all_nodal_degree['joint_bilateral'][subj] = self.J_bl

#            self.nodestr_dir = bct.algorithms.strengths_dir(self.matrix)
#            self.all_nodal_degree['node_strength_dir'][subj] = self.nodestr_dir

            self.nodestr_undir = bct.algorithms.strengths_und(self.matrix)
            self.all_nodal_degree['node_strength_undir'][subj] = self.nodestr_undir

        return self.all_nodal_degree


    def network_level_degree(self, sbj_number, nodes_number, label_dic, make_symmetric=True, binarize=False):
        '''
        computing graph theory node measures specific for each network


        Parameters
        ----------
        sbj_number: int |
                    number of subjects
        nodes_number: int|
                      number of nodes
        label_dic: dict |
                    dictonary computed using files.labels()
        make_symmetric: Boolean|
                        True indicate that the matrix is either upper
                        or lower triangular and need to be symmetrize
                        False indicate that the matrix is a full matrix already
        binarize: Boolean|
                        True will make the connectivity matrix binary
                        Default is False

        Returns
        -------

        dict: : dictonary with the following keys |

        degree: int | Number of links connected to the node
        in_degree: int | Number of inward links
        out_degree: int | Number of outward links
        joint_in_degree: int | number of vertices with in_degree>out_degree
        joint_out_degree: int | number of vertices with out_degree>in_degree
        joint_bilateral: int | number of vertices with in_degree==out_degree
        node_strength_dir: int | node strength (in-strength + out-strength)
        node_strength_undir: int | sum of weights of links connected to the node
        '''

        with open(self.net_label_txt) as f:
            net=f.read().splitlines()
        self.degree = self.nodal_degree(sbj_number, nodes_number, make_symmetric=make_symmetric, binarize=binarize)
        self.values = np.zeros([sbj_number, len(self.degree.keys()), len(net)])
        self.list = list(self.degree.keys())

        for subject in range(sbj_number):
             for key in self.list:
                 for network in net:
                     self.values[subject, self.list.index(key), net.index(network)] = np.mean(self.degree[key][subject][label_dic[network]])

        self.d = {}
        for i in self.degree.keys():
            self.d[i] = self.values[:, self.list.index(i), :]

        return self.d



    def physical_connectivity(self, sbj_number, networks_number, label_dic, make_symmetric=True, binarize=False):
        '''
        Density is the fraction of present connections to possible connections.


        Parameters
        ----------
        sbj_number: int |
                    number of subjects
        networks_number: int|
                      number of networks
        label_dic: dict |
                    dictonary computed using files.labels()
        make_symmetric: Boolean|
                        True indicate that the matrix is either upper
                        or lower triangular and need to be symmetrize
                        False indicate that the matrix is a full matrix already

        binarize= Boolean|
                        True will make the connectivity matrix binary
                        Default is False


        Returns
        -------

        dict: : dictonary with the following keys |

        Density_und: int | Density is the fraction of present connections
                        to possible connections

        '''

        with open(self.net_label_txt) as f:
            net=f.read().splitlines()

        self.physical_connectivity = {
        "Density_und": np.zeros([sbj_number, networks_number]),
        }

        for subj in range(len(self.matrices_files)):
            self.matrix = pd.read_csv(self.matrices_files[subj], sep= ' ', header=None)
            self.matrix = np.array(self.matrix)
            if make_symmetric==True:
                self.matrix = self.matrix + self.matrix.T - np.diag(self.matrix.diagonal())
            else:
                self.matrix = self.matrix

            if binarize==True:
                self.matrix = bct.algorithms.binarize(self.matrix)
            else:
                self.matrix = self.matrix

            np.fill_diagonal(self.matrix,0)

            for network in net:
                self.net_matrix = self.matrix[label_dic[network]]
                self.net_matrix = self.net_matrix[:,label_dic[network]]
                self.kden, self.n, self.k = bct.algorithms.density_und(self.net_matrix)
                self.physical_connectivity['Density_und'][subj, net.index(network)] = self.kden

        return self.physical_connectivity




    def modularity(self, sbj_number, networks_number, label_dic, make_symmetric=True, binarize=False):
        '''
        Computing modularity of the adjencency matrix adjacency matrix


        Parameters
        ----------
        sbj_number: int |
                    number of subjects
        networks_number: int|
                      number of networks
        label_dic:  dict |
                    dictonary computed using files.labels()
        make_symmetric: Boolean|
                        True indicate that the matrix is either upper
                        or lower triangular and need to be symmetrize
                        False indicate that the matrix is a full matrix already
        binarize= Boolean|
                        True will make the connectivity matrix binary
                        Default is False


        Returns
        -------

        dict: : dictonary with the following keys |

        community_louvain: int | Modularity values
        similarity_idx: float32 | Values indicating how much each original
                                 network is similar to the new modules found
                                 with the modularity algorithm
        '''


        with open(self.net_label_txt) as f:
            net=f.read().splitlines()

        self.modularity = {
        "community_louvain": np.zeros([sbj_number, networks_number]),
        "similarity_idx": np.zeros([sbj_number, networks_number])
        }

        for subj in range(len(self.matrices_files)):
            self.matrix = pd.read_csv(self.matrices_files[subj], sep= ' ', header=None)
            self.matrix = np.array(self.matrix)
            if make_symmetric==True:
                self.matrix = self.matrix + self.matrix.T - np.diag(self.matrix.diagonal())
            else:
                self.matrix = self.matrix

            if binarize==True:
                self.matrix = bct.algorithms.binarize(self.matrix)
            else:
                self.matrix = self.matrix
            np.fill_diagonal(self.matrix,0)

            for network in net:
                self.net_matrix = self.matrix[label_dic[network]]
                self.net_matrix = self.net_matrix[:,label_dic[network]]
                self.ci, self.q  = bct.algorithms.community_louvain(self.net_matrix)
                self.modularity['community_louvain'][subj, net.index(network)] = self.q

            self.unico = np.unique(self.ci)
            self.index={}
            for values in self.unico:
                self.index[values]= np.where(self.ci == values)
            self.similarity_matrix = np.zeros([len(net), self.unico.shape[0]])
            for network in net:
                for module in self.unico:
                    self.simil =  np.intersect1d(label_dic[network], self.index[module])
                    self.similarity_matrix[net.index(network), module-1] = self.simil.shape[0]

                self.modularity['similarity_idx'][subj, net.index(network)] = (np.max(self.similarity_matrix[net.index(network)]) - np.min(self.similarity_matrix[net.index(network)])) / label_dic[network].shape[0]

        return self.modularity




    def centrality(self, sbj_number, nodes_number, label_dic, atlas, make_symmetric=True, binarize=False):

        with open(self.net_label_txt) as f:
            net=f.read().splitlines()

        self.atlas = pd.read_excel(atlas, header=None)
        self.atlas = np.array(self.atlas)
        self.ci_original = self.atlas[:,8]

        self.centrality = {
        "edge_betweeness_bin": np.zeros([sbj_number, nodes_number]),
        "edge_betweeness_wei": np.zeros([sbj_number, nodes_number]),
        "eigenvector_centrality_und": np.zeros([sbj_number, nodes_number]),
        "coreness_kcoreness_centrality_bu": np.zeros([sbj_number, nodes_number]),
        "kn_kcoreness_centrality_bu": np.zeros([sbj_number, nodes_number]),
        "module_degree_zscore": np.zeros([sbj_number, nodes_number]),
        "participation_coef": np.zeros([sbj_number, nodes_number]),
        "subgraph_centrality": np.zeros([sbj_number, nodes_number])
        }

        for subj in range(len(self.matrices_files)):
            self.matrix = pd.read_csv(self.matrices_files[subj], sep= ' ', header=None)
            self.matrix = np.array(self.matrix)
            if make_symmetric==True:
                self.matrix = self.matrix + self.matrix.T - np.diag(self.matrix.diagonal())
            else:
                self.matrix = self.matrix

            self.matrix_bin = bct.algorithms.binarize(self.matrix)
            self.matrix_weight = self.matrix

            if binarize==True:
                self.matrix = bct.algorithms.binarize(self.matrix)
            else:
                self.matrix = self.matrix


            np.fill_diagonal(self.matrix,0)
            np.fill_diagonal(self.matrix_bin,0)
            np.fill_diagonal(self.matrix_weight,0)

            self.BC = bct.betweenness_bin(self.matrix_bin)
            self.centrality['edge_betweeness_bin'][subj] = self.BC

            self.BC_w = bct.betweenness_wei(self.matrix_weight)
            self.centrality['edge_betweeness_wei'][subj] = self.BC_w

            self.v = bct.eigenvector_centrality_und(self.matrix)
            self.centrality['eigenvector_centrality_und'][subj] = self.v

            self.coreness, self.kn =  bct.kcoreness_centrality_bu(self.matrix_bin)
            self.centrality['coreness_kcoreness_centrality_bu'][subj] = self.coreness
            self.centrality['kn_kcoreness_centrality_bu'][subj] = self.kn

            self.Z =  bct.module_degree_zscore(self.matrix, ci=self.ci_original)
            self.centrality['module_degree_zscore'][subj] = self.Z

            self.P = bct.participation_coef(self.matrix, ci=self.ci_original)
            self.centrality['participation_coef'][subj] = self.P

            self.Cs = bct.subgraph_centrality(self.matrix_bin)
            self.centrality['subgraph_centrality'][subj] = self.Cs

        return self.centrality
