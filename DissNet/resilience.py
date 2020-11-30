import networkx, numpy, operator, pylab, random, sys
import bct as bct
import pandas as pd
import numpy as np

class Resilience(object):

    def __init__(self, matrices_files, net_label_txt, labels_dic):

        self.matrices_files = matrices_files
        self.net_label_txt = net_label_txt
        self.labels_dic = labels_dic

    def nodal_degree_vulnerability(self, sbj_number, nodes_number, make_symmetric=True, binarize=False,
                                   recalculate = False, attack_type='target'):
        '''
        Performs robustness analysis based on nodal degree.


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
        recalculate: Boolean|
                        It will use sequential (recalculate = True) or
                        simultaneous (recalculate = False) approach.
                        Default is False

        Returns
        -------

        vulnerability: np.array |
                    The overall vulnerability of the network

        '''

        self.all_vulnerability = np.zeros([sbj_number])

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
            if attack_type == 'target':
                self.deg = bct.algorithms.degrees_und(self.matrix)
                self.g = networkx.convert_matrix.from_numpy_array(self.matrix)
                self.m = dict(enumerate(self.deg.flatten(), 0))
                self.l = sorted(self.m.items(), key = operator.itemgetter(1), reverse = True)
                self.x = []
                self.y = []
                self.largest_component = max(networkx.connected_components(self.g), key = len)
                self.n = len(self.g.nodes())
                self.x.append(0)
                self.y.append(len(self.largest_component) * 1. / self.n)
                self.R = 0.0
                for i in range(1, self.n - 1):
                    self.g.remove_node(self.l.pop(0)[0])
                    if recalculate:
                        self.matrix = networkx.convert_matrix.to_numpy_array(self.g)
                        self.g = networkx.convert_matrix.from_numpy_array(self.matrix)
                        self.deg = bct.algorithms.degrees_und(self.matrix)
                        self.m = dict(enumerate(self.deg.flatten(), 0))
                        self.l = sorted(self.m.items(), key = operator.itemgetter(1), reverse = True)
                    self.largest_component = max(networkx.connected_components(self.g), key = len)
                    self.x.append(i * 1. / self.n)
                    self.R += len(self.largest_component) * 1. / self.n
                    self.y.append(len(self.largest_component) * 1. / self.n)

            elif attack_type == 'random':
                self.g = networkx.convert_matrix.from_numpy_array(self.matrix)
                self.l = [(self.node, 0) for self.node in self.g.nodes()]
                random.shuffle(self.l)
                self.x = []
                self.y = []
                self.largest_component = max(networkx.connected_components(self.g), key = len)
                self.n = len(self.g.nodes())
                self.x.append(0)
                self.y.append(len(self.largest_component) * 1. / self.n)
                self.R = 0.0
                for i in range(1, self.n):
                    self.g.remove_node(self.l.pop(0)[0])
                    self.largest_component = max(networkx.connected_components(self.g), key = len)
                    self.x.append(i * 1. / self.n)
                    self.R += len(self.largest_component) * 1. / self.n
                    self.y.append(len(self.largest_component) * 1. / self.n)

            self.all_vulnerability[subj] = np.array(0.5 - self.R / self.n)

        return self.all_vulnerability





    def network_level_vulnerability(self, sbj_number, nodes_number, label_dic, make_symmetric=True, binarize=False,
                                    recalculate = False,  attack_type='target'):
        '''
        Performs network level robustness analysis based on nodal degree


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
        recalculate: Boolean|
                        It will use sequential (recalculate = True) or
                        simultaneous (recalculate = False) approach.
                        Default is False

        Returns
        -------

        vulnerability: np.array |
                    The overall vulnerability for every network

        '''

        with open(self.net_label_txt) as f:
            net=f.read().splitlines()

        self.all_vulnerability = np.zeros([sbj_number, len(net)])

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
            if attack_type == 'target':
                for network in net:
                    self.net2use = self.matrix[label_dic[network]]
                    self.net2use = self.net2use[:,label_dic[network]]
                    self.deg = bct.algorithms.degrees_und(self.net2use)
                    self.g = networkx.convert_matrix.from_numpy_array(self.net2use)
                    self.m = dict(enumerate(self.deg.flatten(), 0))
                    self.l = sorted(self.m.items(), key = operator.itemgetter(1), reverse = True)
                    self.x = []
                    self.y = []
                    self.largest_component = max(networkx.connected_components(self.g), key = len)
                    self.n = len(self.g.nodes())
                    self.x.append(0)
                    self.y.append(len(self.largest_component) * 1. / self.n)
                    self.R = 0.0
                    for i in range(1, self.n - 1):
                        self.g.remove_node(self.l.pop(0)[0])
                        if recalculate:
                            self.net2use = networkx.convert_matrix.to_numpy_array(self.g)
                            self.g = networkx.convert_matrix.from_numpy_array(self.matrix)
                            self.deg = bct.algorithms.degrees_und(self.net2use)
                            self.m = dict(enumerate(self.deg.flatten(), 0))
                            self.l = sorted(self.m.items(), key = operator.itemgetter(1), reverse = True)
                        self.largest_component = max(networkx.connected_components(self.g), key = len)
                        self.x.append(i * 1. / self.n)
                        self.R += len(self.largest_component) * 1. / self.n
                        self.y.append(len(self.largest_component) * 1. / self.n)

                    self.all_vulnerability[subj, net.index(network)] = np.array(0.5 - self.R / self.n)

            elif attack_type == 'random':
                for network in net:
                    self.net2use = self.matrix[label_dic[network]]
                    self.net2use = self.net2use[:,label_dic[network]]
                    self.deg = bct.algorithms.degrees_und(self.net2use)
                    self.g = networkx.convert_matrix.from_numpy_array(self.net2use)
                    self.l = [(self.node, 0) for self.node in self.g.nodes()]
                    random.shuffle(self.l)
                    self.x = []
                    self.y = []
                    self.largest_component = max(networkx.connected_components(self.g), key = len)
                    self.n = len(self.g.nodes())
                    self.x.append(0)
                    self.y.append(len(self.largest_component) * 1. / self.n)
                    self.R = 0.0
                    for i in range(1, self.n):
                        self.g.remove_node(self.l.pop(0)[0])
                        self.largest_component = max(networkx.connected_components(self.g), key = len)
                        self.x.append(i * 1. / self.n)
                        self.R += len(self.largest_component) * 1. / self.n
                        self.y.append(len(self.largest_component) * 1. / self.n)

                    self.all_vulnerability[subj, net.index(network)] = np.array(0.5 - self.R / self.n)

        return self.all_vulnerability
