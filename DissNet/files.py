"""
miscelallaneous functions and classes to work to get .csv structual connectome and labels
Author: Davide Momi, PhD [momi.davide89@gmail.com], https://twitter.com/davemomi
"""

import numpy as np
import pandas as pd
from math import pi
import glob
import seaborn as sns
import matplotlib.pyplot as plt


class Files(object):
    '''
    reading the connectome and nodes' labels files

    Parameters
    ----------

    path: str |
    path where the connectomes are located
    connectome_name: str |
    name of connectome .csv files (e.g. 'new_atlas_Yeo.csv')
    label_txt (optional): str |
    name of nodes' label .txt file (see example 'labels.txt')
    network_txt (optional): str |
    name of networks' label .txt file (see example 'net_labels.txt')

    '''
    def __init__(self, path, connectome_name, label_txt=None, network_txt=None):

        self.path = path
        self.connectome_name = connectome_name
        self.label_txt = label_txt
        self.network_txt = network_txt

    def connectome_files(self):
        self.sc_matrix=glob.glob(self.path + '/*/' + '*' + self.connectome_name)
        return self.sc_matrix


    def labels(self):
        '''
        Returns
        -------

        kwargs : dict |
        Dictionary with keyword arguments to be used for labels
        node and network connectivity extraction.
        For nodes name, key is ['nodes']. For the network index
        of each nodes the key is ['network']

        '''

        with open(self.path + '/' + self.label_txt) as f:
             self.labels=f.read().splitlines()
             with open(self.path + '/' + self.network_txt) as f:
                 self.networks=f.read().splitlines()
                 self.labels_dic = {
                 "nodes": self.labels,
                 "networks": self.networks
                 }
                 for ii in range(len(self.networks)):
                    self.labels_dic.update({self.labels_dic['networks'][ii]:np.array([i for i, s in enumerate(self.labels_dic['nodes']) if self.labels_dic['networks'][ii] in s])})

                 return self.labels_dic
