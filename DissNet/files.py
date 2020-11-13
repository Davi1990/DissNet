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

    path: path where connectome are located
    connectome_name: name of connectome .csv files
    label_txt (optional): name of nodes' label .txt file
    network_txt (optional): name of networks' label .txt file
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
