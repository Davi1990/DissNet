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

    def __init__(self, matrices_files, labels_dic):

        self.matrices_files = matrices_files
        self.labels_dic = labels_dic

    def extract_conn(self, sbj2use, label2use):
        self.matrix = pd.read_csv(self.matrices_files[sbj2use], sep= ' ', header=None)
        self.matrix = np.array(self.matrix)
        np.fill_diagonal(self.matrix,0)
        self.subj_matrix = self.matrix[self.labels_dic[label2use]]
        self.subj_matrix = self.subj_matrix[:,self.labels_dic[label2use]]
        self.streamlines_sum = np.sum(np.sum(self.subj_matrix))
        self.conn_measure = self.streamlines_sum/len(self.labels_dic[label2use])

        return self.conn_measure
