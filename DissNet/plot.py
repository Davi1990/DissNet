import matplotlib.pyplot as plt
import os
import sys
import mayavi
from surfer import Brain
from mayavi import mlab
import numpy as np
import pandas as pd
from math import pi
import networkx as nx
from tvtk.api import tvtk


def plot_surface(subjects_dir, atlas, subject, hemi, surf, opacity=None ,scale_factor=None, threshold=None,
                 connectivity_matrix=None, connectivity_thr=None):
    '''
    plot each nodes using pysurfer

    Parameters
    ----------

    subjects_dir = path to the subjects' folders
    atlas: excel file |
            please se example available in the repo (e.g. new_atlas_coords.xlsx)
    subject = subject folder name that you want to plot
    hemi = str |
        (ie 'lh', 'rh', 'both', or 'split'). In the case
        of 'both', both hemispheres are shown in the same window.
        In the case of 'split' hemispheres are displayed side-by-side
        in different viewing panes
    surf:  str |
        freesurfer surface mesh name (ie 'white', 'inflated', etc.). A cerebellum
        surface model is also available in the subject folder of the repo. For
        plot it simply enter 'cerebellum'
    opacity (optional)= float |
        sets the overall opacity of colors, maintains transparent regions.
        Default is 0.2
    scale_factor (optional) = np.array |
        vector (dim number of regions x 1) representing the size of each node.
        For example, it can be connectivity of that specific node calculated
        using any functions available in metrics. Default is 0.5
    threshold= float |
        threshold to plot only the nodes that overcome that specific values
        Default is None
    connectivity_matrix (optional) = np.array |
        vector (dim number of regions x number of regions) representing
        the connectivity of each node.
    connectivity_thr (optional) = scalar |
    
    '''
    np.set_printoptions(suppress=True) #prevent numpy exponential
    label = pd.read_excel(atlas, header=None)
    matrix = np.array(label)
    coords = np.array(matrix[:,5:8], dtype='float32')
    color = np.array(matrix[:,2:5], dtype='float32')

    def NormalizeData(data):
         return (data - np.min(data)) / (np.max(data) - np.min(data))
    if scale_factor is not None:
        scale_factor_norm = NormalizeData(scale_factor)
    else:
        scale_factor_norm= np.zeros(matrix.shape[0]) + 0.5



    if threshold==None:
     scale_factor_norm= scale_factor_norm
    else:
        scale_factor_norm[scale_factor_norm < threshold] = 0

    if opacity==None:
     opacity= 0.2
    else:
        opacity = opacity

    if surf=='cerebellum':
        hemi= 'rh'
    else:
        hemi=hemi

    brain = Brain(subject, subjects_dir=subjects_dir, hemi=hemi, surf=surf, cortex='ivory', alpha=opacity)
    for roi in range(coords.shape[0]):
       brain.add_foci(coords[roi], color=(color[roi][0]/255.0, color[roi][1]/255.0, color[roi][2]/255.0), scale_factor=scale_factor_norm[roi], hemi='rh')



    def edgelist(df):
       a = df.values
       c = df.columns
       n = len(c)

       c_ar = np.array(c)
       out = np.empty((n, n, 2), dtype=c_ar.dtype)

       out[...,0] = c_ar[:,None]
       out[...,1] = c_ar

       mask = ~np.eye(n,dtype=bool)
       df_out = pd.DataFrame(out[mask], columns=[['Source','Target']])
       df_out['Weight'] = a[mask]
       return df_out

    if connectivity_matrix is not None:
        connectivity_matrix = pd.DataFrame(connectivity_matrix)
        edge_list = edgelist(connectivity_matrix)
        values = np.array(edge_list)
        pos = np.where(values[:,2]>connectivity_thr)
        edges = values[:,0:2][pos]
        weights =  NormalizeData(values[:,2][pos])
        for xx in range(edges.shape[0]):
            mlab.plot3d(coords[edges[xx]][:,0],coords[edges[xx]][:,1], coords[edges[xx]][:,2],
                        tube_radius=weights[xx]+0.25)

    mlab.show()


def spider_plot(labels, values, colour=None, linestyle=None, linecolour=None, linewidth=None, label_colour=None, label_size=None, alpha=None, ylim=None):

    if ylim==None:
        ylim= np.max(values) + (10 - np.max(values)  % 10)
    else:
        ylim=ylim

    if alpha==None:
        alpha = 0.5
    else:
        alpha= alpha

    if label_colour==None:
        label_colour='black'
    else:
        label_colour=label_colour

    if label_size==None:
        label_size=8
    else:
        label_size=label_size

    if linewidth==None:
        linewidth=1
    else:
        linewidth=linewidth

    if linecolour==None:
        linecolour=color
    else:
        linecolour=linecolour

    if linestyle==None:
        linestyle='solid'
    else:
        linestyle=linestyle

    if colour==None:
        colour= 'orange'
    else:
        colour=colour

    values = values * 100 / np.sum(values)
    values = values.tolist()
    values += values[:1]
    N = len(labels)
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]
    ax = plt.subplot(111, polar=True)

    plt.xticks(angles, labels, color=label_colour, size=label_size)

    ax.set_rlabel_position(0)
    plt.ylim(0,ylim)

    ax.plot(angles, values, linecolour , linewidth=linewidth, linestyle=linestyle)

    ax.fill(angles, values, colour, alpha=alpha)
    plt.show()



def bar_plot(network2use, values, align=None, alpha=None, xlabel=None):
    if align==None:
        align='center'
    else:
        align=align

    if alpha==None:
        alpha=0.8
    else:
        alpha=alpha

    if xlabel==None:
        xlabel='Connectivity'
    else:
        xlabel=xlabel

    network2use = pd.read_excel(network2use, header=None)
    Net_label=list(network2use[0])
    barlist= plt.barh(np.arange(len(Net_label)), values, align=align, alpha=alpha)
    plt.yticks(np.arange(len(Net_label)), Net_label)
    plt.xlabel(xlabel)
    plt.xlim(np.min(values) - np.std(values), np.max(values) + np.std(values))
    for net in range(len(Net_label)):
        barlist[net].set_color([network2use[1][net]/255.0,network2use[2][net]/255.0,network2use[3][net]/255.0])

    plt.show()
