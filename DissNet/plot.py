import matplotlib.pyplot as plt
import os
import sys
import mayavi
from surfer import Brain
from mayavi import mlab
import numpy as np
import pandas as pd
from math import pi




def plot_surface(subjects_dir, atlas, subject, hemi, surf, opacity=None ,scale_factor=None, threshold=None):
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
    '''
    np.set_printoptions(suppress=True) #prevent numpy exponential
    label = pd.read_excel(atlas, header=None)
    matrix = np.array(label)
    coords = np.array(matrix[:,5:8], dtype='float32')
    color = np.array(matrix[:,2:5], dtype='float32')
    def NormalizeData(data):
         return (data - np.min(data)) / (np.max(data) - np.min(data))
    if scale_factor.any()==None:
     scale_factor_norm= np.zeros(matrix.shape[0]) + 0.5
    else:
        scale_factor_norm = NormalizeData(scale_factor)

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
    #for hemi in ["lh", "rh"]:
    #    brain = Brain(subject, subjects_dir=subjects_dir, hemi=hemi, surf='pial', cortex='ivory', alpha=0.5)
    #    brain.add_foci(coords, color="green", scale_factor=0.8)
    for roi in range(coords.shape[0]):
        brain.add_foci(coords[roi], color=(color[roi][0]/255.0, color[roi][1]/255.0, color[roi][2]/255.0), scale_factor=scale_factor_norm[roi], hemi='rh')

    #brain.add_foci(coords, color="green", scale_factor=0.8)

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
    for net in range(len(Net_label)):
        barlist[net].set_color([network2use[1][net]/255.0,network2use[2][net]/255.0,network2use[3][net]/255.0])

    plt.show()
