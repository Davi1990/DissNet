import matplotlib.pyplot as plt
import os
import sys
import mayavi
from surfer import Brain
from mayavi import mlab
import numpy as np
import pandas as pd



def plot_surface(subjects_dir, atlas, subject, hemi, surf, opacity=None ,scale_factor=None):
    '''
    plot each nodes using pysurfer

    Parameters
    ----------

    subjects_dir = path to the subjects' folders
    atlas: excel file
    please se example available in the repo (e.g. new_atlas_coords.xlsx)
    subject = subject folder name that you want to plot
    hemi = str
    (ie 'lh', 'rh', 'both', or 'split'). In the case
    of 'both', both hemispheres are shown in the same window.
    In the case of 'split' hemispheres are displayed side-by-side
    in different viewing panes
    surf:  str
    freesurfer surface mesh name (ie 'white', 'inflated', etc.). A cerebellum
    surface model is also available in the subject folder of the repo. For
    plot it simply enter 'cerebellum'
    opacity (optional)= float
    sets the overall opacity of colors, maintains transparent regions.
    Default is 0.2
    scale_factor (optional) = np.array
    vector (dim number of regions x 1) representing the size of each node.
    For example, it can be connectivity of that specific node calculated
    using any functions available in metrics. Default is 0.5
    '''
    np.set_printoptions(suppress=True) #prevent numpy exponential
    label = pd.read_excel(atlas, header=None)
    matrix = np.array(label)
    coords = np.array(matrix[:,5:8], dtype='float32')
    color = np.array(matrix[:,2:5], dtype='float32')
    def NormalizeData(data):
         return (data - np.min(data)) / (np.max(data) - np.min(data))
    if scale_factor==None:
     scale_factor_norm= np.zeros(matrix.shape[0]) + 0.5
    else:
        scale_factor_norm = NormalizeData(scale_factor) + 0.5

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


def spider_plot(labels, values, colour, linestyle, linecolour, linewidth, label_colour, label_size, ylim):
    values = values * 100 / np.sum(values)
    values = values.tolist()
    values += values[:1]
    N = len(labels)
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]
    ax = plt.subplot(111, polar=True)

    plt.xticks(angles, Net, color=label_colour, size=label_size)

    ax.set_rlabel_position(0)
    plt.ylim(0,ylim)

    ax.plot(angles, values, linecolour , linewidth=linewidth, linestyle=linestyle)

    ax.fill(angles, values, colour, alpha=0.3)
    plt.show()



def bar_plot(Net_label, values, align, alpha, xlabel):
    barlist= plt.barh(np.arange(len(Net_label)), values, align=align, alpha=alpha)
    plt.yticks(np.arange(len(Net_label)), Net_label)
    plt.xlabel(xlabel)
    barlist[0].set_color('purple')
    barlist[1].set_color('skyblue')
    barlist[2].set_color('green')
    barlist[3].set_color('violet')
    barlist[4].set_color('wheat')
    barlist[5].set_color('orange')
    barlist[6].set_color('#cb181d')
    plt.show()
