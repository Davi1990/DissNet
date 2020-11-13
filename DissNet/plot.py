import matplotlib.pyplot as plt
import os
import sys
import mayavi
from surfer import Brain
from mayavi import mlab
import numpy as np
import pandas as pd


#os.environ['SUBJECTS_DIR'] = sys.argv[0]

#subjects_dir = sys.argv[0]
#subject = sys.argv[1]


#np.set_printoptions(suppress=True) #prevent numpy exponential
#label = pd.read_excel('/home/davide/Desktop/DAVE/ATLASES/Schaefer/Schafer_for_BrainNet.xlsx', header=None)
#matrix = np.array(label)
#coords = matrix[:,0:3]

#brain = Brain(subject, subjects_dir=subjects_dir, hemi='rh', surf='pial', cortex='ivory', alpha=0.5)


#for hemi in ["lh", "rh"]:
#    brain = Brain(subject, subjects_dir=subjects_dir, hemi=hemi, surf='pial', cortex='ivory', alpha=0.5)
#    brain.add_foci(coords, color="green", scale_factor=0.8)

#brain.add_foci(coords, color="green", scale_factor=0.8)

#mlab.show()
def plot_atlas(subjects_dir, atlas, subject, hemi):
    np.set_printoptions(suppress=True) #prevent numpy exponential
    label = pd.read_excel(atlas, header=None)
    matrix = np.array(label)
    coords = matrix[:,0:3]

    brain = Brain(subject, subjects_dir=subjects_dir, hemi=hemi, surf='pial', cortex='ivory', alpha=0.5)


    #for hemi in ["lh", "rh"]:
    #    brain = Brain(subject, subjects_dir=subjects_dir, hemi=hemi, surf='pial', cortex='ivory', alpha=0.5)
    #    brain.add_foci(coords, color="green", scale_factor=0.8)

    brain.add_foci(coords, color="green", scale_factor=0.8)

    mlab.show()
