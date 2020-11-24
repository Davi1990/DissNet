import DissNet.plot as plotting
import DissNet.metrics as metrics
import DissNet.files as files
import numpy as np

plotting.plot_surface('.', 'examples/new_atlas_coords.xlsx',
                     'fsaverage', 'both', 'cerebellum')

plotting.spider_plot('examples/network_colour.xlsx',
                     np.random.rand(7))

plotting.bar_plot('examples/network_colour.xlsx',
                  np.random.rand(7), percentage_value=True) 
