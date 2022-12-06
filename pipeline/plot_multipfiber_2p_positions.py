# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 17:24:14 2022

@author: qianl
function1 :This function is used to plaot the 3d visualization of multifiber positions
function 2: plot the implatation position of 2p lenses
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict
import os
import random
import matplotlib as mpl
import seaborn as sns

#%% multifiber positions' dict
positions = pd.read_excel('D:/PhD/Photometry/DATA/multifiber_locations.xlsx', 
                          index_col=None, header=0)

color_mapping = {'AMOT':'red',
                'PMOT':'orange',
                'ALOT':'yellow',
                'PLOT':'green',
                'MNacS':'cyan',
                'LNacS':'purple'}
def colorize_region(x):
    return color_mapping[x]

positions['colorized_region'] = positions['region'].apply(colorize_region)

fig = plt.figure(figsize = (5,5))
ax = fig.add_subplot(111, projection='3d')
# sns.scatterplot(data = positions, x = "x", y = "y", z = 'z',hue = 'region', size = 5)
for region in color_mapping.keys():
    
    ax.scatter(positions[positions['region']==region].x, positions[positions['region']==region].y, positions[positions['region']==region].z, 
                     c=color_mapping[region], alpha=1,label = region)
ax.set_xlabel('$x$')
ax.set_xticks([1,1.5])
ax.set_ylabel('$y$')
ax.set_yticks([1.0,1.2,1.4,1.6,1.8])
ax.set_zlabel('$z$')
ax.set_zticks([-6.3,-6.5,-6.7,-6.9,-7.1])
plt.legend(bbox_to_anchor=(1.1, 0.5),frameon = False,loc = 'center left')

ax.view_init(-1, -20)

#%% 2p position plot
#to be implemented