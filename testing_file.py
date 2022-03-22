# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 12:50:45 2022

@author: qianl
"""

import scipy
import numpy as np
import pandas as pd
import matplotlib as plt
from os.path import dirname, join as pjoin
import scipy.io as sio
import math
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict
import os
import random
import matplotlib as mpl
import re
import csv
import pickle

import scipy.stats as stats
from scipy import signal
from statsmodels.gam.api import GLMGam
import statsmodels.gam.smooth_basis as sb
import statsmodels.api as sm
#%%
mu = 0
variance = np.sqrt(3)

sigma = math.sqrt(variance)
x = np.linspace(mu-3*sigma, mu+3*sigma,50)
y = stats.norm.pdf(x,mu,sigma)
plt.plot(x,y)
plt.show()

licks = [0.01,0.05,0.056,0.1,0.3,0.5,2.234,2.4,2.7,3.4,3.8,3.9,5,6,7,8,8.1,8.2,8.3,8.4]
converted_licks = np.histogram(licks,bins = 90,range = (0,9))

convolved = signal.convolve(converted_licks[0],y,mode = 'same')
x = np.linspace(0,9,len(convolved))
plt.plot(x,convolved)
plt.show()
# print(convolved)
#%%

x = np.linspace(0,2,10)
knots = sb.get_knots_bsplines(x,df = 4,)#spacing = 'equal')

basis = sb._eval_bspline_basis(x,degree = 3,knots=knots)[0]

plt.plot(np.convolve(basis[:,0],converted_licks[0],mode = 'full'))
plt.plot(np.convolve(basis[:,1],converted_licks[0],mode = 'full'))
plt.plot(np.convolve(basis[:,2],converted_licks[0],mode = 'full'))