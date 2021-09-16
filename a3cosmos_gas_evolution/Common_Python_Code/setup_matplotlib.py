#!/usr/bin/env python
# 

import os, sys
from matplotlib import pyplot as plt
from matplotlib import ticker as ticker
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from matplotlib import font_manager
import matplotlib as mpl
import numpy as np

def setup_matplotlib():
    
    if os.path.isdir(os.path.expanduser('~')+os.sep+'Library'+os.sep+'Fonts'):
        font_dirs = [os.path.expanduser('~')+os.sep+'Library'+os.sep+'Fonts']
        font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
        font_list = font_manager.createFontList(font_files)
        font_manager.fontManager.ttflist.extend(font_list)
        font_manager.findfont('NGC', rebuild_if_missing=True)
    
    mpl.rcParams['font.family'] = 'NGC'
    
    #mpl.rcParams['text.usetex'] = True
    mpl_version = np.array(mpl.__version__.split('.')).astype(int)
    if mpl_version[0] > 3 or (mpl_version[0] >= 3 and mpl_version[1] >= 1):
        # since matplotlib version 3.1.0, mpl.rcParams['text.latex.preamble'] is a single str not a list
        mpl.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
        mpl.rcParams['text.latex.preamble'] += '\n'
        mpl.rcParams['text.latex.preamble'] += r'\makeatletter \newcommand*{\rom}[1]{\expandafter\@slowromancap\romannumeral #1@} \makeatother'
    else:
        # since matplotlib version 3.1.0, mpl.rcParams['text.latex.preamble'] is a single str not a list
        mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command
        mpl.rcParams['text.latex.preamble'].append(r'\makeatletter \newcommand*{\rom}[1]{\expandafter\@slowromancap\romannumeral #1@} \makeatother')
    
    mpl.rcParams['axes.labelsize'] = '16' # https://matplotlib.org/users/customizing.html
    mpl.rcParams['axes.grid'] = True
    mpl.rcParams['axes.axisbelow'] = True
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['xtick.minor.visible'] = True
    mpl.rcParams['ytick.minor.visible'] = True
    mpl.rcParams['xtick.labelsize'] = '13'
    mpl.rcParams['ytick.labelsize'] = '13'
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True
    #mpl.rcParams['grid.color'] = 'b0b0b0'
    mpl.rcParams['grid.linestyle'] = '--'
    mpl.rcParams['grid.linewidth'] = 0.25
    mpl.rcParams['grid.alpha'] = 0.8
    mpl.rcParams['legend.fontsize'] = '12'
    mpl.rcParams['legend.borderaxespad'] = 0.2 # space between legend border and axis
    #mpl.rcParams['legend.borderpad'] = 0.2 # space between legend content and legend border
    #mpl.rcParams['legend.handletextpad'] = 0.05
    #mpl.rcParams['legend.labelspacing'] = None # in font-size units
    #mpl.rcParams['legend.columnspacing'] = None # in font-size units
    #mpl.rcParams['legend.ncol']


