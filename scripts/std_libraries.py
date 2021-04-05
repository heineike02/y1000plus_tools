### means I need to add it to my python environment

#Indicate operating environment and import core modules
import os 
location_input = input("what computer are you on? a = Ben's laptop, b = gpucluster, c = Ben's desktop, d = other")
location_dict = {'a': "C:\\Users\\BMH_work\\", 'b': "/home/heineike/",
                 'c': "C:\\Users\\Ben\\Documents\\", 'd':'you need to add your location to the location_dict'}
figsave_dict = {'a': "C:\\Users\\BMH_work\\Google Drive\\UCSF\\ElSamad_Lab\\PKA\\Manuscript\\" , 
                'b': "/home/heineike/scratch/",
                'c': "C:\\Users\\Ben\\Google Drive\\UCSF\\ElSamad_Lab\\PKA\\Manuscript\\", 
                'd': 'you need to add your location to the figsave dict'}
figsave_dir = figsave_dict[location_input]
home_dir = location_dict[location_input]
print("home directory is " + home_dir)
base_dir = home_dir + os.path.normpath('github/y1000plus_tools') + os.sep
print("y1000plus_tools dir is " + base_dir )
y1000plus_dir_options = {'a': base_dir + os.path.normpath('genomes/y1000plus') + os.sep,
                         'b':home_dir + os.path.normpath("genomes/y1000plus") + os.sep, 
                         'c': home_dir + os.path.normpath('github/yeast_esr_expression_analysis/expression_data/promoter_phylogenies/y1000plus') + os.sep
                        }
y1000plus_dir = y1000plus_dir_options[location_input]
print("y1000plus data dir is " + y1000plus_dir)

import sys

if not(base_dir in sys.path): 
    sys.path.append(base_dir)
    print("Added " + base_dir + " to path" )

yeast_esr_exp_path = home_dir + os.path.normpath('github/yeast_esr_expression_analysis') + os.sep
#io_library_path_core = io_library_path + 'core' + os.sep
if not(yeast_esr_exp_path in sys.path):
    sys.path.append(yeast_esr_exp_path)
    print("Added " + yeast_esr_exp_path + " to path" )
    
print("Importing y1000plus_tools.py")
import y1000plus_tools

y1000plus_tools.home_dir = home_dir
y1000plus_tools.base_dir = base_dir
y1000plus_tools.y1000plus_dir = y1000plus_dir


y1000plus_tools.yeast_esr_exp.base_dir = yeast_esr_exp_path
y1000plus_tools.yeast_esr_exp.data_processing_dir = yeast_esr_exp_path + os.path.normpath('expression_data') + os.sep


# Since y1000plus_tools loads io_library, just use those functions
# if not(io_library_path_core in sys.path):
#     sys.path.append(io_library_path_core)
#     print("Added " + io_library_path_core + " to path" )


print("importing yeast_esr_exp")
print(sys.path)
import yeast_esr_exp
yeast_esr_exp.base_dir = yeast_esr_exp_path
yeast_esr_exp.data_processing_dir = yeast_esr_exp_path + os.path.normpath('expression_data') + os.sep


print('sys.path : \n')
print(sys.path)

import copy
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
from matplotlib.gridspec import GridSpec

import seaborn as sns
## Add to std library
import pickle
import subprocess

from collections import Counter, OrderedDict
from itertools import chain

import scipy.spatial.distance as spd
#import statsmodels.graphics.gofplots as stats_graph
import scipy.cluster.hierarchy as sch
from statsmodels.distributions.empirical_distribution import ECDF

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, IUPAC
from Bio import SeqIO
from Bio import pairwise2
from Bio import motifs
from Bio import AlignIO
from Bio import Align

import gffutils
from ete3 import Tree, SeqMotifFace, TreeStyle, add_face_to_node, RectFace, NodeStyle, TextFace, AttrFace
#ete3 is not officially supported on windows, and so must be loaded via pip: 
# pip install -U https://github.com/etetoolkit/ete/archive/qt5.zip
# ref: https://groups.google.com/forum/#!topic/etetoolkit/6NblSBPij4o
#20181031: got this error message: twisted 18.7.0 requires PyHamcrest>=1.9.0, which is not installed.

# In order to view ete3 created trees on the gpucluster, you need to use a virtual X server:
### from pyvirtualdisplay import Display
### display = Display(visible=False, size=(1024, 768), color_depth=24)
### display.start()

#for scraping internet data (e.g. ncbi, YGOB)
import requests
from bs4 import BeautifulSoup
#from lxml import etree    #parses xml output
