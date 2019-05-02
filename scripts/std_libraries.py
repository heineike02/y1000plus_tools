#Indicate operating environment and import core modules
import os 
location_input = input("what computer are you on? a = Ben's laptop, b = gpucluster, c = Ben's desktop, d = other")
location_dict = {'a': "C:\\Users\\BMH_work\\", 'b': "/home/heineike/",
                 'c': "C:\\Users\\Ben\\Documents\\", 'd':'you need to add your location to the location_dict'}
home_dir = location_dict[location_input]
print("home directory is " + home_dir)
base_dir = home_dir + os.path.normpath('github/y1000plus_tools') + os.sep
print("y1000plus_tools dir is " + base_dir )
y1000plus_dir_options = {'b':home_dir + os.path.normpath("genomes/y1000plus") + os.sep, 
                         'c': home_dir + os.path.normpath('github/expression_broad_data/expression_data/promoter_phylogenies/y1000plus') + os.sep
                        }
y1000plus_dir = y1000plus_dir_options[location_input]
print("y1000plus data dir is " + y1000plus_dir)

import sys

if not(base_dir in sys.path): 
    sys.path.append(base_dir)
    print("Added " + base_dir + " to path" )

print("Importing y1000plus_tools.py")
from core import y1000plus_tools
y1000plus_tools.home_dir = home_dir
y1000plus_tools.base_dir = base_dir
y1000plus_tools.y1000plus_dir = y1000plus_dir


io_library_path = home_dir + os.path.normpath('github/expression_broad_data') + os.sep
io_library_path_core = io_library_path + 'core' + os.sep
if not(io_library_path_core in sys.path):
    sys.path.append(io_library_path_core)
    print("Added " + io_library_path_core + " to path" )

print("importing io_library.py")
import io_library
io_library.base_dir = io_library_path 
io_library.data_processing_dir = io_library_path + os.path.normpath('expression_data') + os.sep

print('sys.path : \n')
print(sys.path)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
import seaborn as sns
## Add to std library
import pickle
import subprocess


from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from ete3 import Tree, SeqMotifFace, TreeStyle, add_face_to_node, RectFace, NodeStyle  
#ete3 is not officially supported on windows, and so must be loaded via pip: 
# pip install -U https://github.com/etetoolkit/ete/archive/qt5.zip
# ref: https://groups.google.com/forum/#!topic/etetoolkit/6NblSBPij4o
#20181031: got this error message: twisted 18.7.0 requires PyHamcrest>=1.9.0, which is not installed.

# In order to view ete3 created trees on the gpucluster, you need to use a virtual X server:
from pyvirtualdisplay import Display
display = Display(visible=False, size=(1024, 768), color_depth=24)
display.start()
