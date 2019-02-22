#Indicate operating environment and import core modules
import os 
location_input = input("what computer are you on? a = Ben's laptop, b = gpucluster, c = Ben's desktop, d = other")
location_dict = {'a': "C:\\Users\\BMH_work\\", 'b': "/home/heineike/",
                 'c': "C:\\Users\\Ben\\Documents\\", 'd':'you need to add your location to the location_dict'}
home_dir = location_dict[location_input]
print("home directory is " + home_dir)
base_dir = home_dir + os.path.normpath('github/y1000plus_tools') + os.sep
print("y1000plus_tools dir is " + base_dir )
y1000plus_dir = home_dir + os.path.normpath("genomes/y1000plus") + os.sep
print("y1000plus genomes dir is " + y1000plus_dir)

import sys

if sys.path[-1] != base_dir:
    sys.path.append(base_dir)
    print("Added " + base_dir + " to path: " )
    print(sys.path)

print("Importing y1000plus_tools.py")
from core import y1000plus_tools
%load_ext autoreload
%autoreload 2


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from ete3 import Tree
#ete3 is not officially supported on windows, and so must be loaded via pip: 
# pip install -U https://github.com/etetoolkit/ete/archive/qt5.zip
# ref: https://groups.google.com/forum/#!topic/etetoolkit/6NblSBPij4o
#20181031: got this error message: twisted 18.7.0 requires PyHamcrest>=1.9.0, which is not installed.
