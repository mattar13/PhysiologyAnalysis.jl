#This is going to be a really difficult use of this function, but eventually would like to utilize 
#Make this a standalone 
using PyCall
using Conda

#%% Create a new conda environment
Conda.install("cellpose") #Run this to build the CellPose program

py"""

from cellpost import model
"""