#%% use this to plot any 
using PhysiologyAnalysis, PyPlot
path = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Data\ERG\JGP_Files\Traces\Cones_photopic\Adult\NR\2019_24_09_WT_P34_m2_green_photopic\nd1_100p_1ms\Average076.abf"
data = readABF(path) 
data_filter!(data, t_pre = 0.5, t_post = 0.5)
plot_experiment(data)