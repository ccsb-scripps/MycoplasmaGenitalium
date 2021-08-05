# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 01:06:01 2021
@author: MM
Get input data from WC-MG computational simulation for LatticeNucleoid
"""
exec(open('WC-MG-CellPACK-functions.py').read())
#input_dir ='C:\\Users\\marti\\Documents\\wholecellworkflow\\FOR GIT-cellPACKpaper\\input_files\\'
#input_dir ='input_files'+os.sep
input_dir = 'C:/Users/marti/Documents/wholecellworkflow/FOR GIT-cellPACKpaper/input_files/'
output_dir= 'C:/Users/marti/Documents/WC-model/scripts_output_junk/'
#output_folder = 'scripts_output'+os.sep

# Data for LATTICE NUCLEOID

#gets data from getSeriesData_monomers1.json | getSeriesData_complex1.json
#define the frames for monomers and complexes 
frames_m =[150, 1185, 6974] #nb this frame must exist in the json file fname_m (for json the frame is +1)
frames_c=[146, 1190, 6961] #nb this frame must exist in the json file fname_c (for json the frame is +1)

#choromosome binding proteins at a specific frame
for frame in frames_m:
    chromosome_proteins(fname_m, frame, fname_c, frames_c[frames_m.index(frame)],sim)  

#rna polimerase, ribosome and rna data
#gets data from hdf5 file     
sim_c_frames=[145,1189,6960] 
for i in sim_c_frames:
    #RNA polymerase position, lenght of transcript
    rna_poly(sim, i)
    #all mRNAs + positions for ribosomes 
    Ribosome_RNAs(sim,i)