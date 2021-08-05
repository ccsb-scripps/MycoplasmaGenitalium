Script : Python's script integrating data from the WC computational model and structural information for MG proteins (in input_files folder). It generates: csv table used to make recipe in mesoscope, csv table used to make lattice-based models of supercoiled mycoplasma nucleoids with NucleoidLattice by Davaid S. Goodsell

WC-MG-CellPACK-RecipeBuilder-short.py: create csv table with a draft recipe readable by Mesoscope
WC-MG-CellPACK-input_LatticeNucucleoid.py: create csv table with data used by LatticeNucleoids to build MG nucleoid models
WC-MG-CellPACK-functions.py: custom functions to run the scripts above

******
To run WC-MG-CellPACK-RecipeBuilder-short.py you need:
- WC-MG-CellPACK-functions
- Input folder with all the files 
- Create your output folder 

Option 1: change the frame number will change the copy number and retrieve the data from the WC simulation. Retrieving ingredients copy number typically takes few minutes
Option 2: generate the 'curated' recipe or 'auto' recipe by changing 'use_method' at the beginning of the script

******
To run WC-MG-CellPACK-input_LatticeNucucleoid.py you need:
- WC-MG-CellPACK-functions
- Input folder with: protein_data.json, genes_data.json, simulation-1.h5, getSeriesData_monomers1.json, getSeriesData_complex1.json, S3K-Transcription units.csv
- Output folder 

Cell size extracted from simulation-1.h5 in states/Geometry/volume/data

RNA polymerase positions on the genome and RNA polymerase chromosome exploration extracted from simulation-1.h5 in states/RNAPolymerase/

Ribosome and RNAs count and positions extracted from simulation-1.h5 in states/Ribosome/. 
The original ribosome position reported in the WC simulation was edited when: 1) two ribosomes were reported to be the same position onto the same mRNA 2) the ribosome was located less than 60 nucleotides from the end or the beginning of a transcript.

Note that frame number is +1 when using json files (i.e. frame 150 in simulation-1.h5 corresponds to frame 150 in getSeriesData_monomers1.json)

******
Info on Input files:
- protein_data.json : MG protein database downloaded from https://www.wholecellkb.org/
- genes_data.json : MG gene database downloaded from https://www.wholecellkb.org/
- simulation-1.h5: .h5 simulation file of MG wild-type cell simulation set #1-1, downloaded from WholeCellSimDB https://www.wholecellsimdb.org/simulation_batch/1
- getSeriesData_monomers1.json | getSeriesData_complex1.json: json file with positions of DNA binding monomers/complexes calculated in MG wild-type cell simulation set #1-1. Downloaded from WholeCellViz https://www.wholecellviz.org/getSeriesData.php?sim_id=2011_10_19_02_53_45_1&class_name=Chromosome&attr_name=monomerBoundSites | https://www.wholecellviz.org/getSeriesData.php?sim_id=2011_10_19_02_53_45_1&class_name=Chromosome&attr_name=complexBoundSites
- all_dict.json: dictionary holding the list of structural models available for each ingredients. Contains: model name, chain selection, BU selection, membrane orientation (if available), template (only for homology models)
- compartment_updated.json: dictionary with protein localization info updated from the WC computational model
- HHpred_scores.csv: list of Hhscores calculated for each protein ingredients
- ingredient_function.csv: list of ingredient functions. Total of 31 functional categories identified in CYT-MG model (doi: 10.1016/j.jmgm.2015.02.004)  
- method_selection.csv: indicated which structure was chosen for the curated recipe
- WholeCellData.xlsx: from https://www.sciencedirect.com/science/article/pii/S0092867412007763?via%3Dihub#app2, tabs 'S3M-Protein Monomers' added column J, 'S3N-Macromolecular complexes' added column E
- S3K-Transcription units.csv: coming from supplementary materials of WC computational model paper
- uniprot.csv: list of Uniprot Ids for MG genes
updated_function_dictionary.json: list of protein functions from WC-computational model integrated with CYT-MG model (doi: 10.1016/j.jmgm.2015.02.004)