Files form WC computational model ('A Whole-Cell Computational Model Predicts Phenotype from Genotype', https://doi.org/10.1016/j.cell.2012.05.044):

These files were used as input files to retrieve different kind of information about Mycoplasma genitalium and MG simulation data.

- 'protein_data.json' : MG protein database downloaded from https://www.wholecellkb.org/
- 'genes_data.json' = MG gene database downloaded from https://www.wholecellkb.org/
- 'simulation-1.h5': .h5 simulation file of MG wild-type cell simulation set #1-1, downloaded from WholeCellSimDB https://www.wholecellsimdb.org/simulation_batch/1
- 'getSeriesData_monomers1.json' | 'getSeriesData_complex1.json': json file with positions of DNA binding monomers/complexes calculated in MG wild-type cell simulation set #1-1. Downloaded from WholeCellViz http://www.wholecellviz.org/getSeriesData.php?sim_id=2011_10_19_02_53_45_1&class_name=Chromosome&attr_name=monomerBoundSites | http://www.wholecellviz.org/getSeriesData.php?sim_id=2011_10_19_02_53_45_1&class_name=Chromosome&attr_name=complexBoundSites

Cell size extracted from simulation-1.h5 in states/Geometry/volume/data

RNA polymerase positions on the genome and RNA polymerase chromosome exploration extracted from simulation-1.h5 in states/RNAPolymerase/

Ribosome and RNAs count and positions extracted from simulation-1.h5 in states/Ribosome/. 
The original ribosome position reported in the WC simulation was edited when: 1) two ribosomes were reported to be the same position onto the same mRNA 2) the ribosome was located less than 60 nucleotides from the end or the beginning of a transcript.

Note that frame number is +1 when using json files (i.e. frame 150 in simulation-1.h5 corresponds to frame 150 in getSeriesData_monomers1.json)

Scripts available upon request mmaritan@scripps.edu
