# MycoplasmaGenitalium
Modelling of Mycoplasma genitalium

[Maritan M, Autin L, Karr J, Covert MW, Olson AJ, Goodsell DS. Building Structural Models of a Whole Mycoplasma Cell. J Mol Biol. 2022 Jan 30;434(2):167351. doi: 10.1016/j.jmb.2021.167351. Epub 2021 Nov 10. PMID: 34774566; PMCID: PMC8752489.](https://pubmed.ncbi.nlm.nih.gov/34774566/)


Repository covering the ccsb work on the mycoplasma genitalium modelling project based on the [WC computational model](https://doi.org/10.1016/j.cell.2012.05.044). 

![mesoscope](https://github.com/ccsb-scripps/MycoplasmaGenitalium/blob/main/Models/CaptureMesoscope.PNG)
Recipes and models are available on [Mesoscope](https://mesoscope.scripps.edu/beta/) in the list of available example (Load->New Recipe->From Example->Mycoplasma Genitalium


The repository consist in the following folders :
* cellPACK_Data : the copy of the data folder that is used by cellPACKgpu to load and generate the final 3D models.
    * lattices
    * palettes
    * proteins
    * recipes
    * README
* LatticeNucleoids : David Goodsell software to generate initial configuration of DNA and RNA for the model.
    * fortran source
    * input files
    * README
* Models : Exampe of generated models with cellPACKgpu in different file formats. cif and bcif can be directly loaded in [molstar](https://molstar.org/)
    * Model at Beads resolution in pdb format (.pdb) 
    * Model at Atomic resolution in cif (.zip) and bcif format ( convertion done with molstar ). See one example using a [bcif file in molstar](https://molstar.org/viewer/?structure-url=https://raw.githubusercontent.com/ccsb-scripps/MycoplasmaGenitalium/main/Models/cellpack_atom_instances_149_curated.bcif&structure-url-format=mmcif&structure-url-is-binary=1). Switch the representation to spacefill and add a clipping plane. 
    * Model in cellpack binary format (.bin) can be loaded in cellPACKgpu, mesoscope and molstar.
    * Lipids Membrane model in cif and bcif format, with and without holes. Theses models presents overlap between lipids.
    * Bundle in zip format for mesoscope
    * Example of molx filem can be directly drag and drop onto the molstar viewer. You can open the following [premade molstar scene with membrane](https://molstar.org/viewer/?snapshot-url=https://raw.githubusercontent.com/ccsb-scripps/MycoplasmaGenitalium/main/Models/mol-star_state_1189.molx&snapshot-url-type=molx) and [without](https://molstar.org/viewer/?snapshot-url=https://raw.githubusercontent.com/ccsb-scripps/MycoplasmaGenitalium/main/Models/mol-star_state_1189_no_membrane.molx&snapshot-url-type=molx). If you have a powerfull workstation try to load the [model with spacefill representation](https://molstar.org/viewer/?snapshot-url=https://raw.githubusercontent.com/ccsb-scripps/MycoplasmaGenitalium/main/Models/mol-star_state_1189_spacefill.molx&snapshot-url-type=molx).
* Scripts : Python's script integrating data from the WC computational model and structural information for MG proteins (in input_files folder). It generates: csv table used to make recipe in mesoscope, csv table used to make lattice-based models of supercoiled mycoplasma nucleoids with LatticeNucleoids by David S. Goodsell
* cellPACKgpu binary release : standalone windows executable 
## 
