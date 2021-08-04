# MycoplasmaGenitalium
Modelling of Mycoplasma genitalium
Repository covering the ccsb work on the mycoplasma genitalium modelling project based on the [WC computational model](https://doi.org/10.1016/j.cell.2012.05.044).


![mesoscope](https://github.com/ccsb-scripps/MycoplasmaGenitalium/blob/main/Models/CaptureMesoscope.PNG)



The repository consist in the following folders :
* Data : the copy of the data folder that is used by cellPACKgpu to load and generate the final 3D models.
    * lattices
    * palettes
    * proteins
    * recipes
* LatticeNucleoids : David Goodsell software to generate initial configuration of DNA and RNA for the model.
    * fortran source
    * input files
    * README
* Models : Exampe of generated models with cellPACKgpu in different file formats. cif and bcif can be directly loaded in [molstar](https://molstar.org/)
    * Model at Beads resolution in pdb format (.pdb) 
    * Model at Atomic resolution in cif (.zip) and bcif format ( convertion done with molstar )
    * Model in cellpack binary format (.bin)
    * Lipids Membrane model in cif and bcif format
    * Bundle in zip format for mesoscope
    * Example of molx file   
* Script : Python's script integrating all the data from the WC computational model (in input_files folder) which generate csv table used to make recipe in mesoscope
* WholeCellData :  Files form WC computational model
* cellPACKgpu binary release : standalone windows executable 
## 
