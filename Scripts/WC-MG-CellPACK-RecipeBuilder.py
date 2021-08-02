# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 16:08:03 2021

RECIPE BUILDER - ASSEMBLES A MG RECIPE COMPATIBLE WITH MESOSCOPE


"""
import os
#1. MAKE d3 DICTIONARY
#    1.1 chian -  quality - membrane - notes
#        1.1.1 draft chains - membrane orientation - uniprot -hhpred
#        1.1.2 quality monomers: modfold 
#        1.1.2 quality complexes: voromqa
#        1.1.3 membrane
#        1.1.4 MM-notes
#        1.1.5 essentiality
#        1.1.6 DNA binding
#        1.1.7 functional_category 
#        1.1.8 functional_category_clustered
#    1.2 homology models 
#        1.2.1 monomers
#            1.2.1.1 template monomers 
#            1.2.1.2 mapping monomers
#            1.2.1.3 all_mapping 
#        1.2.2 complexes
#            1.2.2.1 mapping complexes
#            1.2.2.2  all_mapping_complex_homology
#        1.2.3. MONOMERS + COMPLEXES : homology_models
#    1.3 homologs
#        1.3.1 PDB-homologs 
#        1.3.2 PDB-homologs edited
#        1.3.3 homologs
#    1.4 solved structures
#        1.4.1 mapping solved
#        1.4.2 solved
#    1.5 all_dict : homology models + homologs + solved   
#    1.6 d1 & d2
#        1.6.1 d1 (monomers)
#        1.6.2 d2 (complexes)
#    1.7 d3: d1 + d2
#        1.7.1 first draft 
#        1.7.2 final d3
#2. RNA (tRNA+rRNA_sRNA)
#
#3. COPY NUMBER 
#
#4. recipe as csv file
exec(open('WC-MG-CellPACK-functions.py').read())

input_dir ='input_files'+os.sep
output_dir = output_folder = 'scripts_output'+os.sep

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@           
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@        1: d3      @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@            
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# ####################################################################################
            
                            # 1.1 : CHAINS - QUALITY - MEMBRANE POSITION 
            
#-###################################################################################
            
#**********************************  
# 1.1.1-3: # CHIAN, OFFSET, AXIS DICTIONARY 
#**********************************
#UNIFIED CHIAN SELECTION AND MEMBRANE POSITION DICTIONARY
e = input_dir+'MG_chains_membrane.csv'                       
#lists with pdbs, ingredient names ect. 
pdbs=[]
ingredients=[]
selection =[]
bu=[]
pcpalAxis=[]
offset=[]
with open(e) as csvfile2:
    spamreader2 = csv.reader(csvfile2)
    headers = next(spamreader2)
    for row in spamreader2:
        #list pdbs
        pdbs.append(row[1])
        #list aname
        ingredients.append(row[0])
        #list chains
        selection.append(row[2])
        #list bu
        bu.append(row[3])    
        #list axis
        pcpalAxis.append(row[5])
        #list offset
        offset.append(row[6])
#made an array with the list of ingredients to be able to search later         
cose= np.array(ingredients)
       
#dictionary with chains and membrane positions
dcm={}
for element in ingredients:
    #if there are more structures for an ingredient it will be index>1
    index = np.where(cose == element)[0]
    dcm[element]=[]
    #for those structures that in this dictionary have only one pdb associated
    if len(index)==1:
        #print (element, index[0], ingredients.index(element))#continue
        dcm[element].append({'pdb_model':pdbs[index[0]], 'chain':selection[index[0]], 'bu':bu[index[0]], 'pcpalAxis':pcpalAxis[index[0]], 'offset':offset[index[0]]})
    else:#more structures for an ingredient
        for i in index:
            dcm[element].append({'pdb_model':pdbs[i], 'chain':selection[i], 'bu':bu[i], 'pcpalAxis':pcpalAxis[i], 'offset':offset[i]})

print ('dcm >dictionary chains and membrane' )

#uniprot 
e = input_dir+'uniprot.csv'                         
uniprot={}
with open(e) as csvfile2:
    spamreader2 = csv.reader(csvfile2)
    headers = next(spamreader2)
    for row in spamreader2:
        aname = row[0]
        mw = row[1]
        uniprot_code =row[2]
        uniprot[aname] = {'uniprot':uniprot_code}
print ('uniprot > uniprot monomers')

#DICTIONARY hhPRED all monomers
c =input_dir+'HHpred_scores.csv'   
hh = {}
with open(c) as csvfile:
    spamreader = csv.reader(csvfile)
    headers = next(spamreader)
    for row in spamreader:
        aname = row[0]
        hh_coverage = float(row[8])
        HH=float(row[9])
        hh[aname]={'HHpred':HH, 'HHcov':hh_coverage}
print ('hh > hhpred monomers')

#**********************************  
# 1.1.2: quality monomers: modfold 
# 1.1.2: quality complexes: voromqa 
#**********************************
# file with the scores of all the monomers as calculated by modfold: this file has the mod fold output for all the best models from each method
csvf = input_dir+'quality_scores_monomers.csv'

#modfold
mod_quality ={}
with open(csvf) as csvfile:
    reader = csv.reader(csvfile)
    headers = next(reader)
    headers = next(reader)
    for row in reader:
        aname = row[0]
        method_selected = row[1]
        phy_quality = row[5]
        sw_quality = row[4]
        rap_quality = row[3]
        int_quality = row[2]
        it_quality = row[6]
        feig_quality = row[7]
        galaxy_quality = row[8]
        mod_quality[aname]={'phy_quality':phy_quality, 
                   'sw_quality':sw_quality, 
                   'int_quality':int_quality, 
                   'rap_quality':rap_quality, 
                   'it_quality':it_quality, 
                   'feig_quality':feig_quality, 
                   'galaxy_quality':galaxy_quality}
print ('mod_quality > modfold scores for monomers')  

#voroMQA 
input_f =input_dir+'voromqa-result-4sep20.txt'
voro_quality ={}
with open(input_f) as f:
    next(f)
    lines = f.readlines()
    #for i in range(2):
    for line in lines:
        if lines.index(line)%2==0:
            pdbid=line.split('/')[2].split(' ')[0].strip()
            score= line.split('/')[2].split(' ')[1].strip()
            #print (pdbid, score)#, line.split(' ')[1].strip())
            voro_quality[pdbid]={'quality':score}
        if lines.index(line)%2:
            #print (line.strip())#.split(' ')[1])#(line.strip())
            continue
print ('voro_quality > voroMQA scores for complexes')

#**********************************  
# 1.1.4: MM-notes =  DICTIONARY WITH COMMENTS ON STRUCTURES 
#**********************************
#notes on biological assembly 
e = input_dir+'MM_notes.csv'                         
notes={}
with open(e) as csvfile2:
    spamreader2 = csv.reader(csvfile2)
    headers = next(spamreader2)
    for row in spamreader2:
        aname = row[0]
        new_notes = row[1]
        notes[aname] = {'MM-notes':new_notes}

#**********************************  
# 1.1.5: essentiality =  which proteins come from essential genes 
#**********************************         
essentiality={}
gene_data = json.load(open(input_dir+'genes_data.json', 'r'))
covert = json.load(open(input_dir+'protein_data.json', 'r'))
for ingredient in covert['data']:
    #1=essential, 0=non essential
    essentiality[ingredient['wid']]={'is_essential':''}
    if 'gene' in ingredient:
        #print (ingredient['wid'], ingredient['gene'])
        for i in gene_data['data']:
            if i['wid']==ingredient['gene']:
                essentiality[ingredient['wid']].update({'is_essential':i['is_essential']['value']})

#**********************************  
# 1.1.6: dsDNA binding 
#**********************************  
#binary dictionary with ONLY dsDNA binding proteins ==1, non binding dsDNA==0
#NB. elongation factors and ribosome do not bind dsDNA
                
DNAbind={}
dsDNA_bp_list, ds_ss_DNA_bp_dictionary = ds_ss_DNA_binding(input_dir+'protein_data.json')
for ingredient in covert['data']:
    #all ingredient have 0 in DNA binding unless they are listed in dsDNA_bp_list
    DNAbind[ingredient['wid']]={'binds_dsdna':'0'}
    if ingredient['wid'] in dsDNA_bp_list:
        DNAbind[ingredient['wid']].update({'binds_dsdna':'1'})  

#**********************************  
# 1.1.7: functional_category 
#********************************** 
#total of 31 categories 
#from CYT-MG model - doi: 10.1016/j.jmgm.2015.02.004
        
functional_category={} 
f = input_dir+'ingredient_function.csv'                         
with open(f) as csvfile:
    spamreader = csv.reader(csvfile)
    #headers = next(spamreader)
    for row in spamreader:
        aname = row[0]
        functional_category[aname] = {'function':row[1].strip()}
print ('functional_category > functions identified in MG-CYT model')

#**********************************  
# 1.1.8: functional_category_clustered
#********************************** 
##total of 12 categories
#clustering MG-CYT function classification in 12 macro categories

metabolism=['glycolysis', 'amino acid metabolism', 'lipid metabolism', 'cofactor metabolism', 'sugar metabolism',
'nucleotide metabolism']
protein_transport_singaling =['membrane transport', 'signaling']
DNA_replication_maintenance =['DNA degradation', 'DNA recombination','DNA remodeling/stabilization','DNA repair','foreign DNA processing','replication']
RNA_synthesis_maturation =['RNA remodeling','RNA degradation','RNA processing','aminoacyl tRNA synthetase']
protein_folding_maturation=['ribosome biogenesis','protein folding','protein degradation',
'post-translational processing',]
cytokinesis_motility=[ 'cytoskeleton', 'cell division', 'cell adhesion']
functional_category_cluster={}
for i in functional_category:
    functional_category_cluster[i]={'function':''}
    if functional_category[i]['function'] in metabolism:#1
        functional_category_cluster[i]['function']='metabolism'
    if functional_category[i]['function'] in protein_transport_singaling:#2
        functional_category_cluster[i]['function']='protein transport/singaling'
    if functional_category[i]['function'] in DNA_replication_maintenance:#3
        functional_category_cluster[i]['function']='DNA replication/maintenance'
    if functional_category[i]['function'] in RNA_synthesis_maturation:#4
        functional_category_cluster[i]['function']='RNA synthesis/maturation'
    if functional_category[i]['function'] in protein_folding_maturation:#5
        functional_category_cluster[i]['function']='protein folding/maturation'
    if functional_category[i]['function'] in cytokinesis_motility:#6
        functional_category_cluster[i]['function']='cytokinesis/motility'
    if functional_category[i]['function']=='host cell interaction':#7
        functional_category_cluster[i]['function']='host cell interaction'
    if functional_category[i]['function']=='MG-specific':#8
        functional_category_cluster[i]['function']='MG-specific'
    if functional_category[i]['function']=='transcription':#9
        functional_category_cluster[i]['function']='transcription'
    if functional_category[i]['function']=='translation':#10
        functional_category_cluster[i]['function']='translation'
    if functional_category[i]['function']=='lipoprotein':#11
        functional_category_cluster[i]['function']='lipoprotein'
    if functional_category[i]['function']=='unknown':#12
        functional_category_cluster[i]['function']='uncharacterized'

print ('functional_category_cluster > macro functional groups')
#-###################################################################################

                    # 1.2 : HOMOLOGY MODELS : homology models

#-###################################################################################

#**********************************  
# 1.2.1:  MONOMERS   
#**********************************
            
# #########################
# 1.2.1.1 : TEMPLATE MONOMERS 
# #########################

#spreadsheet with templates information from all the methods 
myworkbook = load_workbook(input_dir+'summary_templates_info.xlsx') #open file

# MAPPING PHYRE TEMPLATES
p_sheet = myworkbook['phyre_template'] #select sheet
phyre_template = {} # extract template from the file
for row in p_sheet.iter_rows(min_row=2):#skip the header
    #print (row)
    for cell in row:
        if cell.column_letter=='A':
            protId = cell.value
        if cell.column_letter=='E':
            t = cell.value
            template =([t[i:i+4] for i in range(0, len(t), 2)])[1]
            #print (t, template, protId)
            chain = [t[i:i+1] for i in range(0, len(t), 2)][3]        
            phyre_template[protId]=template+'_'+chain

# MAPPING SWISS MODEL TEMPLATES (for monomers and monomers part of homooligomers)
#Templates previously extracted from pdb files
sw_sheet = myworkbook['swiss_template'] #select sheet
swiss_template = {} # extract template from the file
for row in sw_sheet.iter_rows(min_row=2):#skip the header
    #print (row)
    for cell in row:
        if cell.column_letter=='A':
            protId = cell.value
        if cell.column_letter=='B':
            template = cell.value
        swiss_template[protId]=template

#MAPPING RAPTOR TEMPLATES     
raptor_template ={} #extract tempalte from file previously obtained 
r_sheet = myworkbook['raptor_template'] #select sheet
for row in r_sheet.iter_rows(min_row=2):#skip the header
    #print (row)
    for cell in row:
        if cell.column_letter=='A':
            protId = cell.value
            #print (protId)
        if cell.column_letter=='D':
            templ = str(cell.value)
        if cell.column_letter=='E':
            chain = str(cell.value)
            template = templ+'_'+chain
            if protId not in raptor_template:
                raptor_template[protId] = template
           #  when there are more than one template for the same monomer
            else :
                raptor_template[protId]+='; '+template  
        #raptor_template[protId]=template

#MAPPING INTFOLD TEMPLATES
intfold_template ={} #extract tempalte from file 
int_sheet = myworkbook['intfold_template'] #select sheet
for row in int_sheet.iter_rows(min_row=2):#skip the header
    for cell in row:
        if cell.column_letter=='A':
            protId = cell.value
        if cell.column_letter=='C':
            template = str(cell.value)
            if protId not in intfold_template:
                intfold_template[protId] = template
           #  when there are more than one template for the same monomer
            else :
                intfold_template[protId]+='; '+template

# MAPPING ITASSER TEMPLATES 
itasser_template ={}
it_sheet = myworkbook['itasser_template'] #select sheet
for row in it_sheet.iter_rows(min_row=2):#skip the header
    for cell in row:
        if cell.column_letter=='A':
            protId = cell.value
        if cell.column_letter=='B':
            template = cell.value
        itasser_template[protId]=template
        
#MAPPING MG-CYT MONOMERS TEMPLATES
#from CYT-MG model - doi: 10.1016/j.jmgm.2015.02.004     
c = input_dir+'CYT_MG_templates.csv'             
feig_template = {}
with open(c) as csvfile:
    spamreader = csv.reader(csvfile)
    headers = next(spamreader)
    for row in spamreader:
        gene = row[5]
        if gene !='':
            if len(gene)>=3:
                aname = 'MG_'+row[5]+'_MONOMER'
                template = row[9]
                feig_template[aname]=template
            if len(gene)==2:
                aname = 'MG_0'+row[5]+'_MONOMER'
                template = row[9]
                feig_template[aname]=template
            if len(gene)==1:
                aname = 'MG_00'+row[5]+'_MONOMER'
                template = row[9]
                feig_template[aname]=template

# MAPPING GALAXY TEMPLATE 
galaxy_template = {}
g_sheet = myworkbook['galaxy_template'] #select sheet
for row in g_sheet.iter_rows(min_row=2):#skip the header
    for cell in row:
        if cell.column_letter=='A':
            protId = cell.value
        if cell.column_letter=='D':
            template = cell.value
        galaxy_template[protId]=template

# #########################
# 1.2.1.2 : MAPPNIG MONOMERS 
# #########################
myworkbook = load_workbook(input_dir+'summary_monomers_info.xlsx') #open file
wc_workbook = load_workbook(input_dir+'WholeCellData.xlsx') #open file
# MAPPING INTFOLD
mapping_intfold ={}
int_sheet = myworkbook['intfold_models'] #select sheet
for row in int_sheet.iter_rows(min_row=2):#skip the header
    #print (row)
    for cell in row:
        if cell.column_letter=='A':
            protId = cell.value
        if cell.column_letter=='C':
            pdbfile = cell.value
            mapping_intfold[protId]={ 'original_fname': pdbfile, 'pdb_model':protId+'_int.pdb', 'chain':'', 'bu':'', 'template':'', 'quality':'','pcpalAxis':'', 'offset':'' }

#update chains and offsets
updateChainMembrane_new(mapping_intfold, dcm)
#update templates 
updateTemplates(mapping_intfold, intfold_template)
#update quality 
updateQuality(mapping_intfold, mod_quality, 'int_quality')
print('mapping_intfold')

#MAPPING ITASSER
mapping_itasser ={}
it_sheet = myworkbook['itasser_models'] #select sheet
for row in it_sheet.iter_rows(min_row=2):#skip the header
    #print (row)
    for cell in row:
        if cell.column_letter=='A':protId = cell.value
        if cell.column_letter=='B':pdbfile = cell.value
        mapping_itasser[protId]={ 'original_fname': pdbfile, 'pdb_model':protId+'_it.pdb', 'chain':'', 'bu':'', 'template':'', 'quality':'','pcpalAxis':'', 'offset':'' }

#update chains and offsets
updateChainMembrane_new(mapping_itasser, dcm)
#update templates 
updateTemplates(mapping_itasser, itasser_template)
#update quality 
updateQuality(mapping_itasser, mod_quality, 'it_quality')
print('mapping_itasser')  

## MAPPING RAPTOR
mapping_raptor={}
r_sheet = myworkbook['raptor_models'] #select sheet
for row in r_sheet.iter_rows(min_row=2):#skip the header
    #print (row)
    for cell in row:
        if cell.column_letter=='A':protId = cell.value
        if cell.column_letter=='B':pdbfile = cell.value
        mapping_raptor[protId]={ 'original_fname': pdbfile, 'pdb_model':protId+'_rap.pdb', 'chain':'', 'bu':'', 'template':'', 'quality':'','pcpalAxis':'', 'offset':'' }

#update chains and offsets
updateChainMembrane_new(mapping_raptor, dcm)
#update templates 
updateTemplates(mapping_raptor, raptor_template)
#update quality 
updateQuality(mapping_raptor, mod_quality, 'rap_quality')
print('mapping_raptor')

## MAPPING PHYRE
mapping_phyre={}
p_sheet = myworkbook['phyre_models'] #select sheet
for row in p_sheet.iter_rows(min_row=2):#skip the header
    #print (row)
    for cell in row:
        if cell.column_letter=='A':protId = cell.value
        if cell.column_letter=='D':
            pdbf = cell.value
            pdbfile = pdbf.replace(" ","")+".final.pdb"
            mapping_phyre[protId]={ 'original_fname': pdbfile, 'pdb_model':protId+'_phy.pdb', 'chain':'', 'bu':'', 'template':'', 'quality':'','pcpalAxis':'', 'offset':'' }

#update chains and offsets
updateChainMembrane_new(mapping_phyre, dcm)
#update templates 
updateTemplates(mapping_phyre, phyre_template)
#update quality 
updateQuality(mapping_phyre, mod_quality, 'phy_quality')
print('mapping_phyre')

#            mapping SWISS
mapping_swiss = {}
s_sheet = myworkbook['swiss_models'] #select sheet
for row in s_sheet.iter_rows(min_row=2):#skip the header
    #print (row)
    for cell in row:
        if cell.column_letter=='A':
            protId = cell.value
        if cell.column_letter=='B':
            pdbfile = cell.value
        if cell.column_letter=='C':
            original_fname = cell.value
            mapping_swiss[protId]={'original_fname': original_fname, 'pdb_model':pdbfile, 'chain':'', 'bu':'', 'template':'', 'quality':'','pcpalAxis':'', 'offset':'' }

#update chains and offsets
updateChainMembrane_new(mapping_swiss, dcm)
#update templates 
updateTemplates(mapping_swiss, swiss_template)
#update quality 
updateQuality(mapping_swiss, mod_quality, 'sw_quality')
print('mapping_swiss')    

#mapping monomers part of homomers from swiss    
mapping_swiss_homo = {}
sw_sheet = myworkbook['swiss_models_homomers'] #select sheet
for row in sw_sheet.iter_rows(min_row=2):#skip the header
    for cell in row:
        if cell.column_letter=='A':
            protId = cell.value
        if cell.column_letter=='D':
            original_fname = cell.value
        if cell.column_letter=='C':
            pdbfile = cell.value
            mapping_swiss_homo[protId]={ 'original_fname': original_fname, 'pdb_model':pdbfile, 'chain':'', 'bu':'', 'template':'', 'quality':'','pcpalAxis':'', 'offset':'' }

#update chains and offsets
updateChainMembrane_new(mapping_swiss_homo, dcm)
#update templates 
updateTemplates(mapping_swiss_homo, swiss_template)
#update quality 
updateQuality(mapping_swiss_homo, mod_quality, 'sw_quality')
print('mapping_swiss_homo')        

# mapping_CYT-MG MONOMERS
# includes the monomers part of complexes
mapping_feig={}
f_sheet = wc_workbook['S3M-Protein Monomers'] #select sheet
for row in f_sheet.iter_rows(min_row=2):#skip the header
    for cell in row:
        if cell.column_letter=='A':
            protId = cell.value
        if cell.column_letter=='J':
            if cell.value is not None:
                feig = str(cell.value)
                #print (feig)
                if ';' not in feig:
                    mapping_feig[protId]={'pdb_model':feig+'.pdb', 'chain':'', 'bu':'', 'template':'', 'quality':'', 'pcpalAxis':'', 'offset':''}    
                if ';' in feig:
                    mapping_feig[protId]={'1':{'pdb_model':feig.split(';')[0]+'.pdb', 'chain':'', 'bu':'', 'template':'', 'quality':'', 'pcpalAxis':'', 'offset':''}, 
                                         '2':{'pdb_model':feig.split(';')[1].strip()+'.pdb', 'chain':'', 'bu':'', 'template':'', 'quality':'', 'pcpalAxis':'', 'offset':''}}  
                
#update chains/offsets feig  
for i in mapping_feig:
    #one model
    if '1' not in mapping_feig[i]:
        f =  mapping_feig[i]['pdb_model']
        if i in dcm:
            for el in range(len(dcm[i])):
            #update dictionary if the same model is found in dc
                if f == dcm[i][el]['pdb_model']:
                    #print i, f, dc[i]['pdb_model'], dc[i]['chain'], dc[i]['bu']
                    mapping_feig[i].update({'chain':dcm[i][el]['chain'], 'bu':dcm[i][el]['bu'], 'offset':dcm[i][el]['offset'], 'pcpalAxis':dcm[i][el]['pcpalAxis']})
    #two models available
    if '1' in mapping_feig[i]:
        #print i, f
        f1 =  mapping_feig[i]['1']['pdb_model']
        f2 =  mapping_feig[i]['2']['pdb_model']
        if i in dcm:
            for el in range(len(dcm[i])):
                if f1 == dcm[i][el]['pdb_model']:
                    mapping_feig[i]['1'].update({'chain':dcm[i][el]['chain'], 'bu':dcm[i][el]['bu'], 'offset':dcm[i][el]['offset'], 'pcpalAxis':dcm[i][el]['pcpalAxis']})
                if f2 == dcm[i][el]['pdb_model']:
                    mapping_feig[i]['2'].update({'chain':dcm[i][el]['chain'], 'bu':dcm[i][el]['bu'], 'offset':dcm[i][el]['offset'], 'pcpalAxis':dcm[i][el]['pcpalAxis']})

#update templates 
updateTemplates(mapping_feig, feig_template)
for i in mapping_feig:
    #print i
    if '1' not in mapping_feig[i]:
        updateQuality(mapping_feig, mod_quality, 'feig_quality')
    if '1' in mapping_feig[i]:
        if i in mod_quality:
            mapping_feig[i]['1'].update({'quality':mod_quality[i]['feig_quality']})
            mapping_feig[i]['2'].update({'quality':mod_quality[i]['feig_quality']})
print('mapping_feig')

#               mapping GALAXY
mapping_galaxy = {}
g_sheet = myworkbook['galaxy_models'] #select sheet
for row in g_sheet.iter_rows(min_row=2):#skip the header
    #print (row)
    for cell in row:
        if cell.column_letter=='A':
            protId = cell.value
        if cell.column_letter=='B':
            pdbmodel = cell.value+'_g.pdb' #name of the complex
        if cell.column_letter=='E':
            original_fname = cell.value+'.pdb'
            mapping_galaxy[protId]={ 'original_fname': original_fname, 'pdb_model':pdbmodel, 'chain':'', 'bu':'', 'template':'', 'quality':'','pcpalAxis':'', 'offset':'' }

#update chains and offsets
updateChainMembrane_new(mapping_galaxy, dcm)
#update templates 
updateTemplates(mapping_galaxy, galaxy_template)
#update quality 
updateQuality(mapping_galaxy, mod_quality, 'galaxy_quality')
print('mapping_galaxy')

# #########################
# 1.2.1.3 :  all_mapping (monomers)
# #########################
#monomer list is found in S3M-Protein Monomers.csv, from which empty keys are removed: these are the monomers part of heterocomplexes 
all_mapping= {}
mono_sheet = wc_workbook['S3M-Protein Monomers'] #select sheet
for row in mono_sheet.iter_rows(min_row=2):#skip the header
    for cell in row:
        if cell.column_letter=='A':
            aname = cell.value
            all_mapping[aname]={}
            if aname in mapping_intfold:
                all_mapping[aname].update({'intfold':mapping_intfold[aname]})
            if aname in mapping_raptor:
                all_mapping[aname].update({'raptor':mapping_raptor[aname]})
            if aname in mapping_phyre:
                all_mapping[aname].update({'phyre':mapping_phyre[aname]})
            if aname in mapping_itasser:
                all_mapping[aname].update({'itasser':mapping_itasser[aname]})
            if aname in mapping_swiss:
                all_mapping[aname].update({'swiss':mapping_swiss[aname]})
            if aname in mapping_feig:
                if '1' not in mapping_feig[aname]:
                    all_mapping[aname].update({'feig':mapping_feig[aname]})
                if '1' in mapping_feig[aname]:
                    all_mapping[aname].update({'feig': {'1':mapping_feig[aname]['1'], '2':mapping_feig[aname]['2']}})#, mapping_feig[aname]['2']}})
            if aname in mapping_galaxy:
                all_mapping[aname].update({'galaxy':mapping_galaxy[aname]})
            if aname in mapping_swiss_homo:
                all_mapping[aname].update({'swiss':mapping_swiss_homo[aname]})
 
#remove the keys that do not have any homology modeled structure associated with them           
RemoveEmptyKeys(all_mapping)
print ('all_mapping')

#**********************************  
# 1.2.2:  COMPLEXES   
#**********************************
# #########################
# 1.2.2.1 : mapping COMPLEXES 
# #########################

#MAPPING GALAXY COMPLEXES
mapping_galaxy_compl = {}
g_sheet = myworkbook['galaxy_models'] #select sheet
for row in g_sheet.iter_rows(min_row=2):#skip the header
    #print (row)
    for cell in row:
        if cell.column_letter=='B':
            protId = cell.value # ie MG_XXX_DIMER
            pdbmodel = cell.value+'_g.pdb' #name of the complex
        if cell.column_letter=='E':original_fname = cell.value+'.pdb'
        mapping_galaxy_compl[protId]={ 'original_fname': original_fname, 'pdb_model':pdbmodel, 'chain':'', 'bu':'', 'template':'', 'quality':'','pcpalAxis':'', 'offset':'' }

#update chains and offsets
updateChainMembrane_new(mapping_galaxy_compl, dcm)
#update templates using mapping_template
for i in mapping_galaxy_compl:
    p = mapping_galaxy_compl[i]['pdb_model']
    for x in mapping_galaxy:
        if p in mapping_galaxy[x]['pdb_model']:
            mapping_galaxy_compl[i].update({'template':mapping_galaxy[x]['template']})        
#update the qulity with voroMQA scores
updateQualityComplexes(mapping_galaxy_compl, voro_quality)
print('mapping_galaxy_compl')

#MAPPING SWISS HOMOMERS
mapping_swiss_homo_compl = {}
sw_sheet = myworkbook['swiss_models_homomers'] #select sheet
for row in sw_sheet.iter_rows(min_row=2):#skip the header
    for cell in row:
        if cell.column_letter=='B':protId = cell.value
        if cell.column_letter=='D':original_fname = cell.value
        if cell.column_letter=='C':pdbfile = cell.value
        if cell.column_letter=='A':
            monomer_name = cell.value
            #print (monomer_name)
            mapping_swiss_homo_compl[protId]={'original_fname':original_fname, 'pdb_model':pdbfile, 'chain':'', 'bu':'', 'template':mapping_swiss_homo[monomer_name]['template'], 'quality':'', 'pcpalAxis':'', 'offset':''}
    
#update chains and offsets
updateChainMembrane_new(mapping_swiss_homo_compl, dcm)
#update the qulity with voroMQA scores
updateQualityComplexes(mapping_swiss_homo_compl, voro_quality)
print('mapping_swiss_homo_compl')

#CYT-MG COMPLEXES: HETERO AND HOMOCOMPLEXES
feig_complexes={}
f_sheet = wc_workbook['S3N-Macromolecular complexes'] #select sheet
for row in f_sheet.iter_rows(min_row=2):#skip the header
    for cell in row:
        if cell.column_letter=='A':
            protId = cell.value
        if cell.column_letter=='E':
            if cell.value is not None:
                feig = str(cell.value)
                #print (feig)
                if ';' not in feig:
                    feig_complexes[protId] = { 'pdb_model': feig+'.pdb', 'quality':'', 'pcpalAxis':'', 'offset':'', 'chain':'', 'bu':'', 'template':''}
                if ';' in feig: #for selecting a single pdb when it is like this: pros; prs0.pdb
                    feig_complexes[protId] = {'pdb_model': feig.split(';')[0]+'.pdb',  'quality':'', 'pcpalAxis':'', 'offset':'', 'chain':'', 'bu':'', 'template':''}

#update templates feig
a = input_dir+'CYT_MG_templates.csv'            
dft = {} #pdb_name: template
with open(a) as csvfile:
    line1 = csv.reader(csvfile)
    headers = next(line1)
    for row in line1:
        if row[1] !='':
            pdb_feig = row[1]
            #for those models that had the version RNA bound/unbound
            if ';' in pdb_feig:
                pdb_name = pdb_feig.split(';')[0]+'.pdb'
                template = row[9]
                dft[pdb_name]=template
            #grol.pdb is actually grol_noH.pdb
            if 'grol' in pdb_feig:
                pdb_name = pdb_feig+'_noH.pdb'
                template = row[9]
                dft[pdb_name]=template
            #for 'single' models, not RNA bound/unbound            
            if ';' not in pdb_feig:
                pdb_name = pdb_feig+'.pdb'
                template = row[9]
                dft[pdb_name]=template
for i in feig_complexes:
    f = feig_complexes[i]['pdb_model']
    if f in dft:
        #print i, f, dft[f]
        feig_complexes[i].update({'template':dft[f]})

#update chains and offsets
for i in feig_complexes:
    if '1' not in feig_complexes[i]:
        f =  feig_complexes[i]['pdb_model']
        if i in dcm:
            for el in range(len(dcm[i])):
            #update dictionary if the same model is found in dc
                if f == dcm[i][el]['pdb_model']:
                    #print i, f, dc[i]['pdb_model'], dc[i]['chain'], dc[i]['bu']
                    feig_complexes[i].update({'chain':dcm[i][el]['chain'], 'bu':dcm[i][el]['bu'], 'offset':dcm[i][el]['offset'], 'pcpalAxis':dcm[i][el]['pcpalAxis']})
    #two models 
    if '1' in feig_complexes[i]:
        #print i, f
        f1 =  feig_complexes[i]['1']['pdb_model']
        f2 =  feig_complexes[i]['2']['pdb_model']
        if i in dcm:
            for el in range(len(dcm[i])):
                if f1 == dcm[i][el]['pdb_model']:
                    feig_complexes[i]['1'].update({'chain':dcm[i][el]['chain'], 'bu':dcm[i][el]['bu'], 'offset':dcm[i][el]['offset'], 'pcpalAxis':dcm[i][el]['pcpalAxis']})
                if f2 == dcm[i][el]['pdb_model']:
                    feig_complexes[i]['2'].update({'chain':dcm[i][el]['chain'], 'bu':dcm[i][el]['bu'], 'offset':dcm[i][el]['offset'], 'pcpalAxis':dcm[i][el]['pcpalAxis']})
#update the qulity with voroMQA scores
updateQualityComplexes(feig_complexes, voro_quality)
print('feig_complexes')

# #########################
# 1.2.2.2 :     OUTPUT:  all_mapping_complex_homology
# #########################

all_mapping_complex_homology= {}
compl_sheet = wc_workbook['S3N-Macromolecular complexes'] #select sheet
for row in compl_sheet.iter_rows(min_row=2):#skip the header
    for cell in row:
        if cell.column_letter=='A':
            aname = cell.value
            all_mapping_complex_homology[aname]={}
            if aname in feig_complexes:
                all_mapping_complex_homology[aname].update({'feig':feig_complexes[aname]})
            if aname in mapping_swiss_homo_compl:
                all_mapping_complex_homology[aname].update({'swiss':mapping_swiss_homo_compl[aname]})
            if aname in mapping_galaxy_compl:
                all_mapping_complex_homology[aname].update({'galaxy':mapping_galaxy_compl[aname]})
    #this is for all MG_287_MONOMER_ACP, MG_287_MONOMER_ddcaACP etc complexes. I will assign the same structural model of MG_287_MONOMER
            if '287' in aname:
                #print (aname)
                all_mapping_complex_homology[aname].update(all_mapping['MG_287_MONOMER'])
    
#update 'ox' keys using the non-ox keys
for i in all_mapping_complex_homology:
    if 'ox' in i:
        a = i.split('_ox')[0]
        #print a, 'name'
        if a in all_mapping_complex_homology:
            all_mapping_complex_homology[i].update(all_mapping_complex_homology[a])
        if a in all_mapping:
            all_mapping_complex_homology[i].update(all_mapping[a])
            
#remove the keys that do not have any homology modeled structure associated with them           
RemoveEmptyKeys(all_mapping_complex_homology)
print('all_mapping_complex_homology')

#**********************************  
# 1.2.3:  OUTPUT: MONOMERS + COMPLEXES : homology_models 
#**********************************                                                      
homology_models = all_mapping                                                        
homology_models.update(all_mapping_complex_homology)

#-###################################################################################

                    # 1.3 : HOMOLOGS 

#-###################################################################################
                    
#**********************************  
# 1.3.1:  PDB-homologs dictionary
#**********************************  
exp_workbook = load_workbook(input_dir+'Experimental_Structures.xlsx') #open file
#CHAINS: first get them from csv_list, then, update if assinged in the dcm dictionary                      
homologs ={}
PDB_homologs={}
pdb_sheet = exp_workbook['PDB_homologs'] #select sheet
for row in pdb_sheet.iter_rows(min_row=2):#skip the header
    for cell in row:
        if cell.column_letter=='A': 
            protId= cell.value
        if cell.column_letter=='B': 
            pdb= cell.value
        if cell.column_letter=='C': 
            if cell.value is not None:
                chain= cell.value
            else:
                chain=''
        if cell.column_letter=='D': 
            bu= cell.value
        if cell.column_letter=='F': 
            pcpalAxis= cell.value
        if cell.column_letter=='G': 
            offset= cell.value
        if cell.column_letter=='E': 
            #surface= cell.value
            if cell.value is not None: #if true in surface use the axis in that file
                #print (protId, cell.value)
                surface= cell.value
                PDB_homologs[protId] = {'pdb_model': pdb, 'chain':chain, 'template':'-', 'bu':bu, 'quality':'', 'offset':offset, 'pcpalAxis':pcpalAxis}
            else:        
                PDB_homologs[protId] = {'pdb_model': pdb, 'chain':chain, 'template':'-', 'bu':bu, 'quality':'', 'offset': '', 'pcpalAxis': ''}
            homologs[protId]={}

#update chain and membrane
updateChainMembrane_new(PDB_homologs, dcm)

##update quality with voroMQA scores only for COMPLEXES
covert = json.load(open(input_dir+'protein_data.json', 'r'))
for ingredient in covert['data']: 
    #if it's a complex   
    if ingredient['model']=='ProteinComplex':
        complex_id = ingredient['wid']
        #if the complex id is in PDB-homologs
        if complex_id in PDB_homologs:
            pdb = PDB_homologs[complex_id]['pdb_model']+'.pdb'
            #if the structure is in voro_quality score dic
            if pdb in voro_quality:
                #print (complex_id, pdb, voro_quality[pdb]['quality'])
                PDB_homologs[complex_id].update({'quality':voro_quality[pdb]['quality']})

#**********************************  
# 1.3.2:  PDB-homologs-edited dictionary
#**********************************  
PDB_homologs_edited = {}
#monomers
het_sheet = exp_workbook['Heteromers'] #select sheet
for row in het_sheet.iter_rows(min_row=2):#skip the header
    for cell in row:
        if cell.column_letter=='B': protId= cell.value
        if cell.column_letter=='F': pdb_model_complete= cell.value
        if cell.column_letter=='E': 
            chain= ' or :'+str(cell.value)
        if cell.column_letter=='C': 
            pdb_model_single= cell.value
            #print (pdb_model_single)
            if pdb_model_single is not None and len(pdb_model_single)>4:
                #print (pdb_model_single)
                PDB_homologs_edited[protId]={'pdb_model':pdb_model_single, 'chain':chain, 'template':'-', 'bu':'', 'quality':'',  'offset':'', 'pcpalAxis':''}

#update chain and membrane
updateChainMembrane_new(PDB_homologs_edited, dcm)

#complexes
PDB_homologs_edited2 = {}
for row in het_sheet.iter_rows(min_row=2):#skip the header
    for cell in row:
        if cell.column_letter=='A': protId= cell.value
        if cell.column_letter=='D': chain= ' or :'+str(cell.value)
        if cell.column_letter=='F': 
            pdb_model_complete= cell.value
            if pdb_model_complete is not None and len(pdb_model_single)>4:
                PDB_homologs_edited2[protId]={'pdb_model':pdb_model_complete, 'chain':'', 'template':'-', 'bu':'', 'quality':'', 'offset':'', 'pcpalAxis':''}

#update chain and membrane
updateChainMembrane_new(PDB_homologs_edited2, dcm)
#merge monomers and complexes dictionary
PDB_homologs_edited.update(PDB_homologs_edited2)       
#update the qulity with voroMQA scores
updateQualityComplexes(PDB_homologs_edited2, voro_quality)

#**********************************  
# 1.3.3:  OUTPUT: homologs
#**********************************  
for i in homologs:
    if i in PDB_homologs:
        homologs[i].update({'PDB-homolog':PDB_homologs[i]})
    if i in PDB_homologs_edited:
        homologs[i].update({'PDB-homolog-edit':PDB_homologs_edited[i]})
        
#-###################################################################################

                    # 1.4: SOLVED: dictionary for solved structures 

#-###################################################################################
#**********************************  
# 1.4.1:  MAPPING_SOLVED
#**********************************                     
solved= {}
mapping_solved={}
solv_sheet = exp_workbook['MG_solved'] #select sheet
for row in solv_sheet.iter_rows(min_row=2):#skip the header
    for cell in row:
        if cell.column_letter=='A': 
            protId= cell.value
            solved[protId]={}
        if cell.column_letter=='B': 
            pdbfile= cell.value
        if cell.column_letter=='C': 
            chain= cell.value
        if cell.column_letter=='D': 
            other_solved_structures= cell.value
            mapping_solved[protId] = {'pdb_model':pdbfile,'chain':chain, 'bu':'', 'template':'-', 'quality':'','pcpalAxis':'', 'offset':'', 'other_pdbs': other_solved_structures} #'notes':notes ,      

#update chain and membrane if previously assinged
updateChainMembrane_new(mapping_solved, dcm)

#update quality with voroMQA scores only for COMPLEXES
covert = json.load(open(input_dir+'protein_data.json', 'r'))
for ingredient in covert['data']: 
    #if it's a complex   
    if ingredient['model']=='ProteinComplex':
        complex_id = ingredient['wid']
        #if the complex id is in PDB-homologs
        if complex_id in mapping_solved:
            pdb = mapping_solved[complex_id]['pdb_model']+'.pdb'
            #if the structure is in voro_quality score dic
            if pdb in voro_quality:
                #print (complex_id, pdb, voro_quality[pdb]['quality'])
                mapping_solved[complex_id].update({'quality':voro_quality[pdb]['quality']})

#**********************************  
# 1.4.2: output: solved
#**********************************     
#add 'solved'
for i in solved:
    if i in mapping_solved:
        solved[i].update({'solved':mapping_solved[i]})
print('solved')

#-###################################################################################

                    # 1.5: all_dict : homology models + homologs + solved

#-###################################################################################
                    
all_dict = homologs
for i in all_dict:
    if i in homology_models:
        #print i
        #l.append(i)
        all_dict[i].update(homology_models[i])
        #break
    if i in solved:
        #print i 
        all_dict[i].update(solved[i])
Dic2json(all_dict, 'all_dict')
print ('all dict')

#-###################################################################################

                    # 1.6:                 d1 & d2

#-###################################################################################
#**********************************  
# 1.6.1: d1 :  monomers
#**********************************
d1={}
function_updated = json.load(open(input_dir+'updated_function_dictionary.json', 'r'))
compartment_updated = json.load(open(input_dir+'compartment_updated.json', 'r'))
mono_sheet = wc_workbook['S3M-Protein Monomers'] 
for row in mono_sheet.iter_rows(min_row=2):#skip the header
    for cell in row:
        if cell.column_letter=='A': 
            protId= cell.value
        if cell.column_letter=='F': 
            seq= cell.value
        if cell.column_letter=='E': 
            leng= float(cell.value)
        if cell.column_letter=='G': 
            mw= cell.value
            d1[protId] = {'sequence':seq, 'length':leng, 'function': function_updated[protId]['function'], 'functional_category':functional_category_cluster[protId]['function'] , 'MM-notes': '', 'compartment':compartment_updated[protId]['compartment'], 'mw':mw} 

#complexation by data from covert
protdata = json.load(open(input_dir+'protein_data.json', 'r'))
#biosynthesis = {}
biosythesis_participants = {}
for i in  protdata['data']:   
    if i['model']=='ProteinMonomer':
        participation = []
        participation_empty = []
        #when monomers are not involved in any complex write '0'
        if len(i['protein_complex_biosythesis_participants'])== 0:
            biosythesis_participants[i['wid']]={'participation_in_complexes': participation_empty, 'model':'ProteinMonomer'}
        #when monomers are involved in complexes
        if len(i['protein_complex_biosythesis_participants'])!= 0:
            for c in range(len(i['protein_complex_biosythesis_participants'])):
                participation.append(i['protein_complex_biosythesis_participants'][c]['protein_complexes'])
                biosythesis_participants[i['wid']]={'participation_in_complexes': participation, 'model':'ProteinMonomer'}
#print ('biosythesis_participants > participation to complex yes/no')

UpdateDictionary(d1, hh)
UpdateDictionary(d1, uniprot)
UpdateDictionary(d1, biosythesis_participants)
UpdateDictionary(d1, notes)
UpdateDictionary(d1, essentiality)
UpdateDictionary(d1, DNAbind)

#**********************************  
# 1.6.2: d2 :  complexes
#**********************************
d2={}
compl_sheet = wc_workbook['S3N-Macromolecular complexes']
for row in compl_sheet.iter_rows(min_row=2):#skip the header
    for cell in row:
        if cell.column_letter=='A': 
            protId= cell.value
        if cell.column_letter=='B': 
            name= cell.value
        if cell.column_letter=='F': 
            mw= cell.value
            d2[protId] = {'function':name,  'functional_category':functional_category_cluster[protId]['function'], 'MM-notes':'', 'compartment': compartment_updated[protId]['compartment'], 'mw':mw}#, 'offset':'', 'pcpalAxis':''} 

#dictionary for composition and stoichiometrty
protdata = json.load(open(input_dir+'protein_data.json', 'r'))
#dictionary with the components and stochiometries
biosynthesis = {}
#dictionary with the list of the other complexes in which single complexes are involved
biosythesis_participants = {}
for i in  protdata['data']: 
    if i['model']=='ProteinComplex':
        #print 'ID: '+i['wid']
        mol = []
        coeff = []
        # molecules that make the complex
        for x in range(len(i['biosynthesis'])):
            mol.append(i['biosynthesis'][x]['molecule'])
            coeff.append(i['biosynthesis'][x]['coefficient'])
            biosynthesis[i['wid']]={'molecules': mol, 'coefficient':coeff}
        # other complexes in which a complex might be involved
        participation = []
        for c in range(len(i['protein_complex_biosythesis_participants'])):
            if len(i['protein_complex_biosythesis_participants'][c]['protein_complexes'])== 0:
                not_other_complex = []
                biosythesis_participants[i['wid']]={'participation_in_complexes': not_other_complex, 'model':'ProteinComplex'}
            else:
                participation.append(i['protein_complex_biosythesis_participants'][c]['protein_complexes'])
                biosythesis_participants[i['wid']]={'participation_in_complexes': participation, 'model':'ProteinComplex'}
                
#merge biosynthesis and biosythesis_participants
UpdateDictionary(biosynthesis, biosythesis_participants)               
#update d2
UpdateDictionary(d2, biosynthesis)                
print ('d2: + biosynthesis and components')    

# HHPRED complexes
HHpred1 ={} #HHpred1 = complex with not only monomers as components >>> HHpred left empty
HHpred2 ={} #HHpred2 = complex built by n monomers 

#HHpred1 = complex with not only monomers as components >>> HHpred left empty
#HHpred2 = complex built by n monomers 
#HHpred3 = complex built witho other complexes >>> 1) update d2 with HHpred2, 2) calculate HHpred using info from d2-updated

for aname in d2:
    hh_x_coeff =[]
    sum_coeff = []
    #for every element in 'molecules'
    for i in d2[aname]['molecules']:
        #if element in 'molecule' does not end in 'MONOMER' (to include the 'molecules' that are made by other complexes
        #if element in 'molecule' is not the name of the key
        #if element in 'molecule' is longer that 8, to exclude 'H', 'MGrrnA16S', MG_0001 and stuff like this 
        if i.endswith('MONOMER')==False and i!=aname and len(i)>8:
            #print aname, i 
            HHpred1[aname]={'HHpred':''}
   #if element in 'molecule' has not the same attribute listed above, do this
        else:
    # if this element is in d1 (meaning that must be a MG_XXX_MONOMER'):
            if i in d1:
                #get HHscore monomer part of the complex from d1 
                HHsingle = float(d1[i]['HHpred'])
                #get the position of a monomer in the list 'molecules' in d2
                index_molecules = d2[aname]['molecules'].index(i)
                #get the coefficient of the monomer searched above using its position in the list 'coefficient' in d2
                index_coeff = float(d2[aname]['coefficient'][index_molecules].split('-')[1])
                #multiply HHpred score with the coefficient
                hh_coeff_single = index_coeff*HHsingle
                #make a list with HHpred*coefficient for each monomer of a complex
                hh_x_coeff.append(hh_coeff_single)
                #make a list with all the coefficient for each monomer
                sum_coeff.append(index_coeff)
                #calculate the weighted average of HHpred scores
                weighted_average = sum(hh_x_coeff)/sum(sum_coeff) 
        if len(hh_x_coeff) !=0:
            HHpred2[aname]={'HHpred':weighted_average}
#update HHpred2 with HHpred1 : leave empty for those complexes that are made by other complexes and not only monomers             
HHpred2.update(HHpred1)
#update d2 with the HHpred data calculated so far
UpdateDictionary(d2, HHpred2)    
print ('d2: + HHpred2')     

HHpred3 ={}
for aname in d2:
    #quando HHpred e' vuoto:
    # if aname=='RIBOSOME_30S':
     if (d2[aname]['HHpred'] == ''):
        hh_x_coeff_d1 = []
        hh_x_coeff_d2 = []
        a = []
        c = []
        HH_d1 = 0
        index_molecules_d1 = 0
        index_coeff_d1 =0
        HH_d2 = 0
        index_molecules_d2 = 0
        index_coeff_d2 =0
        #print aname+'<<<<<'
    #per tutti gli elementi in molecules che 
        for i in d2[aname]['molecules']:
            if i in d1 and i!=aname:
                HH_d1 = float(d1[i]['HHpred'])
                index_molecules_d1 = d2[aname]['molecules'].index(i)
                index_coeff_d1 = float(d2[aname]['coefficient'][index_molecules_d1].split('-')[1])
                hh_coeff_single_d1 = index_coeff_d1*HH_d1
                a.append(hh_coeff_single_d1)
                c.append(index_coeff_d1)
                #print  i, hh_coeff_single_d1
    #questo se i componenti di un complesso sono complessi
            if i in d2 and i!=aname and d2[i]['HHpred']!='':
                HH_d2 = float(d2[i]['HHpred'])
                index_molecules_d2 = d2[aname]['molecules'].index(i)
             #get the coefficient of the monomer searched above using its position in the list 'coefficient' in d2
                index_coeff_d2 = float(d2[aname]['coefficient'][index_molecules_d2].split('-')[1])
             #multiply HHpred score with the coefficient
                hh_coeff_single_d2 = index_coeff_d2*HH_d2
             #make a list with HHpred*coefficient for each monomer of a complex
                a.append(hh_coeff_single_d2)
                c.append(index_coeff_d2)
                #print i, hh_coeff_single_d2
            # this i is a pain so i custom made its search
            if i=='MG_224_9MER_GTP':
                i='MG_224_MONOMER'
                HH_d1 = float(d1[i]['HHpred'])
                hh_coeff_single_d1 = 1*HH_d1
                a.append(hh_coeff_single_d1)
                c.append(1) 
        #all the i that do not respect the previous criteria
            else:
                continue   
        #calculate the weighted average of HHpred scores
        if sum(c)!=0:
            weighted_average = sum(a)/sum(c)
            #print weighted_average
            HHpred3[aname]={'HHpred':weighted_average}
        if sum(c)==0:
            HHpred3[aname]={'HHpred':''}
            #print aname
        if aname=='RIBOSOME_30S':
           # print HHpred3['RIBOSOME_30S_IF3']['HHpred']
            hh_x_coeff_d1 = []
            l = []
            m = []
            for i in d2[aname]['molecules']:
                if i in d1 and i!=aname:
                    HH_d1 = float(d1[i]['HHpred'])
                    index_molecules_d1 = d2[aname]['molecules'].index(i)
                    index_coeff_d1 = float(d2[aname]['coefficient'][index_molecules_d1].split('-')[1])
                    hh_coeff_single_d1 = index_coeff_d1*HH_d1
                    l.append(hh_coeff_single_d1)
                    m.append(index_coeff_d1)
                    weighted_average2 = (sum(l)+float(d1['MG_196_MONOMER']['HHpred']))/(21)
                    HHpred3['RIBOSOME_30S_IF3']={'HHpred':weighted_average2}
#updare for RIBOSOME 70, need a second round                     
for aname in HHpred3:
   if aname=='RIBOSOME_70S':#        
       HHpred3[aname].update({'HHpred':(HHpred3['RIBOSOME_30S']['HHpred']+HHpred3['RIBOSOME_50S']['HHpred'])/2})
print ('d2: + HHpred3')    

UpdateDictionary(d2, HHpred3)
UpdateDictionary(d2, notes)
UpdateDictionary(d1, essentiality)
UpdateDictionary(d2, DNAbind)

#-###################################################################################

                    # 1.7:                 d3
                    
#-###################################################################################
#**********************************  
# 1.7.1:  draft of d3
#**********************************
#make d3
monomer = d1
multimer = d2
#merge dic
d3 = dict(monomer)
d3.update(multimer)

#assign 'method' in each ingredient of d3
seq_file= input_dir+'method_selection.csv'
with open(seq_file) as csvfile:
    spamreader = csv.reader(csvfile)
    headers = next(spamreader)
    for row in spamreader:
        aname = row[0]
        method = row[1]
#uncomment this for a normal recipe 
        d3[aname].update({'method':method})
        
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
#this is for when you want a recipe only with the homology models derived from homologous structures 
#no PDB files are generated since all the pdb models come from the PDB
#        d3[aname].update({'method':'PDB-homolog'})
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#**********************************  
# 1.7.2: update d3 with other dictionaries based on the method specified 
#**********************************            
# 0 = phyre 1=intfold 2=raptor 3=swiss 4=galaxyy 5=itasser 6 =pdb-homolog 7 = pdb-homolog-edit 8 =solved 9 =feig 10 = no method 
method = ['phyre', 'intfold', 'raptor', 'swiss', 'galaxy', 'itasser', 'PDB-homolog', 'PDB-homolog-edit', 'solved', 'feig','' ]
    
for i in d3: 
    #print (i)
    if d3[i]['method']==method[0]:
        d3[i].update(all_dict[i][method[0]])
        #shutil.copy(all_dict[i][method[0]]['dir'], PDBs_folder+os.sep+all_dict[i][method[0]]['pdb_model'])
    if d3[i]['method']==method[1]:
        d3[i].update(all_dict[i][method[1]])
        #shutil.copy(all_dict[i][method[1]]['dir'], PDBs_folder+os.sep+all_dict[i][method[1]]['pdb_model'])
    if d3[i]['method']==method[2]:
        d3[i].update(all_dict[i][method[2]])
        #shutil.copy(all_dict[i][method[2]]['dir'], PDBs_folder+os.sep+all_dict[i][method[2]]['pdb_model'])
    if d3[i]['method']==method[3]:
        d3[i].update(all_dict[i][method[3]])
        #shutil.copy(all_dict[i][method[3]]['dir'], PDBs_folder+os.sep+all_dict[i][method[3]]['pdb_model'])
    if d3[i]['method']==method[4]:
        d3[i].update(all_dict[i][method[4]])
        #shutil.copy(all_dict[i][method[4]]['dir'], PDBs_folder+os.sep+all_dict[i][method[4]]['pdb_model'])
    if d3[i]['method']==method[5]:
        d3[i].update(all_dict[i][method[5]])
        #shutil.copy(all_dict[i][method[5]]['dir'], PDBs_folder+os.sep+all_dict[i][method[5]]['pdb_model'])
    if d3[i]['method']==method[6]:
        d3[i].update(all_dict[i][method[6]])
    if d3[i]['method']==method[7]:
        d3[i].update(all_dict[i][method[7]])
        #shutil.copy(all_dict[i][method[7]]['dir'], PDBs_folder+os.sep+all_dict[i][method[7]]['pdb_model'])
    if d3[i]['method']==method[8]:
        d3[i].update(all_dict[i][method[8]])
    if d3[i]['method']==method[9]:
        #some feig monomers have 2 models 
        if len(all_dict[i][method[9]])==2:
            d3[i].update(all_dict[i][method[9]]['1'])
            #shutil.copy(all_dict[i][method[9]]['1']['dir'], PDBs_folder+os.sep+all_dict[i][method[9]]['1']['pdb_model'])
        else:# len(all_dict[i][method[9]])==2:
            d3[i].update(all_dict[i][method[9]])
            #shutil.copy(all_dict[i][method[9]]['dir'], PDBs_folder+os.sep+all_dict[i][method[9]]['pdb_model'])
    if d3[i]['method']==method[10]:
        d3[i].update({'bu': '','chain': '', 'offset': '','pcpalAxis': '','pdb_model': '','quality': '', 'template':''})

jsonf = 'd3.json'
f= open(output_dir+jsonf,"w")
f.write(json.dumps(d3))
f.close()
print ('d3')   

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@     2 : COPY NUMBER     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
frame_mono =149 #1184 #6973       frames_m json =[150, 1185, 6974] 
frame_complex =150 #1189 #6960    frames_c json =[146, 1190, 6961]
#set the array for  a specific simulation specified above 'sim'

array_names_mono, array_names_compl, data_mono, data_compl = set_np_arrays(sim)

l =[]
for i in covert['data']:
    l.append(i['wid'])

#this is an absolutely outrageous number of lines for the task is doing. 
#However today i am tired and even if this is slow af i will keep it and call the day. amen

states_c = []
for i in array_names_compl:
    states_c.append(i.decode().split('-')[1])
states_c = set(states_c)
#buond proteins count is defined in the nucleoid recipe so here we do not consider them 

states_m = []
for i in array_names_mono:
    states_m.append(i.decode().split('R-')[1])
states_m = set(states_m)
#buond proteins count is defined in the nucleoid recipe so here we do not consider them 

copy_numb = {}   
for i in l:   
    #initialize everything with 0
    copy_numb[i]={'count':0, 'compartment':0} 
    for state in states_c:
        aname = i+'-'+state
        #complexes
        if aname.encode() in array_names_compl: 
            position = np.where(array_names_compl == aname.encode())[0][0]
            molecule_per_compartment = []
            for index in range(6):
                count=data_compl[position, index , frame_complex]
                molecule_per_compartment.append(count)
                #for all states that are inactive, processes etc take only the one >0
                if count!=0 and state!='mature':
                    #print (aname, index, count)
                    #molecules in cytoplasm
                    if molecule_per_compartment.index(max(molecule_per_compartment))==float('0'):
                        copy_numb[aname]={'count':max(molecule_per_compartment), 'compartment':'c'}
                        #print (aname, max(molecule_per_compartment))
                    #molecules in DNA compartment
                    if molecule_per_compartment.index(max(molecule_per_compartment))==float('1'):
                        copy_numb[aname]={'count':max(molecule_per_compartment), 'compartment':'d'}
                        #print (aname)
                    #molecules in extracellular space
                    if molecule_per_compartment.index(max(molecule_per_compartment))==float('2'):
                        copy_numb[aname]={'count':max(molecule_per_compartment), 'compartment':'e'}
                        #print (aname)
                    #molecules in membrane 
                    if molecule_per_compartment.index(max(molecule_per_compartment))==float('3'):
                        copy_numb[aname]={'count':max(molecule_per_compartment), 'compartment':'m'}
                        #print (aname)
                    #molecules in transmembrane organelle cytoplasm
                    if molecule_per_compartment.index(max(molecule_per_compartment))==float('4'):
                        copy_numb[aname]={'count':max(molecule_per_compartment), 'compartment':'tc'}
                        #print (aname)
                    #molecules in transmembrane organelle membrane 
                    if molecule_per_compartment.index(max(molecule_per_compartment))==float('5'):
                        copy_numb[aname]={'count':max(molecule_per_compartment), 'compartment':'tm'}
                        #print (aname)
                #mature proteins do not get the 'matrue' tag in the dictionary key 
                if state=='mature':
                                    #molecules in cytoplasm
                    if molecule_per_compartment.index(max(molecule_per_compartment))==float('0'):
                        copy_numb[i]={'count':max(molecule_per_compartment), 'compartment':'c'}
                        #print (aname, max(molecule_per_compartment))
                    #molecules in DNA compartment
                    if molecule_per_compartment.index(max(molecule_per_compartment))==float('1'):
                        copy_numb[i]={'count':max(molecule_per_compartment), 'compartment':'d'}
                        #print (aname)
                    #molecules in extracellular space
                    if molecule_per_compartment.index(max(molecule_per_compartment))==float('2'):
                        copy_numb[i]={'count':max(molecule_per_compartment), 'compartment':'e'}
                        #print (aname)
                    #molecules in membrane 
                    if molecule_per_compartment.index(max(molecule_per_compartment))==float('3'):
                        copy_numb[i]={'count':max(molecule_per_compartment), 'compartment':'m'}
                        #print (aname)
                    #molecules in transmembrane organelle cytoplasm
                    if molecule_per_compartment.index(max(molecule_per_compartment))==float('4'):
                        copy_numb[i]={'count':max(molecule_per_compartment), 'compartment':'tc'}
                        #print (aname)
                    #molecules in transmembrane organelle membrane 
                    if molecule_per_compartment.index(max(molecule_per_compartment))==float('5'):
                        copy_numb[i]={'count':max(molecule_per_compartment), 'compartment':'tm'}
                        #print (aname)continue
   #MONOMERS
    for state in states_m:
        aname = i+'-'+state 
        if aname.encode() in array_names_mono: 
            position = np.where(array_names_mono == aname.encode())[0][0]
            molecule_per_compartment = []
            for index in range(6):
                count=data_mono[frame_mono, index , position]
                molecule_per_compartment.append(count)
                if count!=0 and state!='mature':
                    if molecule_per_compartment.index(max(molecule_per_compartment))==float('0'):
                        copy_numb[aname]={'count':max(molecule_per_compartment), 'compartment':'c'}
                        #print (aname, max(molecule_per_compartment))
                    if molecule_per_compartment.index(max(molecule_per_compartment))==float('1'):
                        copy_numb[aname]={'count':max(molecule_per_compartment), 'compartment':'d'}
                        #print (aname, max(molecule_per_compartment))
                    if molecule_per_compartment.index(max(molecule_per_compartment))==float('2'):
                        copy_numb[aname]={'count':max(molecule_per_compartment), 'compartment':'e'}
                        #print (aname, max(molecule_per_compartment))
                    if molecule_per_compartment.index(max(molecule_per_compartment))==float('3'):
                        copy_numb[aname]={'count':max(molecule_per_compartment), 'compartment':'m'}
                        #print (aname, max(molecule_per_compartment))
                    if molecule_per_compartment.index(max(molecule_per_compartment))==float('4'):
                        copy_numb[aname]={'count':max(molecule_per_compartment), 'compartment':'tc'}
                        #print (aname, max(molecule_per_compartment))
                    if molecule_per_compartment.index(max(molecule_per_compartment))==float('5'):
                        copy_numb[aname]={'count':max(molecule_per_compartment), 'compartment':'tm'}
                        #print (aname, max(molecule_per_compartment))
                if state=='mature':
                    if molecule_per_compartment.index(max(molecule_per_compartment))==float('0'):
                        copy_numb[i]={'count':max(molecule_per_compartment), 'compartment':'c'}
                    if molecule_per_compartment.index(max(molecule_per_compartment))==float('1'):
                        copy_numb[i]={'count':max(molecule_per_compartment), 'compartment':'d'}
                    if molecule_per_compartment.index(max(molecule_per_compartment))==float('2'):
                        copy_numb[i]={'count':max(molecule_per_compartment), 'compartment':'e'}
                    if molecule_per_compartment.index(max(molecule_per_compartment))==float('3'):
                        copy_numb[i]={'count':max(molecule_per_compartment), 'compartment':'m'}
                    if molecule_per_compartment.index(max(molecule_per_compartment))==float('4'):
                        copy_numb[i]={'count':max(molecule_per_compartment), 'compartment':'tc'}
                    if molecule_per_compartment.index(max(molecule_per_compartment))==float('5'):
                        copy_numb[i]={'count':max(molecule_per_compartment), 'compartment':'tm'}

#this sums all the copy number of the same ingredient in the same id
#it DOES NOT conut the DNA binding prot bound copy number to its mature one 
#the copy number here considers all protein states apart from DNA bindign protein 'bound'  
copy_numb = edit_copy_numb_cluster_all_states(copy_numb)
UpdateDictionary(copy_numb, d3)    
for i in copy_numb:
    if '-'  in i:
        aname = i.split('-')[0]
        copy_numb[i].update(d3[aname])  
        
d3=copy_numb
print ('d3 w copy number')

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@     3 : RNA dictionary      @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#list of rna ingredients
rna_list =[]
rRNA_l = ['MGrrnA16S', 'MGrrnA23S', 'MGrrnA5S']
tRNA_l = ['tRNA', 'tRNA-aminoacylated']#all tRNAs will be collapsed in two categories: mature, aminoacylated
sRNA_l = ['MG_0001', 'MG_0002', 'MG_0003', 'MG_0004']
rna_list= rRNA_l+tRNA_l+sRNA_l

#count RNAs
#TOTAL NUMBER OF RNAS not only free
tRNA, rRNA, sRNA, mRNA = free_rna_count(sim, frame_complex)

#sum of all tRNAs 
tRNA_aminmoacylated =[]
tRNA_others =[] #mature, processed, damaged etc
sum_tRNA=[]
for i in tRNA:
    if 'aminmoacylated' in i[1]:
        tRNA_aminmoacylated.append(i[2])
    else: 
        tRNA_others.append(i[2])
sum_tRNA_aminmoacylated =sum(tRNA_aminmoacylated)
sum_tRNA_others=sum(tRNA_others)
sum_tRNA=[sum_tRNA_others,sum_tRNA_aminmoacylated]
#sum_trna = sum(tRNA['f2']) #sum of all tRNAs

#could be done faster but i suck 
MG_0001=[]
MG_0002=[]
MG_0003=[]
MG_0004=[]
for i in sRNA:
    if 'MG_0001' in i[1]:
        MG_0001.append(i[2])
    if 'MG_0002' in i[1]:
        MG_0002.append(i[2])
    if 'MG_0003' in i[1]:
        MG_0003.append(i[2])
    if 'MG_0004' in i[1]:
        MG_0004.append(i[2])
sum_MG_0001=sum(MG_0001)
sum_MG_0002=sum(MG_0002)
sum_MG_0003=sum(MG_0003)
sum_MG_0004=sum(MG_0004)
sum_sRNA = [sum_MG_0001,sum_MG_0002, sum_MG_0003, sum_MG_0004] 

#sum of all rRNAs
#sum_rrna = sum(rRNA['f2']) not sure why i did this before 

MGrrnA16S=[]
MGrrnA23S=[]
MGrrnA5S=[]
for i in rRNA:
    if 'MGrrnA16S' in i[1]:
        MGrrnA16S.append(i[2])
    if 'MGrrnA23S' in i[1]:
        MGrrnA23S.append(i[2])
    if 'MGrrnA5S' in i[1]:
        MGrrnA5S.append(i[2])
sum_MGrrnA16S=sum(MGrrnA16S)
sum_MGrrnA23S=sum(MGrrnA23S)
sum_MGrrnA5S=sum(MGrrnA5S)
sum_rRNA = [sum_MGrrnA16S,sum_MGrrnA23S, sum_MGrrnA5S] 

#mRNA (for other purposes)
sum_mRNA=[]
for i in mRNA:
    sum_mRNA.append(i[2])

#first draft

rna ={}
for i in rna_list:
    rna[i]={'function':'', 'MM-notes':'',  'compartment':'c', 'model': 'RNA',  
       'method': '', 'pdb_model': '', 'chain': '', 'bu': 'AU',  
       'offset': '', 'pcpalAxis': '', 'count':'', 'function_category':''}

#rRNAs 
rRNA_chains = [' or :AA', ' or :BB', ' or :BA']
for i in rRNA_l:
    rna[i].update({'function': 'rRNA', 'MM-notes': 'structure from E. coli ribosome structure 4V69',
       'method': 'PDB-homolog', 'pdb_model': '4V69', 'chain': rRNA_chains[rRNA_l.index(i)], 'bu': 'AU', 'count':sum_rRNA[rRNA_l.index(i)], 'function_category':'translation'})
#tRNA
for i in tRNA_l:
    rna[i].update({'function': 'tRNA', 'MM-notes': 'Same structure has been used for all tRNAs', 'method': 'PDB-homolog',
 'pdb_model': '6TNA', 'chain':'', 'bu': 'AU', 'count':sum_tRNA[tRNA_l.index(i)], 'function_category':'translation'})

#sRNAs
sRNA_pdb=['4WFM', '4V2S', '3DHS', '2CZJ']
sRNA_chiains=[' or :A', ' or :Q', ' or :A', ' or :B']
sRNA_function=['RNA synthesis/maturation', 'RNA synthesis/maturation', 'RNA synthesis/maturation', 'translation']
for i in sRNA_l:
    rna[i].update({'function': 'sRNA', 'MM-notes': 'selected structure from the PDB', 'method': 'PDB-homolog', 
 'pdb_model': sRNA_pdb[sRNA_l.index(i)], 'chain':sRNA_chiains[sRNA_l.index(i)], 'bu': 'AU', 'count':sum_sRNA[sRNA_l.index(i)], 'function_category':sRNA_function[sRNA_l.index(i)]})
    

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@     3 : recipe in csv format      @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

make_csv_recipe(d3,rna, 'root')

