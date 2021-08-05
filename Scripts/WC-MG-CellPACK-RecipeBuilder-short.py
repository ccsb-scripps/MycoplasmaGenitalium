# -*- coding: utf-8 -*-
"""
@author: MM
Created on Wed Jul 28 16:08:03 2021
ASSEMBLES A MG DRAFT RECIPE COMPATIBLE WITH MESOSCOPE
"""
exec(open('WC-MG-CellPACK-functions.py').read())
#input_dir ='C:\\Users\\marti\\Documents\\wholecellworkflow\\FOR GIT-cellPACKpaper\\input_files\\'
#input_dir ='input_files'+os.sep
input_dir = 'C:/Users/marti/Documents/wholecellworkflow/FOR GIT-cellPACKpaper/input_files/'
output_dir= 'C:/Users/marti/Documents/WC-model/scripts_output_junk/'
#output_folder = 'scripts_output'+os.sep

#1. make d3: final recipe as a dictionary
#    1.1 general info on ingredients
#        essentiality
#        DNA binding
#        functional_category 
#        functional_category_clustered  
#    1.2 d1 (dictionary monomers) & d2 (dictionary complexes)
#        d1 (monomers)
#        d2 (complexes)
#    1.3 d3: d1 + d2
#3. COPY NUMBER
#2. RNA (tRNA+rRNA_sRNA)
#4. recipe as csv file

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@      
## 1.1 : GENERAL INFO ON INGREDIENTS        
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
#**********************************  
# uniprot
#**********************************
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
##print ('uniprot > uniprot monomers')
#**********************************  
# HHscore
#**********************************
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
#**********************************  
# essentiality =  which proteins come from essential genes 
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
#dsDNA binding 
#**********************************  
# dictionary with  dsDNA binding proteins ==1, non binding dsDNA==0
#NB. elongation factors and ribosome do not bind dsDNA       
DNAbind={}
dsDNA_bp_list, ds_ss_DNA_bp_dictionary = ds_ss_DNA_binding(input_dir+'protein_data.json')
for ingredient in covert['data']:
    #all ingredient have 0 in DNA binding unless they are listed in dsDNA_bp_list
    DNAbind[ingredient['wid']]={'binds_dsdna':'0'}
    if ingredient['wid'] in dsDNA_bp_list:
        DNAbind[ingredient['wid']].update({'binds_dsdna':'1'})  
#**********************************  
#functional_category 
#********************************** 
#total of 31 categories from CYT-MG model - doi: 10.1016/j.jmgm.2015.02.004  
functional_category={} 
f = input_dir+'ingredient_function.csv'                         
with open(f) as csvfile:
    spamreader = csv.reader(csvfile)
    #headers = next(spamreader)
    for row in spamreader:
        aname = row[0]
        functional_category[aname] = {'function':row[1].strip()}
#print ('functional_category > functions identified in MG-CYT model')
#**********************************  
#functional_category_clustered
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

#-###################################################################################
# 2. d1 & d2
#-###################################################################################
#this file holds info about all the availabe structural models 
all_dict = json.load(open(input_dir+'all_dict.json', 'r'))

#**********************************  
# d1 :  monomers
#**********************************
d1={}
function_updated = json.load(open(input_dir+'updated_function_dictionary.json', 'r'))
compartment_updated = json.load(open(input_dir+'compartment_updated.json', 'r'))
wc_workbook = load_workbook(input_dir+'WholeCellData.xlsx') #open file
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
            d1[protId] = {'sequence':seq, 'length':leng, 'function': function_updated[protId]['function'], 'functional_category':functional_category_cluster[protId]['function'] , 'compartment':compartment_updated[protId]['compartment'], 'mw':mw} 
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
UpdateDictionary(d1, essentiality)
UpdateDictionary(d1, DNAbind)
#**********************************  
#d2 :  complexes
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
            d2[protId] = {'function':name,  'functional_category':functional_category_cluster[protId]['function'], 'compartment': compartment_updated[protId]['compartment'], 'mw':mw}#, 'offset':'', 'pcpalAxis':''} 
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
#print ('d2: + biosynthesis and components')    
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
#print ('d2: + HHpred2')     
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
#print ('d2: + HHpred3')    
UpdateDictionary(d2, HHpred3)
UpdateDictionary(d1, essentiality)
UpdateDictionary(d2, DNAbind)
#-###################################################################################
#d3                    
#-###################################################################################
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
 
#update d3 with other dictionaries based on the method specified             
method = ['phyre', 'intfold', 'raptor', 'swiss', 'galaxy', 'itasser', 'PDB-homolog', 'PDB-homolog-edit', 'solved', 'feig','' ]
for i in d3: 
    #print (i)
    if d3[i]['method']==method[0]:
        d3[i].update(all_dict[i][method[0]])
    if d3[i]['method']==method[1]:
        d3[i].update(all_dict[i][method[1]])
    if d3[i]['method']==method[2]:
        d3[i].update(all_dict[i][method[2]])
    if d3[i]['method']==method[3]:
        d3[i].update(all_dict[i][method[3]])
    if d3[i]['method']==method[4]:
        d3[i].update(all_dict[i][method[4]])
    if d3[i]['method']==method[5]:
        d3[i].update(all_dict[i][method[5]])
    if d3[i]['method']==method[6]:
        d3[i].update(all_dict[i][method[6]])
    if d3[i]['method']==method[7]:
        d3[i].update(all_dict[i][method[7]])
    if d3[i]['method']==method[8]:
        d3[i].update(all_dict[i][method[8]])
    if d3[i]['method']==method[9]:
        if len(all_dict[i][method[9]])==2:
            d3[i].update(all_dict[i][method[9]]['1'])
        else:
            d3[i].update(all_dict[i][method[9]])
    if d3[i]['method']==method[10]:
        d3[i].update({'bu': '','chain': '', 'offset': '','pcpalAxis': '','pdb_model': '','quality': '', 'template':''})
#jsonf = 'd3.json'
#f= open(output_dir+jsonf,"w")
#f.write(json.dumps(d3))
#f.close()
#print ('d3')   
        
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@      
## 2 : COPY NUMBER        
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
frame_mono =149 #1184 #6973       frames_m json =[150, 1185, 6974] 
frame_complex =150 #1189 #6960    frames_c json =[146, 1190, 6961]
#set the array for  a specific simulation specified above 'sim'
array_names_mono, array_names_compl, data_mono, data_compl = set_np_arrays(sim)
l =[]
for i in covert['data']:
    l.append(i['wid'])
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
#it DOES NOT count the DNA binding prot bound copy number to its mature one 
#the copy number here considers all protein states apart from DNA bindign protein 'bound'  
copy_numb = edit_copy_numb_cluster_all_states(copy_numb)
UpdateDictionary(copy_numb, d3)    
for i in copy_numb:
    if '-'  in i:
        aname = i.split('-')[0]
        copy_numb[i].update(d3[aname])  
        
d3=copy_numb
#print ('d3 w copy number')
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@      
## 3 : RNA dictionary        
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 

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
    rna[i]={'function':'', 'compartment':'c', 'model': 'RNA',  
       'method': '', 'pdb_model': '', 'chain': '', 'bu': 'AU',  
       'offset': '', 'pcpalAxis': '', 'count':'', 'function_category':''}
#rRNAs 
rRNA_chains = [' or :AA', ' or :BB', ' or :BA']
for i in rRNA_l:
    rna[i].update({'function': 'rRNA',
       'method': 'PDB-homolog', 'pdb_model': '4V69', 'chain': rRNA_chains[rRNA_l.index(i)], 'bu': 'AU', 'count':sum_rRNA[rRNA_l.index(i)], 'function_category':'translation'})
#tRNA
for i in tRNA_l:
    rna[i].update({'function': 'tRNA', 'method': 'PDB-homolog',
 'pdb_model': '6TNA', 'chain':'', 'bu': 'AU', 'count':sum_tRNA[tRNA_l.index(i)], 'function_category':'translation'})
#sRNAs
sRNA_pdb=['4WFM', '4V2S', '3DHS', '2CZJ']
sRNA_chiains=[' or :A', ' or :Q', ' or :A', ' or :B']
sRNA_function=['RNA synthesis/maturation', 'RNA synthesis/maturation', 'RNA synthesis/maturation', 'translation']
for i in sRNA_l:
    rna[i].update({'function': 'sRNA',  'method': 'PDB-homolog', 
 'pdb_model': sRNA_pdb[sRNA_l.index(i)], 'chain':sRNA_chiains[sRNA_l.index(i)], 'bu': 'AU', 'count':sum_sRNA[sRNA_l.index(i)], 'function_category':sRNA_function[sRNA_l.index(i)]})

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@      
## 4: RECIPE IN CSV FORMAT  
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
make_csv_recipe(d3,rna, 'root')

