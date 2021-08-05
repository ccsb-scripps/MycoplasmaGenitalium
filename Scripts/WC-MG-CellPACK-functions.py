# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 13:14:28 2021

@author: MM
"""
import re
import os
import json
import h5py
import numpy as np
import csv
import shutil
import h5py
from openpyxl import load_workbook

#----------------------------------------------------------------------------
# FUNCTIONS TO GET DATA FROM WC-SIMULATIONS and ASSEMBLE A DRAFT RECIPE REDABLE ON MESOSCOPE
#----------------------------------------------------------------------------
#INPUT FILES INFO
    #to run all the following code you need:
    # 1) protein database and gene database from cover wholecellkb: 'protein_data.json' | 'genes_data.json'
    # 2) simulation file .h5 format: simulation-1.h5 (downloaded from http://www.wholecellsimdb.org/simulation_batch/1))
    # 3) json file with positions of DNA binding monomers/complexes: getSeriesData_monomers1.json | getSeriesData_complex1.json
    # downloaded from: http://www.wholecellviz.org/getSeriesData.php?sim_id=2011_10_19_02_53_45_1&class_name=Chromosome&attr_name=monomerBoundSites
    # and http://www.wholecellviz.org/getSeriesData.php?sim_id=2011_10_19_02_53_45_1&class_name=Chromosome&attr_name=complexBoundSites 
    # 4) frame number (note that when extracting data from getSeriesData_monomers1.json | getSeriesData_complex1.json or sim 
    #the frame number is different due to the file type but they refer to the same frame in the simulation.
    #Complexes and monomers do not have the same time stamps because 
    
#frames json: 150(mono)/146(complex); 1185(mono)/1190 (complex); 6974 (mono)/6961(complex)
#frames when working with the hdf5 files is: (frame json) -1
#frames sim: 149(mono)/145(complex); 1184(mono)/1189(complex); 6973(mono)/6960(complex)

#input files:
input_dir ='input_files'+os.sep

#outfput folder:
output_dir = output_folder = 'scripts_output'+os.sep

#holds DNA binding prot positions for monomers 
fname_m = input_dir+'getSeriesData_monomers1.json'
#holds DNA binding prot positions for complexes 
fname_c = input_dir+'getSeriesData_complexes1.json'
#simulation #1-1
sim = input_dir+'simulation-1.h5' 

#TRANSCRIPTION UNIT
    #creates a dictionary where every monomer is associated with its transcription unit, gene length, coords, etc 
    #attention: the function uses some actual file and refers to actual file locations (ie csvf, input_dir, 'genes_data.json')
    
    #you need: 
    #       hdf5 file
    #       beware of 'S3K-Transcription units.csv'
    #       'protein_data.json' downloaded from http://www.wholecellkb.org/
    #       'genes_data.json' downloaded from http://www.wholecellkb.org/
    
def prepare_dictionary_index(sim):  
    x = h5py.File(sim)
    #create the variable for the labels 
    labels_mono= x.get("states/ProteinMonomer/counts/labels/0")
    #make a numpy array with the monomer labels 
    array_mono = np.array(labels_mono)
    #make initial dictionary index to protein dictionary and to gene ie 482: {'protID': 'MG_470_MONOMER'}
    index_dict = {} 
    for i in array_mono:
        if '-nascent'.encode() in i:
            index = np.where(array_mono == i)[0][0]+1
            #nb. i use decode because of p3 issues with ascii 
            index_dict[index]={'protID':i.decode().split('-')[0]}    
    #make a ditionary for trancription units ie. 'TU_335': {'type': 'polycistronic', 'length_TU': 2766, 'start_tu': 577268}
    tu = {}
    #input is a csv file 
    csvf = input_dir+'S3K-Transcription units.csv'
    with open(csvf) as csvfile:
        csvreader = csv.reader(csvfile)
        header = next(csvreader)
        header = next(csvreader)
        header = next(csvreader)
        for row in csvreader:
            TU = row[0] #TU column 
            gen = row[3] #genes column
            leng = row[5] #length genes column 
            start_coord = row[4] #start coordinate column 
            product = row[2]
            if 'mRNA' in row[2]:
                if ',' in gen: #',' means there are multiple genes in the same TU -> polycistronic
                    tu[TU]={'type':'polycistronic', 'length_TU':int(leng)+1, 'start_tu':int(start_coord), 'product':product}
                else: #no ',' means a single gene in the TU -> monocistronic
                    tu[TU]={'type':'monocistronic', 'length_TU':int(leng)+1, 'start_tu':int(start_coord), 'product':product}
            else:
                if ',' in gen:
                    tu[TU]={'type':'polycistronic', 'length_TU':int(leng)+1, 'start_tu':int(start_coord), 'product':product}
                else:
                    tu[TU]={'type':'monocistronic', 'length_TU':int(leng)+1, 'start_tu':int(start_coord), 'product':product}
    genes = json.load(open(input_dir+'genes_data.json', 'r'))    
    #make initial dictionary index for genes 
    gene_dic={}
    for agene in genes['data']:
        gene_dic[agene['wid']]={'TU':agene['transcription_units'][0], 'gene_length':agene['length'], 'type':agene['type'][0], 'length_TU':''}
    for i in gene_dic:
        if gene_dic[i]['TU'] in tu:
            gene_dic[i].update({'length_TU':tu[gene_dic[i]['TU']]['length_TU']})
    #update the index_dict to include tu, gene len, coord, direction, type etc
    #ie. 482: {'protID': 'MG_470_MONOMER', 'gene': 'MG_470', 'gene_length': 810, 'TU': 'TU_335',  'coord_gene': 579224,
    #          'direction': 'Reverse', 'type': 'polycistronic', 'length_TU': 2766, 'coord_tu': 577268}
    covert = json.load(open(input_dir+'protein_data.json', 'r'))
    for i in index_dict:
        name = index_dict[i]['protID']
        for aname in covert['data']:
            if name == aname['wid']:
                index_dict[i].update({'gene':aname['gene']})
                break
        gene = index_dict[i]['gene']
        for agene in genes['data']:
            if gene == agene['wid']:
                index_dict[i].update({'gene_length':int(agene['length']), 'TU':agene['transcription_units'][0], 'coord_gene':int(agene['coordinate']), 'direction':agene['direction']})
        unit = index_dict[i]['TU']
        for u in tu:
            if unit == u:
                index_dict[i].update({'type':tu[u]['type'],'length_TU':tu[u]['length_TU'], 'coord_tu':tu[u]['start_tu']})
    return index_dict, gene_dic

#usage
index_dict, genes_dict =prepare_dictionary_index(sim)

   
def dna_binding_list(json_data):
    covert = json.load(open(input_dir+json_data, 'r'))
    l =[]
    d={}
    for i in covert['data']:
        if i['model']=='ProteinMonomer':
            if i['dna_footprint']!= None :
            #print i['dna_footprint']
                l.append(i['wid'])
                d[i['wid']]={'dna_footprint':i['dna_footprint']['length']}
        if i['model']=='ProteinComplex':
            if i['dna_footprint']!= None :
           #print i['dna_footprint']
               l.append(i['wid'])
               d[i['wid']]={'dna_footprint':i['dna_footprint']['length']}
    return l, d
##usage
##this is a list of all DNA binding proteins in the simulation
dna_bp_list, dna_bp_dictionary = dna_binding_list('protein_data.json')

#PROTEIN COORD CHROMOSOME
    #makes a csv file with the proteinid+position+strand from json file      
    #NB if you want this to work, write the PRECISE frame number found in the json file 
    #this is based on json so the frame you use here corresponds to -1 in the sim
    
    #you need: hdf5 file (sim)
    #         json file for monomer bound to dna (getSeriesData_monomers1.json)
    #         json file for complexes bound to dna (getSeriesData_complexes1.json)
    #         frame number for monomers 
    #         frame number for complexes 
    
def chromosome_proteins(fname_m, frame_number_m, fname_c, frame_number_c,sim):
    z = []
    csv = open(output_folder+'chromosome_proteins_'+str(frame_number_m-1)+'.csv', 'w')
    csv.write('#DNA_BINDING PROTEINS\n')
    csv.write('MONOMERS FRAME '+' '+str(frame_number_m-1)+'\n')
    csv.write('COMPLEXES FRAME '+' '+str(frame_number_c-1)+'\n')
    csv.write('protId, coordinates, strand, dna_footprint\n')
    #open json file with list of dictionaries: each dictionary represents a DNA-bound protein molecule at an instance of simulation time
    mono = json.load(open(fname_m, 'r'))
    compl = json.load(open(fname_c, 'r'))
    #open the simulation file format hdf5
    x = h5py.File(sim)
    #array for the monomers nasme <> identification number
    labels_mono= x.get("states/ProteinMonomer/counts/labels/0")
    labels_compl = x.get("states/ProteinComplex/counts/labels/0")
    #make a list with the frame number
    frame_m = [frame_number_m]
    frame_c= [frame_number_c]
    #for a specific frame look into the json file and find: 
    for i in range(len(mono)):
        if mono[i]['time']==frame_m:
            for el in range(len(mono[i]['pos'])):
                index = mono[i]['val'][el]-1 #-1 BECAUSE IT COMES FROM MATLAB
                pos = str(mono[i]['pos'][el])
                strand = str(mono[i]['strand'][el])
                #print index, labels_mono[index], pos, strand
                #print mono[i]['val'][el]-1, labels_mono[index], mono[i]['pos'][el], mono[i]['strand'][el]
                protid=labels_mono[index].decode().split('-nascent')[0]
                footprint= dna_bp_dictionary[protid]['dna_footprint']
                csv.write(protid+','+pos+','+strand+','+footprint+'\n')
    for i in range(len(compl)):
        if compl[i]['time']==frame_c:
            for el in range(len(compl[i]['pos'])):
                index = compl[i]['val'][el]-1  #-1 BECAUSE IT COMES FROM MATLAB
                pos = str(compl[i]['pos'][el])
                strand = str(compl[i]['strand'][el])      
                protid=labels_compl[index].decode().split('-nascent')[0]
                footprint= dna_bp_dictionary[protid]['dna_footprint']
                #print compl[i]['val'][el], labels_compl[index], compl[i]['pos'][el], compl[i]['strand'][el] 
                #print labels_compl[index], pos, strand
                if 'RNA_POLYMERASE' not in labels_compl[index].decode():
                    csv.write(labels_compl[index].decode().split('-nascent')[0]+','+pos+','+strand+','+footprint+'\n')
                if 'RNA_POLYMERASE' in labels_compl[index].decode():
                    z.append(labels_compl[index])
                    pos = str(compl[i]['pos'][el]+37) #+37 because it is half of the rnapoly foot print of 75. It was encoded differently in the h5 and json files so adding +37 gives the numbers of the h5
                    csv.write(labels_compl[index].decode().split('-nascent')[0]+','+pos+','+strand+','+footprint+'\n')
            print('rna polymerases number', len(z))
    csv.close()
    
#RNA POLYMERASE
    #to run this use the number of the json -1 
    #make a csv file starting from the simulation file and listing the RNApoly count 
    #in the state column the actively transcribing polymerases are indicated by any number different from 0, -1, -2, -3
    
    #you need: hdf5 file (sim)
    #          frame
    #          outputfolder  
    
def rna_poly(sim, frame):
    x = h5py.File(sim)
    active = x.get("states/RNAPolymerase/nActive/data")
    free = x.get("states/RNAPolymerase/nFree/data")
    promoter = x.get("states/RNAPolymerase/nSpecificallyBound/data")
    nonspecificallybound = x.get("states/RNAPolymerase/nNonSpecificallyBound/data")
    #define the total number of rna_poly
    tot_rna_poly = active[0,0,frame]+free[0,0,frame]+promoter[0,0,frame]+nonspecificallybound[0,0,frame]    
    #print (active[0,0,frame], free[0,0,frame], promoter[0,0,frame], nonspecificallybound[0,0,frame], (tot_rna_poly))
    #position and strand 
    pos_data = x.get("states/RNAPolymerase/positionStrands/data")
    #transcribing/ free / non specifically bound / promoter bound. ie. state_data[0][0][500]=state_data[protid][state][frame]
    state_data = x.get("states/RNAPolymerase/states/data")
    #sotre the length of transcript
    TranscriptProgress_data = x.get("states/Transcript/boundTranscriptProgress/data")   
    #stores the transcription unit index that is being transcribed. 
    tu_data = x.get("states/Transcript/boundTranscriptionUnits/data")
    #make an array of 5 columns with integers full of only 0
    #a = np.zeros(shape=(pos_data.shape[0], 5), dtype=np.int)
    a = np.zeros(shape=(tot_rna_poly, 5), dtype=np.int)
    #update the array witht the two columns of coordinate and strand (pos_data)
    a[:,:2] = pos_data[:tot_rna_poly,:,frame]
    #fill the 3rd column with state_data. NB (-3=>promoter bound, -2=>free, -1=>non-specifically bound, 0=>not exist, 1=>actively transcribing)
    a[:,2] = np.array(state_data[:tot_rna_poly,:,frame]).flatten()
    #fill the 4th column with length of transcript TranscriptProgress_data
    a[:,3] = np.array(TranscriptProgress_data[:tot_rna_poly,:,frame]).flatten()
    #fill the 5th column with the transcription unit index (0 to 334) tu_data
    a[:,4] = np.array(tu_data[:tot_rna_poly,:,frame]).flatten()
    #write a csv file with this array
    np.savetxt(output_folder+'RNApoly_'+str(frame)+'.csv', a, 
               delimiter=",", fmt="%d", 
               header='coordinate (nt), strand, state, transcript length (nt), TranscriptionUnit id', 
               comments='#coordinate: refers to the START of RNApoly\n #state: -3=>promoter bound; -2=>free; -1=>non-specifically bound; 0=>not exist; 1=>actively transcribing\n #TU not for mRNA=> TU_088: rRNA; TU_178/TU_181/TU_226/TU_227: sRNA;TU_007/TU_037/TU_123/TU_134/TU_135/TU_136/TU_137/TU_138/TU_139/TU_140/TU_141/TU_146/TU_166/TU_168/TU_172/TU_187/TU_188/TU_189/TU_194/TU_220/TU_221/TU_223/TU_248/TU_250/TU_252/TU_276: tRNA\n ')
    #comments='in state: -3=>promoter bound, -2=>free, -1=>non-specifically bound, 0=>not exist, everything =or> 1 =>actively transcribing')    
#    #this below is just to count the polymerases
#    poly = []
#    promoter = []
#    free = []
#    nonspecific = []
#    active = []
#    for i in pos_data[:, :, frame]:
#        coordinate = i[0]
#        strand = i[1]
#        if strand !=0:
#            poly.append(strand)
#            #print 'RNA-POLYMERASE', coordinate, strand
#            #csv.write('RNA_POLYMERASE'+','+str(coordinate)+','+str(strand)+'\n')
#    for state in state_data[:, :, frame]:
#       if state[0]==[-1]:
#           nonspecific.append(state[0])
#       if state[0]==[-3]:
#           promoter.append(state[0])
#       if state[0]==[-2]:
#           free.append(state[0])
#       #every state that is not 0, -1, -2, -3 is ACTIVELY TRANSCRIBING
#       if state[0]!=[0] and state[0]!=[-2] and state[0]!=[-1] and state[0]!=[-3]:
#           active.append(state[0])        
    print('tot polymerases number:', (tot_rna_poly))
    print ('active', (active[0,0,frame]), 'free', (free[0,0,frame]), 'promoter', (promoter[0,0,frame]), 'Non specifically bound', (nonspecificallybound[0,0,frame]))
#    print('actively transcribing (1 or 0): ', len(active))
#    print('free (-2 or 3): ', len(free))
#    print('non specific bound(-1 or 2) :', len(nonspecific))
#    print('promoter bound (-3 or 1): ', len(promoter))    

#RNAs
    #output 1 (now commented):csv wit a list of sRNA and mRNA 
    #output 2:arrays of mRNA, tRNA, sRNA, rRNA
    #TOTAL NUMBER OF RNAS, NOT FREE RNAS

def free_rna_count(sim, frame):
    d, genes_dict =prepare_dictionary_index(sim)
    f = h5py.File(sim)
    labels0 = f.get("states/Rna/counts/labels/0")
    labels_tu = np.array(labels0)
    data_rna = f.get("states/Rna/counts/data")
    #array with only compartment 0 at a specific frame
    z =data_rna[frame, 0]
        #data for the non zero RNA molecules
    non_zero =z[z>0]   
        #find the id for the RNA molecules > 0 and make and array
    idx_non_zero = np.where(data_rna[frame, 0]>0) #this is a tuple 
    idx_RNA = np.array(idx_non_zero)
        #make a new array with zero
    a=np.zeros(non_zero.shape,
                   dtype = 'i4, U20, i4, i4, U20')
        #column 0: RNA ids 
    a['f0'] = np.array(idx_RNA[:]).flatten()
        #column 2: count of RNA molecules 
    a['f2'] = np.array(non_zero[:]).flatten()
    for i,idx in enumerate(a['f0']):
        #column 1: name of RNA molecule
        a['f1'][i]= labels_tu[idx]
        rna_id =labels_tu[idx].decode().split('-')[0]
        #print (rna_id)
        if rna_id in genes_dict:#for rRNA, tRNA, sRNA
            #column 3: gene len
            #column 4: RNA type
            a['f3'][i]= genes_dict[rna_id]['gene_length']
            a['f4'][i]= genes_dict[rna_id]['type']
        if rna_id not in genes_dict:  #mRNA 
            #print (rna_id)
            for x in genes_dict:
                if rna_id==genes_dict[x]['TU']:
                    a['f3'][i]= genes_dict[x]['length_TU']
                    a['f4'][i]= genes_dict[x]['type'] #reorder the array on mRNA, tRNA, rRNA, sRNA column
        #a_sorted = np.sort(a, order=['f3'], axis=0)
        #depending on what you want you can create csv with the different types of free RNAs
    mRNA = a[a['f4']=='mRNA']
    sRNA = a[a['f4']=='sRNA']
    rRNA = a[a['f4']=='rRNA']
    tRNA_1 = a[a['f4']=='tRNA'] #careful here because there are also initiator_tRNA
    tRNA_2 = a[a['f4']=='initiator_tRNA']
    tRNA = np.concatenate((tRNA_1, tRNA_2))
    arr = np.concatenate((sRNA, mRNA))
    arr_ordered = np.sort(arr, order=['f4'])
    return tRNA, rRNA, sRNA, mRNA

#RIBOSOME 
    #to run this use the number of the json -1 
    #make a csv file starting from the simulation file and listing the ribosomes position, directions, states etc
    #you need: 
    #       hdf5 file (sim)
    #       frame
    #outputfolder
    #run ribosome 
    #run prepare_dictionary_index(sim) 
    
#(PART 1)
    #make an array with ribosome positions, SHIFTED positons, directions, states etc
def ribosome(sim, frame, output_folder):
    d, genes_dict =prepare_dictionary_index(sim)
    x = h5py.File(sim)
    #ribosome state
    active = x.get("states/Ribosome/nActive/data")
    stalled = x.get("states/Ribosome/nStalled/data")
    nonexist = x.get("states/Ribosome/nNotExist/data")
    #define the total number of ribosomes and take this info from the simulation avoiding the 0
    tot_rib = active[0,0,frame]+nonexist[0,0,frame]+stalled[0,0,frame]
    #-1=>stalled, 0=>not exist, 1=>actively translating
    state = x.get("states/Ribosome/states/data")
    #gene that is being transcribed /  transcript bound to a ribosome / 
    #Karr: protein coding gene indices indicate the indices of the proteins in the array of proteins.
    boundRNA = x.get("states/Ribosome/boundMRNAs/data")
    #mRNA position in codons of each 70s ribosome from the start codon. equal to the length of the polypeptide
    mRNAPosition = x.get("states/Ribosome/mRNAPositions/data")   
    #position of the ribosome (in codons) along the tmRNA template
    tmRNAPosition = x.get("states/Ribosome/tmRNAPositions/data")
    #make an array of 4 columns with integers full of only 0
    #a = np.zeros(boundRNA.shape[0], 
    #add a 15th column for the number of Ribosome bound
    a=np.zeros(tot_rib,
                 dtype = 'U20, i4, i4, U20, U20, i4, i4, U20, i4, i4, U20, U20, i4, i4, i4, i4, i4')
                 #dtype =i4 for numbers, S20=strings 
    #update the array. 1st column state, 2nd column gene translated, 3rd column rib position on mRNA, 4th column rib 
    # 0 state
    a['f1'] = np.array(state[:tot_rib,:,frame]).flatten()
    # 1 boundmRNA -  prot index
    a['f2'] = np.array(boundRNA[:tot_rib,:,frame]).flatten()
    # 2 ribosome position in codons 
    a['f12'] = np.array(mRNAPosition[:tot_rib,:,frame]).flatten()
    # 3 ribosome position in codons if ribosome stalled
    a['f14'] = np.array(tmRNAPosition[:tot_rib,:,frame]).flatten()    
    for i, idx in enumerate(a['f2']):
        a['f0'][i] = 'RIBOSOME'
        try:
            # 4 gene name 
            a['f4'][i] = d[idx]["gene"]
            # 5 protein name dtype S20
            a['f3'][i] = d[idx]["protID"]
            # gene length 
            a['f6'][i] = d[idx]['gene_length']
            # transcription unit name 
            a['f7'][i] = d[idx]['TU']
            # length transcription unit
            a['f8'][i] = d[idx]['length_TU']
            # mono/polycistronic
            a['f11'][i] = d[idx]['type']
            # 11 direction
            a['f10'][i] = d[idx]['direction']
            # starting coord forward gene /end reverse genes 
            a['f5'][i] = d[idx]['coord_gene']
            # starting coordinate TU
            a['f9'][i] = d[idx]['coord_tu']
        except:
            #print (i,idx)
            continue     
    #messing around in the next line for debugging
    #a_sorted =a
    a_sorted = np.sort(a, order=['f2'], axis=0)
    #rinv is the mRNA id
    unq,rinv, count = np.unique(a_sorted, axis=0,return_inverse=True, return_counts=True)
    #makes another array with only the lines that are repeated > 1
    repeated_groups = unq[count > 1]
    #for each line that has duplicates make an array with the indices of those lines in a
    for repeated_group in repeated_groups:
        #this 'continue' is to make work ribosome_mRNA functions
        #so it DOES NOT check for repeted postions in this moment
        continue
        repeated_idx = np.argwhere(a_sorted == repeated_group) 
        #print (repeated_idx)
        for i in range(len(repeated_idx)): 
            n=20 #number of codons to shift the ribosome on the same gene. Produces a 60nt shift in ribosome position
            if a_sorted['f1'][repeated_idx[i]]!=0: #if the state is not 'non existent'
                #calculate the new ribosome position by adding 20*0, 20*1, 20*2... 
                #i is the index so n*20=0 => the position in the first repeated line will not change
                new_position=a_sorted['f12'][repeated_idx[i]]+n*i
                #calculated the ribosome positin in nt => we need this to check that the ribosome position is not longer than the gene
                new_position_in_nt = new_position*3
                #if the ribosome new ribosome position isnt longer than the gene 
                if new_position_in_nt<a_sorted['f6'][repeated_idx[i]]: #gene length
                    a_sorted['f12'][repeated_idx[i]]=new_position 
                    #print ('new_position '+str(new_position))
                #if the previous condition fail the ribosome position 
                else:
                    print ('error ribosome shift')
                    #shift back 20 codons
                    a_sorted['f12'][repeated_idx[i]]=a_sorted['f12'][repeated_idx[i]]-n*i
    #mrna id
    a_sorted['f16'] = rinv#id
    # 10 ribosome position in nt =  position in codons * 3
    a_sorted['f13'] = a_sorted['f12']*3
    #### RIBOSOME POSITION onto fulle length mRNA CALCULATIONS    
    #ribosome fw=  (coord gene -  coord tu) + rib position
    i = np.where(a_sorted['f10'] == 'Forward')
    #a['f14']==RIBOSOME POSITION
    a_sorted['f15'][i] = (a_sorted['f5'][i] - a_sorted['f9'][i]) + a_sorted['f13'][i]
    #ribosome rv = [(tu coord + tu length) - (gene coord + gene length)] +rib position
    z = np.where(a_sorted['f10'] =='Reverse')
    a_sorted['f15'][z] = ((a_sorted['f9'][z] + a_sorted['f8'][z]) - (a_sorted['f5'][z]+ a_sorted['f6'][z])) + a_sorted['f13'][z]  
    return a_sorted    
    
#(PART 2)
def Ribosome_RNAs(sim,frame):
    a = ribosome(sim, frame, output_folder)
    #array removing non-existing ribosomes 
    non_zero_rib = a[a['f1']!=0]
    #f0='RIBOSOME', f16= mRNAid, f6 = length Gene, f7= TU, f8= length TU, 
    #f14=tmRNA position - I USE IT ONLY BECAUSE IT IS A LIST OF 0 f15=ribosome position
    bound = non_zero_rib[['f0', 'f16', 'f6','f7', 'f8', 'f15', 'f14']]
    for i in bound:
        #if the ribosome is too close to the beginning of the transcript OR it was in postion 0
        if i[5]<30: 
            idx = np.argwhere(bound == i)
            print ('hey here ribosome position was too close to the beginning, 30nt shift applied', idx, i[5])
            #add 30nt
            bound['f15'][idx[0][0]]=bound['f15'][idx[0][0]]+30
        #if the position + 30 is longer than the actual TU
        if i[5]+30 > i[4]:
            difference=(i[5]+30)-i[4]
            idx = np.argwhere(bound == i)
            #print (i[5]+30, i[4], difference)
            bound['f15'][idx[0][0]]=bound['f15'][idx[0][0]]- difference
            print ('hey here ribosome position was too close to the end of the mRNA, shifted back', idx, i[5])
        #replace f14 (column of 0s) with a copy of the ribosome positions
    bound['f14']=non_zero_rib['f15']
    nRIB= len(bound)
    #ARRAY FOR all  MRNAs
    tRNA, rRNA, sRNA, mRNA = free_rna_count(sim, frame)  
    mRNA['f1'] = [n[:6] for n in mRNA['f1']] 
    #array that repeats the TU whose copy number is >1, number of repetitions is in f2=count
    rep_mRNA = np.repeat(mRNA,mRNA['f2'])
    nRNA= len(rep_mRNA)
    #label=mRNA / mRNAid  / TU / mRNA length / ribosome position 
    tot_mRNA = np.zeros((nRIB+nRNA,), dtype = 'U20, i4, U20, i4, i4, i4')#, i4')
    indices = np.arange(0,nRNA)
    tot_mRNA['f0'][:nRNA]= rep_mRNA['f4'] #label: 'mRNA'
    tot_mRNA['f1'][:nRNA]= indices
    tot_mRNA['f2'][:nRNA]= rep_mRNA['f1'] #TU
    tot_mRNA['f3'][:nRNA]= rep_mRNA['f3'] #mRNA length
    tot_mRNA['f4'][:nRNA]= rep_mRNA['f2'] #count --- this is temporary
    mRNA_lastpos={}
    #RIBOSOMES
    for i in range(nRIB):
        #i=ribosome id
        tu = bound['f7'][i]
        mRNA_ids = np.where(tot_mRNA['f2'] == tu)[0]
        #pick one randomly --------------------------NB this means it can change
        mRNA_id = np.random.choice(mRNA_ids,1)[0]
        bound['f16'][i] = mRNA_id
        bound['f6'][i]=i
        #bound = np.sort(bound, order=['f16'], axis=0)
        #mRNA_lastpos[mRNA_id] = bound['f14'][i]
        if (mRNA_id in mRNA_lastpos) :
            #calculate the difference
            diff = bound['f15'][i] - mRNA_lastpos[mRNA_id]  #0 same, 60 expected different
            #print (i, diff)
            if abs(diff) <= 60 : #if the difference is <60
                #if shifting the ribosome is still (f15 ribosome position)< than gene length (f8) -30
                if (mRNA_lastpos[mRNA_id] + 60)<(bound['f8'][i]-30):
                    bound['f14'][i] = mRNA_lastpos[mRNA_id] + 60
                    #rint (i, mRNA_id, bound['f15'][i], mRNA_lastpos[mRNA_id], mRNA_lastpos[mRNA_id] + 60)
                else: #if shifting the ribosome over the gene length (f8)
                    print ('error ribosome shift', i, bound['f7'][i] , bound['f15'][i])
                   #new position for the ribosome at least 60nt before the end of the mRNA
                    bound['f14'][i] = mRNA_lastpos[mRNA_id] - 60
                    #print (mRNA_lastpos[mRNA_id] - 60)
            else:#if the difference is > or = 60 keep the position
                bound['f14'][i] = bound['f15'][i]
        else :#if (mRNA_id NOT in mRNA_lastpos) :
            bound['f14'][i] = bound['f15'][i]
        mRNA_lastpos[mRNA_id] = bound['f14'][i] #last position on that given mRNA
        #break
    indices = np.arange(0,nRIB)
    tot_mRNA['f0'][nRNA:]= bound['f0'] #label: 'RIBOSOME'
    tot_mRNA['f1'][nRNA:]= indices
    tot_mRNA['f2'][nRNA:]= bound['f7'] #TU
    tot_mRNA['f3'][nRNA:]= bound['f8'] #TU length
    tot_mRNA['f4'][nRNA:]= bound['f16'] #mRNA id
    #tot_mRNA['f5'][nRNA:]= bound['f15'] #original position
    tot_mRNA['f5'][nRNA:]= bound ['f14'] #new position 
    np.savetxt(output_folder+'mRNA_ribosome_'+str(frame)+'.csv', tot_mRNA, 
                                  #'U20, i4, U20, i4, i4, i4'
               delimiter=",",  fmt='%s,%d,%s,%d,%d, %d', 
               header='#name, id, TU, TU length, copy id, ribosome_position (in RIBOSOME)')    
#CELL SIZE
def cell_size(sim, frame):
    f = h5py.File(sim)
    #shape of states/Geometry
    data_geom = f.get('states/Geometry/volume/data')
    #volume unit is Liters
    volume_liters = data_geom[0, 0, frame]
    #convert volume from liters (L) to cubic meters (m3)
    volume_m3 = volume_liters/1000 
    #calculate radius of the hypotetic sphere with this volume
    r = ((3 * volume_m3) / (4 * np.pi))**(1/3)
    print ('volume (m3)=', volume_m3)
    print ('sphere radius (m)=', r)

#HOW MUCH CHROMOSOME HAS BEEN BOUND BY RNAPOLY  
   #caclulates how much chromosome is explored by RNApolymerase
   #RNA poly footprint is 75nt and it binds ssDNA and the region is dsDNA
   #this function returns a list (a set actually) of nt that are touched by RNApoly up to a certain frame
   #It assumes that RNA poly moves always in direction forward 
   #TAKES A LOT to complete 
   
def rnapoly_exploration(frame, sim):
    x = h5py.File(sim)
    pos_data = x.get("states/RNAPolymerase/positionStrands/data")
    #pos_data.shape --> (156, 2, 29649)
    l1 =[]
    l2=[]
    #like saying for i in range 156: per ogni RNApoly
    for i in range(len(pos_data)):
        a=pos_data[i,0,:frame]
        b=pos_data[i,1,:frame]
        # a= array([ 17928, 266088,  84028, 196991,  82073])
        # b= array([2, 1, 2, 1, 2])
        idx1 = np.where(b == 1)
        idx2 = np.where(b == 2)
        #print(a[idx])
        for z in a[idx1]:
            if z!=0:
                #footprtnt +75 (strand 1 directionality 5>3)
                l1.append(z+range(75))
        for z in a[idx2]:
            if z!=0:
                #footprint -75 (strand 1 directionality 3>5)
                l2.append(z+range(75))
        #to array
        array1=np.array(l1).ravel()
        array2=np.array(l2).ravel()
        #to list
        all_positions1 = list(array1)
        all_positions2 = list(array2)
        #make a set to remove redundancy
        #break
    set_positions1 = set(all_positions1)
    set_positions2 = set(all_positions2)
    set_positions_sum = set_positions1.union(set_positions2)
        #break#
    percent_explored_strand1=len(set_positions1)*100/580076
    percent_explored_strand2=len(set_positions2)*100/580076
    percent_explored_sum=len(set_positions_sum)*100/580076
    print ('how many positions have been explored up to this frame:', len(set_positions_sum))
    print ('percentage of explored chromosome up to this frame (%):', percent_explored_sum )
    return set_positions_sum

#list_frames=[145,1189,6960] corresponding to 149, 1184, 6973

#usage
#rnapoly_expl_f149 = rnapoly_exploration(145, sim)
#rnapoly_expl_f1184 = rnapoly_exploration(1189, sim)
#rnapoly_expl_f6973 = rnapoly_exploration(6960, sim)

#make a dictionary ,then a json. 
def rnapoly_exploration_dic(dic, json_fname):
    #basically a list of indexes
    list_keys=[]
    for i in range(len(dic)):
        list_keys.append(i)
    list_key = set(list_keys)
    #the set coming from rnapoly_exploration has numpy.int32, that do not go on json, so we convert in float
    frames_float =[]
    for i in dic:
        frames_float.append(float(i))
    frames_float =  set(frames_float)
    dictionary_rnapoly_frames =dict(zip(list_keys, frames_float))
    #Dic2json(dictionary_rnapoly_frames, json_fname, output_folder)
    return dictionary_rnapoly_frames

#################### FUNCTIONS FOR RECIPE BUILDER 

#function to update the chain: if a specific pdb file has been assigned with a chain, only that PDB FILE with that specific name will update 
def updateChain(dic1, dic2):
    for i in dic1:
        f = dic1[i]['pdb_model']
        if i in dic2:
  #     update dictionary if the same model is found in dic2
            if f in dic2[i]['pdb_model']:
                dic1[i].update({'chain':dic2[i]['chain'], 'bu':dic2[i]['bu']})

def updateChainMembrane_new(dic1, dic2):
    for i in dic1:
        f= dic1[i]['pdb_model']
        if i in dic2:
            for el in range(len(dic2[i])):
                if f == dic2[i][el]['pdb_model']:
                    dic1[i].update({'chain':dic2[i][el]['chain'], 'bu':dic2[i][el]['bu'], 'offset':dic2[i][el]['offset'], 'pcpalAxis':dic2[i][el]['pcpalAxis']})
                    
#update the mapping dictionary with a template dictionary
def updateTemplates(dic1, dic2):
    for i in dic1:
        if i in dic2:
            dic1[i].update({'template':dic2[i]})
            
#use mod_quality dictionary to update all_mapping dictionary with the quality information. 
#dic1 will be updated with info in dic2 if they share a key 
def updateQuality (dic1, dic2, akey):
    for i in dic1:
        if i in dic2:
            dic1[i].update({'quality':dic2[i][akey]})
            
def updateQualityComplexes (dic1, dic2):
    for i in dic1:
        pdb = dic1[i]['pdb_model']
        if pdb in dic2:
            dic1[i].update({'quality':dic2[pdb]['quality']})

#NB IF A SPECIFIC PDB HAS BEEN POSITIONED AND THE POSITION IS REPORTED IN 'membrane' this script will overwrite the dictionary
#woks using the pdb_models, not using the protId            
def updateOffset(dic1, dic2):
    for i in dic1:
        f = dic1[i]['pdb_model']
        if i in dic2:
            if f in dic2[i]['pdb_model']:
                dic1[i].update({'offset':dic2[i]['offset'], 'pcpalAxis':dic2[i]['pcpalAxis']})

#update dictionary with dictionary
def updateDict(dic1, dic2):
    for i in dic1:
            dic1[i].update({'pdb_model':dic2[i]['pdb_model'], 'chain':dic2[i]['chain']})
            
#in a dictionary find all the keys that are empty, make a list and delete them            
def RemoveEmptyKeys(adic):
    l = []
    for i in adic:
        if len(adic[i])==0:
            l.append(i)
    for i in l:
        if i in adic:
            if len(adic[i])==0:
                del adic[i]
                
#write json file from a dictionary
def Dic2json(adic, fname):
    jsonf = fname+'.json'
    f= open(output_dir+jsonf,"w")
    f.write(json.dumps(adic))
    f.close()
    
#merge dictionaries
def MergeDictionaries(dic1, dic2, dic3):
    dic3 = dict(dic1)
    dic3.update(dic2)
    
#update dictionary general
def UpdateDictionary(dic1, dic2):
    for i in dic1:
        if i in dic2:
            dic1[i].update(dic2[i])

#make numpy array with a specific simulation 
def set_np_arrays(sim):  
    f = h5py.File(sim)
    #make the arrays for the protein complexes names
    labels0_compl = f.get("states/ProteinComplex/counts/labels/0")
    array_names_compl = np.array(labels0_compl)
    #make the arrays for the protein monomers names
    labels0_mono= f.get("states/ProteinMonomer/counts/labels/0")
    array_names_mono = np.array(labels0_mono)
    #data_mono shape is (frames, compartment, id )
    #data_compl shape is (id, compartment, frame)
    data_compl = f.get("states/ProteinComplex/counts/data")
    data_mono = f.get("states/ProteinMonomer/counts/data")
    return array_names_mono, array_names_compl, data_mono, data_compl

def ds_ss_DNA_binding(json_data): #list with dsDNA binding protein, dictionary with all DNA binding (ss + ds)
    covert = json.load(open(json_data, 'r'))
    l =[] #list only with dsDNA binding 
    d={} #dictionary with all DNA binding proteins 
    for i in covert['data']:
        if i['dna_footprint']!= None:
            d[i['wid']]={'dna_footprint':i['dna_footprint']['length'],'dna_binding':i['dna_footprint']['binding']}
            if i['dna_footprint']['binding']=='dsDNA':
                l.append(i['wid'])
    return l, d 

def edit_copy_numb_cluster_all_states(dic_copy_numb):    
    ##create a list with the process, folded etc proteins in varius states, sum their copy number to the mature ones and then remove them from the dictionary                        
    ##the only proteins left will be RNA bound 
    readme=open(output_dir+'READ_ME_PROTEIN_INFO.txt', 'w')
    to_remove =[]
    for i in copy_numb:
        if '-' in i:#this are the ingredients that have some other kind of states other than mature
            if d3[i.split('-')[0]]['binds_dsdna']=='0':#not DNA binding 
                #for all proteins not binding DNA and not ribosome or poly
                #inactive, prcessed, folded, bound etc
                if i !='RIBOSOME_70S-bound' and i !='RNA_POLYMERASE-bound' and i !='RNA_POLYMERASE_HOLOENZYME-bound':
                    #print (i)
                    to_remove.append(i)
                    #print (i, copy_numb[i], copy_numb[i.split('-')[0]])
                    new_copy_numb = copy_numb[i]['count']+copy_numb[i.split('-')[0]]['count']
                    #print ('sum this copy number to the mature one', i, copy_numb[i]['count'])
                    readme.write(i.split('-')[0]+' copy number is the sum of the mature state and '+i+'\n')
                    copy_numb[i.split('-')[0]].update({'count':new_copy_numb})
                else: #for ribosome and rna poly
                    to_remove.append(i)
                    #print ('this goes to lattice nucleoid', i, copy_numb[i]['count'])
                    readme.write(i+' copy number is accounted in the lattice nucleoid '+i.split('-')[0]+'\n')
            elif d3[i.split('-')[0]]['binds_dsdna']=='1': 
                if '-bound' in i:#remove DNA bound count - it is in the lattice
                    to_remove.append(i)
                    #print ('this goes to lattice nucleoid', i, copy_numb[i]['count'])
                    readme.write(i+' copy number is accounted in the lattice nucleoid '+i.split('-')[0]+'\n')
                else: 
                    #DNA binding protein not bound to the dna are considered cytoplasmic proteins when they are inactivated, processed etc
                    new_copy_numb = copy_numb[i]['count']+copy_numb[i.split('-')[0]]['count']
                    #print ('sum this copy number to the mature one', i, copy_numb[i]['count'])
                    readme.write(i.split('-')[0]+' copy number is the sum of the mature state and '+i+'\n')
                    copy_numb[i.split('-')[0]].update({'count':new_copy_numb})       
                    to_remove.append(i)
    readme.close()
    #remove DNA bound and proteins in other states than mature
    #the only protein left with the bound state are the bound to RNA, the structures of them will be the same as 
    for i in to_remove:
        if i in copy_numb:
            del copy_numb[i]   
    return copy_numb
  
#proteins count at one frame 
def all_state_proteins_at_one_frame(frame_m, frame_c, json_data):   
    all_data = json.load(open(json_data, 'r'))
    array_names_mono, array_names_compl, data_mono, data_compl, lprotein, lcomplex, l, adic = set_arrays(sim, json_data)
    #tot monomers = sum_over_compartments ( nascent + mature + bound + ... )
    tot_individual_proteins = np.sum(data_mono[frame_m,:,:],axis=0)
    #tot_complexes = sum_over_compartments ( nascent + mature + bound + ... )
    tot_complexes = np.sum(data_compl[:,:,frame_c],axis=1)
    anb = np.zeros(len(all_data['data']))
    for i in range(len(tot_individual_proteins)):
        prot_name = array_names_mono[i].decode().split("-")[0]
        state= array_names_mono[i].decode().split("-")[1]
        #if state=='mature' or state=='bound':
        #print (i, l.index(prot_name, tot_individual_proteins[i]))
        indice = l.index(prot_name)
        anb[indice] += tot_individual_proteins[i]    
    for i in range(len(tot_complexes)):
        complex_name = array_names_compl[i].decode().split("-")[0]
        state= array_names_compl[i].decode().split("-")[1]
        #complex_count = tot_complexes[i]
        #if state=='mature' or state=='bound':
        #print (i, l.index(prot_name, tot_individual_proteins[i]))
        indice = l.index(complex_name)
        anb[indice] += tot_complexes[i]
    adictionary={}
    for i in range(len(l)):
        adictionary[l[i]] = anb[i]
    return adictionary    

#write a recipe in csv format
def make_csv_recipe(dic,rna,fname):
    dictionary2csv = open (output_folder+fname+'.csv', 'w')
    dictionary2csv.write('name, molecular_weight, confidence (HHpred),pdb, selection, bu, uniprot, label, surface, compartment, comments, complexation, template, quality (modfold=monomers/voroMQA=complexes), offset, pcpalAxis, protein copy number, method, dsDNA-binding, function \n')
    list_homology = ['phyre', 'intfold', 'raptor', 'swiss', 'galaxy', 'itasser']
    dic==d3
    for aname in d3:
            quality='-1'
            #quality_modfold='-1'
            #quality_voroMQA='-1
            complexation='0.0' #default
            membrane = ''#false
            tree=fname+".mge"
            uniprot = ''
            mw =''
            offset = '0,0,0'
            axis = '0,0,1'
            if d3[aname]['offset']!='':
                offset = d3[aname]['offset'] 
                axis = d3[aname]['pcpalAxis']  
             #if the quality was assinged take it from d3, otherwise use the standard '-1'
             #for the monomers the quality is ModFOLD score
             #for the complexes te quality is voroMQA
            if d3[aname]['quality']!='':
                quality = d3[aname]['quality']
            if d3[aname]['mw']!=None:
                mw = d3[aname]['mw']
            #if the uniprot was assinged take it from d3, otherwise use the standard '-1'
            if 'uniprot' in d3[aname]:
                uniprot = d3[aname]['uniprot']     
            if d3[aname]['method']in list_homology:
                data_source='Homology Modeling'
            if d3[aname]['method']=='feig':
                data_source='CYT-MG-model'
            if d3[aname]['method']=='solved':
                data_source='PDB (MG experimentally solved)'
                #pdb=d3[aname]['pdb_model'].upper()
            if d3[aname]['method']=='PDB-homolog':
                data_source='PDB-homolog'
                #pdb=d3[aname]['pdb_model'].upper()
            if d3[aname]['method']=='PDB-homolog-edit':
                data_source='PDB-homolog (manually assembled)'
   #for monomers
            if d3[aname]['model']=='ProteinMonomer':
                #if d3[aname]['quality']!='':
                #    quality_modfold= d3[aname]['quality']
                # for 'real' monomers, aka monomers that do not participate in protein  complexes
                if len(d3[aname]['participation_in_complexes'])==0:
                    complexation='0.0' #monomers complexation = 0
                    #place them in the membrane
                    if d3[aname]['compartment']=='m' or d3[aname]['compartment']=='tm':
                         membrane='TRUE'
                    #place them in the extracell
                    if d3[aname]['compartment']=='e':
                        tree=fname
   # for monomers in complexes
            #place them in the membrane even when they are not complexed. 
                if len(d3[aname]['participation_in_complexes'])!=0:
                    complexation = '0.5'
                    if d3[aname]['compartment']=='m' or d3[aname]['compartment']=='tm':
                         membrane='TRUE'#Membrane protein complexation occurs following insertion [PUB_0018]
                         complexation = '0.5'
    # for complexes
            if d3[aname]['model']=='ProteinComplex':
             #   if d3[aname]['quality']!='':
             #       quality_voroMQA= d3[aname]['quality']
                complexation='1.0'
                #membrane placement
                if d3[aname]['compartment']=='m' or d3[aname]['compartment']=='tm':
                    membrane='TRUE'
                #extracell placement
                if d3[aname]['compartment']=='e':
                    tree=fname
            astr=aname+","
            astr+=str(mw)+','
            astr+=str(d3[aname]['HHpred'])+','
            astr+=str(d3[aname]['pdb_model'])+','
            astr+=str(d3[aname]['chain'])+','
            astr+=str(d3[aname]['bu'])+','
            astr+=str(uniprot)+','
            astr+=str(d3[aname]['function'])+','
            astr+=membrane+','
            astr+=tree+','
            astr+=str(d3[aname]['compartment'])+','
            astr+=complexation+','
            astr+=str(d3[aname]['template'])+','
            astr+=str(quality) +','
            #astr+=str(quality_modfold) +','
            #astr+=str(quality_voroMQA) +','
            astr+='"'+str(offset) +'",'
            astr+='"'+str(axis) +'",'
            astr+=str(d3[aname]['count'])+','
            astr+=str(data_source)+','
            astr+=str(d3[aname]['binds_dsdna'])+','
            astr+=str(d3[aname]['functional_category'])+'\n'
            dictionary2csv.write(astr)
    #after all the proteins add the RNA structures from rna dictionary
    for i in rna:
        dictionary2csv.write(str(i)+','+
                             ''+','+#mw
                             str('-1')+','+#confidence
                             str(rna[i]['pdb_model'])+','+
                             str(rna[i]['chain'])+','+
                             str(rna[i]['bu'])+','+
                             ''+','+ #uniprot
                             str(rna[i]['function'])+',,'+#surface
                             str('root.mge')+','+#compartment
                             str(rna[i]['compartment'])+#comments
                             ',,,'+
                             str('-1')+#quality
                             ',,,'+#
                             str(rna[i]['offset'])+','+
                             str(rna[i]['pcpalAxis'])+','+
                             str(rna[i]['count'])+','+
                             str(rna[i]['method'])+','+
                             str('0')+','+
                             str(rna[i]['function_category'])+'\n')
    dictionary2csv.close()










