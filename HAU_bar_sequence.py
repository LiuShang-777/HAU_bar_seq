#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 16:42:53 2019

@author: ls
"""

import pandas as pd
import os
import argparse
import warnings
warnings.filterwarnings('ignore')
parser=argparse.ArgumentParser(description='scripts for sequence extract for HAU barbadense genome')
parser.add_argument('-gn','--genename',type=str,help='the gene name file')
parser.add_argument('-gm','--genome',type=str,help='the genome file')
parser.add_argument('-gf','--gfffile',type=str,help='the gff file')
parser.add_argument('-sc','--selectclass',type=str,help='select class',choices=['gene','three_prime_UTR','five_prime_UTR',
                                                                              'exon','mRNA'])
parser.add_argument('-rf','--resultfile',type=str,help='result file')
parser.add_argument('-up','--upstream',type=int,help='the upstream')
args=parser.parse_args()
genome_file=args.genome
gff_file=args.gfffile
select_class=args.selectclass
result_file=args.resultfile
upstream=args.upstream
gene_name_file=args.genename

with open(gene_name_file,'r') as file:
    target_gene_name=[]
    for line in file:
        line=line.strip()
        target_gene_name.append(line)
with open (genome_file,'r') as file:
    name=[]
    seq=[]
    s=''
    seq_all=[]
    for line in file:
        line=line.strip()
        if line[0]=='>':
            name.append(line)
            s=''.join(seq)
            seq_all.append(s)
            s=''
            seq=[]
        else:
            line=line.strip()
            seq.append(line)
        
    s=''.join(seq)
    seq_all.append(s)
    s=''  
    seq=[]
    del seq_all[0]
print('genome file has been loaded')
if os.path.isdir('barbadense'):
    print('the barbadense directory has been exist, just use it')
else:    
    os.mkdir('barbadense')
    for i in range(26):
        with open('barbadense/%s.fa'%name[i][1::],'w') as file:
            file.write(name[i]+'\n')
            file.write(seq_all[i])
    print('the barbadense directory containing 26 chromosomes has been created')
def get_sequence(s,start,down,strand):
    dic=str.maketrans('ACTG','TGAC')
    if strand=='+':
        sequence=s[start-1:down]
    elif strand=='-':
        sequence=s[start-1:down]
        sequence=sequence[::-1]
        sequence=sequence.translate(dic)
    else:
        print('the format is wrong')
    return sequence
with open(gff_file,'r') as file:
	file_list=[]
	for line in file:
		line=line.strip()
		if line[0]=='#':
			pass
		else:
			file_list.append(line)
with open('barbadense_modified_gff.gff3','w') as file:
	for i in file_list:
		file.write(i+'\n')
	del file_list
gff_file='barbadense_modified_gff.gff3'
all_gff=pd.read_csv(gff_file,sep='\t',header=None)
print('the genome gff file has been imported')
new_8=[]
for i in all_gff[8]:
    new_8.append(i[3:18])
all_gff[8]=new_8
del new_8
chrome=[]
for i in all_gff[0]:
    if i[0]!='S':
        chrome.append(i)
clean_gff=all_gff[all_gff[0].isin(chrome)]
clean_gff=clean_gff[clean_gff[8].isin(target_gene_name)]
del chrome
gene_gff=clean_gff[clean_gff[2]==select_class]
if len(gene_gff)==0:
    print('genes in list have no select class you choose')
dic={}
for i in os.listdir('barbadense'):
    with open('barbadense/%s'%i,'r') as file:
        for line in file:
            if line[0]=='>':
                pass
            else:
                dic[i[0:-3]]=line
gene_name=[]
sequence=[]
c_sequence=[]
for i in gene_gff[8]:
    gene_name.append(i)
    temp=gene_gff[gene_gff[8]==i]
    s=dic[temp.iloc[0,0]]
    sequence.append(get_sequence(s,temp.iloc[0,3],temp.iloc[0,4],
                                 temp.iloc[0,6]))
    if temp.iloc[0,6]=='+':        
        c_sequence.append(get_sequence(s,temp.iloc[0,3]-upstream,temp.iloc[0,3],
                                 temp.iloc[0,6]))
    else:
        c_sequence.append(get_sequence(s,temp.iloc[0,4],temp.iloc[0,4]+upstream,
                                 temp.iloc[0,6]))
with open(result_file,'w') as file:
    for i,j in zip(gene_name,sequence):
        file.write('>'+i+'\n')
        file.write(j+'\n')
with open('upstream.txt','w') as file:
    for i,j in zip(gene_name,c_sequence):
        file.write('>'+i+'\n')
        file.write(j+'\n')
if select_class!='gene':   
        gene_list=[]
        for i in clean_gff[clean_gff[2]=='gene'][8]:
            if i not in gene_list:
                gene_list.append(i)
        utr_3_list=[]
        for i in gene_gff[8]:
            if i not in utr_3_list:        
                utr_3_list.append(i)
        no_3_utr_list=[]
        for i in gene_list:
            if i not in utr_3_list:
                no_3_utr_list.append(i)
        print('in target list %d has utr'%len(utr_3_list))
        print('in target list %d has no utr'%len(no_3_utr_list))
        print('in target list %d in it'%len(gene_list))
        utr_3_list=list(set(utr_3_list))
        with open('distribution_utr.txt','w') as file:
            file.write('utr: %d'%len(utr_3_list)+'\n')
            file.write('no_utr: %d'%len(no_3_utr_list)+'\n')
        uniq=[]
        ugene_name=[]
        repeat=[]
        for i in gene_gff[8]:
            ugene_name.append(i)
        for i in ugene_name:
            if i not in uniq:
                uniq.append(i)
            else:
                repeat.append(i)
        with open('utr_count.txt','w') as file:
            for i in uniq:
                file.write(i+'\t'+str(ugene_name.count(i))+'\n')
