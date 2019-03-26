# -*- coding: utf-8 -*-
import os
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
import glob
import sys
import re
import csv
import shutil

def create_folder(output_path):
    if os.path.exists(output_path):
        if query_yes_no('[INFO] This folder already exists. Do you want to replace it?\n\t{}'.format(output_path)):
            os.system("rm -r "+ output_path)
        else:
            sys.exit()
    os.makedirs(output_path)

def read_fasta_to_dic(file):
    '''read fasta file into a dictionary
    '''
    with open(file, "r") as head_f:
        seq_dic = defaultdict(str)
        for line in head_f:
            line = line.rstrip()
            if line.startswith('>'):
                seq_name = line[1:]
            else:
                seq_dic[seq_name] += line
    return seq_dic

def read_fasta_to_dic_replace_old(file):
    '''read fasta file into a dictionary
    replace the old sequence if a name is seen multiple times
    '''
    with open(file, "r") as head_f:
        seq_dic = defaultdict(str)
        for line in head_f:
            line = line.rstrip()
            if line.startswith('>'):
                seq_name = line[1:]
                if seq_name in seq_dic.keys():
                    seq_dic[seq_name]=''
            else:
                seq_dic[seq_name] += line
    return seq_dic


def read_fasta_to_dic_with_names(file):
    names = []
    with open(file, "r") as head_f:
        seq_dic = defaultdict(str)
        for line in head_f:
            line = line.rstrip()
            if line.startswith('>'):
                seq_name = line[1:]
                names.append(seq_name)
            else:
                seq_dic[seq_name] += line
    return seq_dic,names

def read_head_tail_fasta_to_test_dic(head_file,tail_file):
    test_dic = defaultdict(lambda: [0,'',''])

    head_dic,heads = read_fasta_to_dic_with_names(head_file)
    tail_dic,tails = read_fasta_to_dic_with_names(tail_file)

    if len(heads)!=len(tails):
        print('[ERROR] The number of Heads should be the same as the number of Tails.')
        sys.exit()
    for x_h,x_t in zip(heads,tails):
        name = "{}:{}".format(x_h,x_t)
        test_dic[name][1] = head_dic[x_h]
        test_dic[name][2] = tail_dic[x_t]
    return test_dic


def read_fasta_to_df(file):

    with open(file, "r") as head_f:
        seq_dic = defaultdict(list)
        for line in head_f:
            line = line.rstrip()
            if line.startswith('>'):
                seq_name = line[1:]
            else:
                seq_dic[seq_name] += list(line)
    seq_df = pd.DataFrame.from_dict(seq_dic, orient='index')
    return seq_df

def read_bgc_to_test_dic(folder):
    '''Read the antiSMASH output
    INPUT: antismash folder
    OUTPUT: gene_dic
    gene_dic[seq_id][seq] = protein_sequence
    gene_dic[seq_id]['H'] = Head DD sequence (C-terminal of protein)
    gene_dic[seq_id]['T'] = Tail DD sequence (N-terminal of protein)
    gene_dic[seq_id]['mono'] = monomer
    '''
    C_length = 100
    N_length = 50 

    if glob.glob(os.path.join(folder,'*.region001.gbk')):
        file = glob.glob(os.path.join(folder,'*.region001.gbk'))[0]
        gene_dic = read_antismash5(file,C_length,N_length)
    else:
        try:
            file = glob.glob(os.path.join(folder,'*.cluster001.gbk'))[0]
        except:
            print('[ERROR] Input file not valid.')
        gene_dic = read_antismash4(file,C_length,N_length)

    return gene_dic

def check_file_type(file):
    name = shutil._basename(file)
    with open(file) as f:
        first_line = f.readline().rstrip()
        second_line = f.readline().rstrip()

    if re.search('gbk$|gb$',name,re.IGNORECASE):
        return 'antiSMASH'
    if re.search(',',second_line) and re.search('csv$',name,re.IGNORECASE):
        return 'csv'
    if first_line.startswith('>') and re.search('fa$|fasta$',name,re.IGNORECASE):
        return 'fa'
    
    print('[ERROR] Invalid file format. Accepted file formats are ".gbk", ".csv" and ".fasta"/".fa".')
    sys.exit()


def read_fa_to_test_dic(file,C_length,N_length):
    '''Read the fa pathway input
    INPUT: fasta file
    OUTPUT: gene_dic
    gene_dic[seq_id][seq] = protein_sequence
    gene_dic[seq_id]['H'] = Head DD sequence (C-terminal of protein)
    gene_dic[seq_id]['T'] = Tail DD sequence (N-terminal of protein)
    gene_dic[seq_id]['mono'] = monomer
    '''

    seq_dic = read_fasta_to_dic(file)
    gene_dic = defaultdict(lambda: defaultdict())

    for seq_name in seq_dic.keys():
        seq = seq_dic[seq_name]
        head,tail = get_N_C_seq(seq,C_length,N_length)
        gene_dic[seq_name]['seq'] = seq
        gene_dic[seq_name]['H'] = head
        gene_dic[seq_name]['T'] = tail
    return gene_dic

def read_csv_to_test_dic(file,C_length,N_length):
    '''Read the csv pathway input
    INPUT: csv file
    OUTPUT: gene_dic
    gene_dic[seq_id][seq] = protein_sequence
    gene_dic[seq_id]['H'] = Head DD sequence (C-terminal of protein)
    gene_dic[seq_id]['T'] = Tail DD sequence (N-terminal of protein)
    gene_dic[seq_id]['mono'] = monomer
    '''
    f = open(file,'r')
    reader = csv.reader(f)
    gene_dic = defaultdict(lambda: defaultdict())
    for row in reader:
        if len(row)==1:
            print('[ERROR] input csv file must have at least 2 columns (column1: sequence name; column2: sequence)')
            sys.exit()
        seq_name = row[0]
        seq = row[1]
        gene_dic[seq_name]['seq'] = seq
        head,tail = get_N_C_seq(seq,C_length,N_length)
        gene_dic[seq_name]['H'] = head
        gene_dic[seq_name]['T'] = tail
        if len(row)==3:
            gene_dic[seq_name]['mono'] = row[2]
    
    return gene_dic


def read_antiSMASH_gbk_to_test_dic(file,C_length,N_length):
    '''Read the antiSMASH output
    INPUT: antismash gbk file
    OUTPUT: gene_dic
    gene_dic[seq_id][seq] = protein_sequence
    gene_dic[seq_id]['H'] = Head DD sequence (C-terminal of protein)
    gene_dic[seq_id]['T'] = Tail DD sequence (N-terminal of protein)
    gene_dic[seq_id]['mono'] = monomer
    '''

    if re.search('region',file):
        gene_dic = read_antismash_region(file,C_length,N_length)
    elif re.search('cluster',file):
        gene_dic = read_antismash_cluster(file,C_length,N_length)
    else:
        print('[ERROR] Input file not valid.')
        sys.exit()
    if not gene_dic:
        print('[ERROR] Input .gbk is not a PKS region.')
        sys.exit()
    if len(gene_dic.keys())==1:
        print('[ERROR] Input .gbk contains only one PKS gene. Need at least 2 PKS genes for pathway prediction.')
        sys.exit()
    return gene_dic
    

def read_antismash_cluster(file,C_length,N_length):
    seq_record = list(SeqIO.parse(file, "genbank"))[0]
    # genes_in_gbk = [x for x in seq_record.features if x.type == 'CDS']
    gene_dic = defaultdict(lambda: defaultdict())

    for seq in [x for x in seq_record.features]: 
        if seq.type == 'CDS' and \
            'sec_met' in seq.qualifiers.keys() and \
                any([re.search('Type: t1pks|PKS_',x) for x in seq.qualifiers['sec_met']]):
                # seq.qualifiers['protein_id'][0] in gene_dic.keys():
                # 'gene_functions' in seq.qualifiers.keys() and \
                # re.search('biosynthetic',seq.qualifiers['gene_functions'][0]):

            if  'gene' in seq.qualifiers:
                seq_id = seq.qualifiers['gene'][0]
            elif 'locus_tag' in seq.qualifiers:
                seq_id = seq.qualifiers['locus_tag'][0] 
            else:
                seq_id = seq.qualifiers['protein_id'][0] 
            protein_seq = seq.qualifiers['translation'][0] 
            gene_dic[seq_id]['seq'] = protein_seq
            # print(seq_id,gene_dic[seq_id]['seq'])

            head,tail = get_N_C_seq(protein_seq,C_length,N_length)
            gene_dic[seq_id]['H'] = head
            gene_dic[seq_id]['T'] = tail
            try:
                gene_dic[seq_id]['mono'] = seq.qualifiers['aSProdPred'][0].split('-')
            except:
                pass

    return gene_dic
    
def get_N_C_seq(protein_seq,C_length,N_length):
    '''Get the C-terminal DD and N-termianl DD given the protein sequences
    '''
    C_seq = min(C_length,len(protein_seq))
    N_seq = min(N_length,len(protein_seq))
    head = protein_seq[-C_seq:]
    tail = protein_seq[:N_seq]
    return head,tail


def read_antismash_region(file,C_length,N_length):
    seq_record = list(SeqIO.parse(file, "genbank"))[0]
    # genes_in_gbk = [x for x in seq_record.features if x.type == 'CDS']
    gene_dic = defaultdict(lambda: defaultdict())

    for seq in [x for x in seq_record.features]: 
        if seq.type == 'CDS' and \
            'gene_functions' in seq.qualifiers.keys() and \
                any([re.search('t1pks|mod_KS|PKS_',x) for x in seq.qualifiers['gene_functions']]):
                # seq.qualifiers['protein_id'][0] in gene_dic.keys():
                # 'gene_functions' in seq.qualifiers.keys() and \
                # re.search('biosynthetic',seq.qualifiers['gene_functions'][0]):

            if  'gene' in seq.qualifiers:
                seq_id = seq.qualifiers['gene'][0]
            elif 'locus_tag' in seq.qualifiers:
                seq_id = seq.qualifiers['locus_tag'][0] 
            else:
                seq_id = seq.qualifiers['protein_id'][0] 
            protein_seq = seq.qualifiers['translation'][0] 
            gene_dic[seq_id]['seq'] = protein_seq
            # print(seq_id,gene_dic[seq_id]['seq'])

            C_seq = min(C_length,len(protein_seq))
            N_seq = min(N_length,len(protein_seq))
            gene_dic[seq_id]['H'] = protein_seq[-C_seq:]
            gene_dic[seq_id]['T'] = protein_seq[:N_seq]

    for seq in [x for x in seq_record.features]: 
        if seq.type == 'aSDomain' and \
            seq.qualifiers['locus_tag'][0] in gene_dic.keys() and\
            'specificity' in seq.qualifiers.keys() and \
            any([re.search('consensus:',x) for x in seq.qualifiers['specificity']]):

            seq_id = seq.qualifiers['locus_tag'][0]

            for x in  seq.qualifiers['specificity']:
                if re.search('consensus:',x):
                    mono = x.replace('consensus: ','')
                    break

            if mono and 'mono' in gene_dic[seq_id].keys():
                gene_dic[seq_id]['mono'] += [mono]
            elif mono:
                gene_dic[seq_id]['mono'] = [mono]

            # print(seq_id,mono)

    return gene_dic

def query_yes_no(question):
    """Ask a yes/no question and return their answer.
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    prompt = " [y/n] "
    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()

        if choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")