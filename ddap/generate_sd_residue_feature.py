from collections import defaultdict
import re
import itertools
import pandas as pd
import numpy as np 
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import os
from ddap.utils import create_folder, read_fasta_to_dic,read_fasta_to_dic_replace_old


def generate_sd_residue_feature_for_TEST(test_dic_permu,cv_out_path,all_SD_pairs,train_head_ali,train_tail_ali,THREADS):
    SD_FOLDER = os.path.join(cv_out_path,'test_files_for_SD_residues')
    create_folder(SD_FOLDER)

    head_fa,tail_fa = \
        write_head_tail_fasta(test_dic_permu,SD_FOLDER)

    test_head_ali, test_tail_ali = \
        align_with_clustalo_with_profile(head_fa,tail_fa,train_head_ali,train_tail_ali,THREADS)

    all_SD_pairs_in_test = \
        match_residue_position_and_extract_feature(train_head_ali,train_tail_ali,test_head_ali, test_tail_ali,all_SD_pairs)

    test_head_ali_dic = read_fasta_to_dic_replace_old(test_head_ali)   
    test_head_num_dic = digitalize_aa(test_head_ali_dic)
    test_tail_ali_dic = read_fasta_to_dic_replace_old(test_tail_ali)  
    test_tail_num_dic = digitalize_aa(test_tail_ali_dic)
        # using read_fasta_to_dic() will cause an error when 
    # the testing sequences have the same name as the training sequences

    test_data_feautre = \
        extract_residue_by_position(test_dic_permu,all_SD_pairs_in_test,test_head_num_dic,test_tail_num_dic)
    write_features_to_file(test_data_feautre,test_dic_permu,SD_FOLDER)
    return test_data_feautre

def write_head_tail_fasta(train_test_dic,SD_FOLDER):
    head_seq_dic =  defaultdict(str)
    tail_seq_dic = defaultdict(str)

    for k, v in train_test_dic.items():
        head,tail = k.split(':')
        head_seq_dic[head] = v[1]
        tail_seq_dic[tail] = v[2]
    
    head_fa = os.path.join(SD_FOLDER,'heads.fa')
    write_fasta_from_dic(head_seq_dic,head_fa)

    tail_fa = os.path.join(SD_FOLDER,'tails.fa')
    write_fasta_from_dic(tail_seq_dic,tail_fa)

    return head_fa,tail_fa

def align_with_clustalo(head_fa,tail_fa,THREADS):
    head_ali = head_fa.replace('.fa','.afa')
    tail_ali = tail_fa.replace('.fa','.afa')

    # cmd = "clustalo -i {} -o {} --outfmt=a2m --force --threads={}".format(head_fa, head_ali,THREADS)
    cmd = "clustalo -i {} -o {} --outfmt=a2m --force ".format(head_fa, head_ali)
    os.system(cmd)

    cmd = "clustalo -i {} -o {} --outfmt=a2m --force ".format(tail_fa, tail_ali)
    os.system(cmd)
    return head_ali, tail_ali

def align_with_clustalo_with_profile(head_fa,tail_fa,train_head_ali,train_tail_ali,THREADS):
    head_ali = head_fa.replace('.fa','.afa')
    tail_ali = tail_fa.replace('.fa','.afa')
    cmd = "clustalo -i {} --p1={} -o {} --outfmt=a2m --force ".format(head_fa,train_head_ali,head_ali)
    # cmd = "clustalo -i {} --p1={} -o {} --outfmt=a2m --force --threads={}".format(head_fa,train_head_ali,head_ali,THREADS)
    os.system(cmd)
    
    cmd = "clustalo -i {} --p1={} -o {} --outfmt=a2m --force ".format(tail_fa,train_tail_ali,tail_ali)
    # cmd = "clustalo -i {} --p1={} -o {} --outfmt=a2m --force --threads={}".format(tail_fa,train_tail_ali,tail_ali,THREADS)
    os.system(cmd)
    return head_ali, tail_ali

def write_fasta_from_dic(seq_dic, out_file):
    out_f = open(out_file,"w")
    for i in seq_dic.keys():    
        print(">{}".format(i),file=out_f)
        print(seq_dic[i],file=out_f)
    out_f.close()

def digitalize_aa(seq_dic):
    AA_DIC = generate_AA_map_dic()
    seq_dic_num = defaultdict(list)
    for seq_name, seq in seq_dic.items():
        aa_list = list(seq)
        aa_list=[int(AA_DIC[aa]) for aa in aa_list]
        seq_dic_num[seq_name]=aa_list
        #seq_df = pd.DataFrame.from_dict(seq_dic_num)
    return seq_dic_num

def k_clustering(data_array,N_CLUSTER,SEED):
    # Number of clusters
    kmeans = KMeans(n_clusters = N_CLUSTER, random_state = SEED)
    # Fitting the input data
    kmeans = kmeans.fit(data_array)
    # Getting the cluster labels
    labels = kmeans.predict(data_array)
    # Centroid values
    centroids = kmeans.cluster_centers_

    return labels,centroids

def pca_visualization(array,labels,N_CLUSTER,SD_FOLDER):
    pca = PCA(n_components=3)
    projected = pca.fit_transform(array)
    fig = plt.figure(1, figsize=(9.5,9.5))
    ax = fig.add_subplot(2,2,1)
    #plt.subplot(221)
    ax.scatter(projected[:, 0], projected[:, 1],
            c=labels, edgecolor='none', alpha=0.8,
            cmap=plt.cm.get_cmap('Accent', N_CLUSTER))
    ax.set_xlabel('PC 1')
    ax.set_ylabel('PC 2')
    #fig.colorbar()
    
    ax = fig.add_subplot(2,2,2)
    ax.scatter(projected[:, 0], projected[:, 2],
            c=labels, edgecolor='none', alpha=0.8,
            cmap=plt.cm.get_cmap('Accent', N_CLUSTER))
    ax.set_xlabel('PC 1')
    ax.set_ylabel('PC 3')

    ax = fig.add_subplot(2,2,3)
    ax.scatter(projected[:, 1], projected[:, 2],
            c=labels, edgecolor='none', alpha=0.8,
            cmap=plt.cm.get_cmap('Accent', N_CLUSTER))
    ax.set_xlabel('PC 2')
    ax.set_ylabel('PC 3')
    #fig.colorbar()
    #plt.show()

    ax = fig.add_subplot(224, projection='3d')
    ax.scatter(projected[:, 0], projected[:, 1], projected[:, 2],
            c=labels, edgecolor='none', alpha=0.8,
            cmap=plt.cm.get_cmap('Accent', N_CLUSTER))
    ax.set_xlabel('PC 1')
    ax.set_ylabel('PC 2')
    ax.set_zlabel('PC 3')
    #fig.colorbar()

    fig.savefig(os.path.join(SD_FOLDER,'docking_domain_{}_clusters.pdf'.format(N_CLUSTER)))

def generate_AA_map_dic():
    AA_DIC = {}
    HYDROPHOBIC_AA = ['F', 'M', 'W', 'I', 'V', 'L', 'P', 'A']  
    HYDROPHILIC_AA = ['N', 'C','Q', 'G', 'S', 'T', 'Y']        
    CHARGED_P_AA = ['K','R','H']                               
    CHARGED_N_AA = ['D','E']                                    

    for aa in HYDROPHOBIC_AA:
        AA_DIC[aa] = 0 
    for aa in HYDROPHILIC_AA:
        AA_DIC[aa] = 0 
    for aa in CHARGED_P_AA:
        AA_DIC[aa] = 1
    for aa in CHARGED_N_AA:
        AA_DIC[aa] = -1
    AA_DIC['X'] = 0
    AA_DIC['-'] = 0
    return AA_DIC

def extract_residue_by_position(train_dic,all_SD_pairs,head_dic_num,tail_dic_num):
    data_feautre = defaultdict(list)
    for k in train_dic.keys():
        head,tail = k.split(':')
        for pair in all_SD_pairs:
            head_i, tail_i = pair.split(';')
            add_this = sorted([int(head_dic_num[head][int(head_i)]) , int(tail_dic_num[tail][int(tail_i)])])
            data_feautre[k] += add_this
    return data_feautre

def match_residue_position_and_extract_feature(train_head_ali,train_tail_ali,test_head_ali,test_tail_ali,all_SD_pairs):
    all_SD_pairs_in_test = []

    train_head_list = read_names_from_fa(train_head_ali)
    train_tail_list = read_names_from_fa(train_tail_ali)

    test_head_ali_dic = read_fasta_to_dic_replace_old(test_head_ali)
    test_tail_ali_dic = read_fasta_to_dic_replace_old(test_tail_ali)

    test_head_seq_len = len(list(test_head_ali_dic.values())[0])
    test_tail_seq_len = len(list(test_tail_ali_dic.values())[0])

    loc = -1
    map_head_loc = defaultdict(int)
    map_tail_loc = defaultdict(int)
    for i in range(test_head_seq_len):
        this_column = []
        for name in train_head_list:
            this_column.append(test_head_ali_dic[name][i])
        if not all(x == '-' for x in this_column):
            loc+=1
            map_head_loc[int(loc)] = int(i)
    
    loc = -1
    for i in range(test_tail_seq_len):
        this_column = []
        for name in train_tail_list:
            this_column.append(test_tail_ali_dic[name][i])
        if not all(x == '-' for x in this_column):
            loc+=1
            map_tail_loc[int(loc)] = int(i)
    
    for p in all_SD_pairs:
        head_loc, tail_loc = p.split(';')
        head_new_loc = map_head_loc[int(head_loc)]
        tail_new_loc = map_tail_loc[int(tail_loc)]
        all_SD_pairs_in_test.append('{};{}'.format(head_new_loc,tail_new_loc))

    return all_SD_pairs_in_test

def read_names_from_fa(file):
    seq_names = []
    with open(file, "r") as head_f:
        for line in head_f:
            line = line.rstrip()
            if line.startswith('>'):
                seq_names.append(line[1:])
    return seq_names
def write_features_to_file(data_feautre,train_test_dic,SD_FOLDER):
    out_f = open(os.path.join(SD_FOLDER,'train_sd_residue_feature.txt'),'w')
    print('linker_id','linker_label','features',sep='\t',file= out_f)

    for k,v in train_test_dic.items():
        print(k,v[0],'\t'.join([str(x) for x in data_feautre[k]]),sep='\t',file= out_f)
    out_f.close()

def write_seq_profile(head_list,seq_profile,SD_FOLDER):
    out_f = open(os.path.join(SD_FOLDER,'seq_profile_for_clustering.txt'),'w')
    for head, profile in zip(head_list,seq_profile):
        tail = head.replace('_HC','_TN')
        print("{}:{}".format(head,tail),'\t'.join([str(x) for x in profile]),sep='\t',file = out_f)
    out_f.close()

def write_cluster_result(head_list,labels,SD_FOLDER):
    out_f = open(os.path.join(SD_FOLDER,'clustering_results.txt'),'w')
    for head, label in zip(head_list,labels):
        tail = head.replace('_HC','_TN')
        print("{}:{}".format(head,tail),label,sep='\t',file = out_f)
    out_f.close()