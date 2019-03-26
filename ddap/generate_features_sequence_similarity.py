import os 
import glob
from collections import defaultdict
import pandas as pd
import numpy as np
from pandas.io.common import EmptyDataError
from ddap.utils import create_folder

import time

def generate_features_sequence_similarity(train_or_test,pair_dic,test_dic,output_path,number_top,strategy):
    # dic["sfdafd_H1-T1_HC:sfdafd_H1-T1_TN"] = [0,'atcgpoks','skdjugdodnbf']
    # 0 for "linkers that are not interacting with each other", 
    # 1 for "linkers that are interacting with each other")

    train_test_out_folder = os.path.join(output_path,train_or_test+'_files_for_BLAST')
    create_folder(train_test_out_folder)
    
    pair_df = preprocess_dic_to_df(pair_dic)
    test_df = preprocess_dic_to_df(test_dic)

    # run blast on head/tail sequences separately and generate a sequence similarity matrix
    # in blastp, search test_dic(query) against pair_dic(subject)
    
    head_query_df = test_df[['head_name', 'head_seq']].copy().sort_values('head_name')\
        .drop_duplicates().set_index('head_name', inplace=False)
    head_subject_df = pair_df[['head_name', 'head_seq']].copy().sort_values('head_name')\
        .drop_duplicates().set_index('head_name', inplace=False)
    head_sim_df = generate_similarity_matrix(head_query_df, head_subject_df, 'head',train_or_test,train_test_out_folder)

    tail_query_df = test_df[['tail_name', 'tail_seq']].copy().sort_values('tail_name')\
        .drop_duplicates().set_index('tail_name', inplace=False)
    tail_subject_df = pair_df[['tail_name', 'tail_seq']].copy().sort_values('tail_name')\
        .drop_duplicates().set_index('tail_name', inplace=False)
    tail_sim_df = generate_similarity_matrix(tail_query_df, tail_subject_df, 'tail',train_or_test,train_test_out_folder)
    
    known_pair_heads = head_subject_df.index.tolist()  
    
    score_dic = defaultdict(list)
    # record top n similarity features
    columns = ['head_score','tail_score','head_name','tail_name']
    start_time = time.time()
    
    for unknown_pair in test_dic.keys():  
        query_h, query_t = unknown_pair.split(":")
        score_df = pd.DataFrame(columns=columns,index = known_pair_heads)
        
        known_pair_heads_loop = [t for t in known_pair_heads if "_".join(t.split("_")[0:2]) != "_".join(query_h.split("_")[0:2])]
        # for training samples:
        # do not look for linkers in the same pathway
        for subject_h in known_pair_heads_loop: 
            subject_t = subject_h.replace('HC','TN')
            head_sim = head_sim_df.loc[query_h,subject_h]
            tail_sim = tail_sim_df.loc[query_t,subject_t]
            score_df.loc[subject_h] = [head_sim, tail_sim, subject_h, subject_t]
        
        score_df[['head_score', 'tail_score']] = score_df[['head_score', 'tail_score']].astype(float)
        score_df['mean_score'] = score_df[["head_score", "tail_score"]].mean(axis=1) 
        score_df['min_score'] = score_df[["head_score", "tail_score"]].min(axis=1) 

        score_df = score_df[['min_score','mean_score','head_score','tail_score','head_name','tail_name']]
        sort_col = "{}_score".format(strategy)
        top_scores = score_df.nlargest(number_top,[sort_col]).values.ravel()   
        score_dic[unknown_pair] = top_scores 
    
    end_time = time.time()    
    # print(end_time - start_time)  

    out_name = 'similarity_feature_for_{}_data.txt'.format(train_or_test)

    with open(os.path.join(train_test_out_folder,out_name),"w") as f:
        print("pair_ID","pair_label", "\t".join(columns*number_top) ,sep="\t",file=f)
        for unknown_pair in score_dic.keys():
            print(unknown_pair,test_dic[unknown_pair][0],"\t".join(list(str(x) for x in score_dic[unknown_pair])),sep="\t",file=f)
    return score_dic
    # score_dic["name of head:name of tail"] = [mean_similarity, head similarity,tail_similarity]

def preprocess_dic_to_df(dic):
    df = pd.DataFrame.from_dict(dic, orient='index')
    df.columns = ['label', 'head_seq','tail_seq']
    df['head_name'] = df.index.map(lambda s:s.split(':')[0])
    df['tail_name'] = df.index.map(lambda s:s.split(':')[1])
    return df

def generate_similarity_matrix(query_df,subject_df,h_or_t,train_or_test,directory):

    run_blastp(query_df,subject_df,h_or_t,directory)

    all_result_df = read_df(os.path.join(directory,h_or_t+"_results_blast.txt"))
    out_file = os.path.join(directory, "{}_{}_similarity_matrix.txt".format(train_or_test,h_or_t))
    sim_df = write_norm_matrix(all_result_df, query_df,subject_df,out_file )
    return sim_df

def run_blastp(query_df,subject_df, h_or_t, path):  
    subject_file =  os.path.join(path,h_or_t +"_blast_subject.fa")
    write_fasta_from_df(subject_df, subject_file)

    query_file =  os.path.join(path,h_or_t +"_blast_query.fa")
    write_fasta_from_df(query_df, query_file)
    blast_out = os.path.join(path,h_or_t + "_results_blast.txt")
    # run all v all BLAST
    cmd = "blastp -query {} -subject {} -out {} -outfmt 7".format(query_file, subject_file, blast_out)
    os.system(cmd)

def read_df(file):
    try:
        df = pd.read_csv(file, sep='\t', header=None,skip_blank_lines=True, comment='#', usecols=[0,1,11])
    except EmptyDataError:
        df = pd.DataFrame()
    return df 

def write_fasta_from_df(df,out_file):
    # This function takes in a dataframe 
    # index = sequence names;  'seq' = sequences
    # and write this sequence into .fasta files
    out_f = open(out_file,"w")
    for index, row in df.iterrows(): 
        print(">"+index,file=out_f)
        print(row.iloc[0],file=out_f)
    out_f.close()

def write_norm_matrix(df,query_df,subject_df,out_file):
    # This function does several things:
    # (1) change the dataframe from "long" format to "wide" format
    # (2) normalize the dataframe, and write this normalized df to output
    # (3) return a dictionary

    FILL = 0
    df.columns = ['query', 'subject', 'bit-score']
    query_list = df['query'].tolist()
    subject_list = df['subject'].tolist()

    full_query_list  = query_df.index.tolist()
    full_subject_list  = subject_df.index.tolist()

    # some subject sequences may not have similarity score with any 
    # query sequence. we need to use the following code to make sure 
    # the similarity matrix contains all query sequences and 
    # subject sequences.
    for name in full_query_list:
        if name not in query_list:
            df.loc[len(df)]=([name,subject_list[0],FILL])
    for name in full_subject_list:
        if name not in subject_list:
            df.loc[len(df)]=([query_list[0],name,FILL])

    df_pivot = pd.pivot_table(df,index=['query'],
        columns = 'subject', values = "bit-score", fill_value = FILL) # fill the empty bit_scores with 0.
    #df_pivot_normalized = df_pivot.apply(lambda x: x/x.max(), axis=1)
    df_pivot.to_csv(out_file, header=True,index=True)
    return df_pivot