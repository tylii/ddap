#!/usr/bin/env python
# -*- coding: utf-8 -*-

import warnings
warnings.filterwarnings("ignore")
import re
import random
import sklearn
from sklearn import ensemble
from sklearn import metrics
from sklearn.model_selection import KFold
from collections import defaultdict
import argparse
import os
from os.path import dirname, abspath 
import sys
import glob
import numpy as np 
import itertools
import pickle
from rdkit import Chem
from rdkit.Chem import Draw, AllChem,rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from ddap.generate_features_sequence_similarity import generate_features_sequence_similarity
from ddap.generate_sd_residue_feature import generate_sd_residue_feature_for_TEST  
from ddap.utils import create_folder,read_head_tail_fasta_to_test_dic,read_antiSMASH_gbk_to_test_dic,check_file_type,read_fa_to_test_dic,read_csv_to_test_dic,query_yes_no
from __init__ import get_data_path
import matplotlib.pyplot as plt
from cairosvg import svg2png
import subprocess
import shutil


def main():
    ap = argparse.ArgumentParser()
    # ------- Output File Arguments ----------
    ap.add_argument("-o","--output", type=str,required=True,
      help="output folder")
    # ------- Input File Arguments ----------
    ap.add_argument("--head", type=str,
      help="the .fasta file that contains Heads (C-terminal DD sequences)")
    ap.add_argument("--tail", type=str,
      help="the .fasta file that contains Tails (N-terminal DD sequences)")
    ap.add_argument('-i',"--input", type=str,
      help='tnput file for pathway prediction; accepted file formats are ".gbk", ".csv" and ".fasta"/".fa".')

    args = vars(ap.parse_args())

    PATH_INPUT= args['input'] 
    OUTPUT_PATH = args["output"]
    C_length = 100
    N_length = 50 
    DATA_PATH = get_data_path()  # abs path of data folder
    create_folder(OUTPUT_PATH)

    MODE = check_input(args)

    if MODE == 'dd':
        # ------ run the program ------------------
        print('[INFO] Reading DD sequences ...')
        test_dic = read_head_tail_fasta_to_test_dic(args["head"],args["tail"])
        print('[INFO] Predicting DD affinity ...')
        predict_dd_affinity(test_dic,OUTPUT_PATH)
       
    if MODE == 'pathway':
        # -------------- preprocess output of antismash --------------------
        print('[INFO] Reading antiSMASH output ...')
        input_type = check_file_type(PATH_INPUT)
        if input_type == 'antiSMASH':
            gene_dic = read_antiSMASH_gbk_to_test_dic(PATH_INPUT,C_length,N_length)
        elif input_type == 'fa':
            gene_dic = read_fa_to_test_dic(PATH_INPUT,C_length,N_length)
        elif input_type == 'csv':
            gene_dic = read_csv_to_test_dic(PATH_INPUT,C_length,N_length)
        else:
            print('[ERROR] invalid input file!')
            sys.exit()
        
        if len(gene_dic.keys())>=10:
            if not query_yes_no('[INFO] Input PKS has {} genes. This may result in large (>1G) output files. Do you want to proceed?'.format(len(gene_dic.keys()))):
                sys.exit()

        print('[INFO] Predicting DD affinity ...')
        test_dic_permu = get_test_dic_from_bgc_info(gene_dic)
        Ypred2,Ypred3,Ypred_stacked,test_ids = predict_dd_affinity(test_dic_permu,OUTPUT_PATH)
        print('[INFO] Predicting pathway scores ...')
        path_score_dic = predict_pathway(Ypred2,Ypred3,Ypred_stacked,test_ids,gene_dic,OUTPUT_PATH,DATA_PATH)
    
    print('[INFO] Prediction finished! Your result is in this folder:\n\t{}'.format(OUTPUT_PATH))

        
def check_input(args):
    if args['head'] and args['tail'] and (not args['input']):
        if os.path.isfile(args['head']) and os.path.isfile(args['tail']):
            return 'dd'
        else:
            print('[ERROR] Input file(s) does not exist.')
            sys.exit()
    elif args['input'] and (not args['head']) and (not args['tail']):
        if os.path.isfile(args['input']):
            return 'pathway'
        else:
            print('[ERROR] Input file does not exist.')
            sys.exit()   
    else:
        print('[ERROR] Please specify input file with --head HEAD.fa --tail TAIL.fa or --input your_folder') 
        sys.exit()

def predict_dd_affinity(test_dic,OUTPUT_PATH):
    TOP_N = 3
    STRATEGY = 'mean'
    WEIGHT3 = 0.6
    THREADS = 1
    DATA_PATH = get_data_path()

    OUTPUT_PATH = os.path.join(OUTPUT_PATH,'Docking_Domain_Affinity_Prediction_Result')
    create_folder(OUTPUT_PATH)

    train_dic = pickle.load(open(os.path.join(DATA_PATH,'train_dic.pkl'), 'rb')) 
    # ------------ prepare train and test dataset ------------
    cv_out_path = OUTPUT_PATH
    test_ids = list(test_dic.keys()) 
    Ypred_list = []
    # ---------- Generate Feature 2 (Model1) For Testing Data ------------- 
    Xtest2 = get_seq_sim_feature(train_dic,test_dic,cv_out_path, TOP_N ,STRATEGY)
    model_path = os.path.join(DATA_PATH,'model_1.pkl')
    Ypred2 = make_prediction(model_path,Xtest2)
    score_out_file = os.path.join(OUTPUT_PATH,'{}_DD_affinity_prediction.csv'.format('Model1'))
    write_dd_predicted_scores(score_out_file,test_ids,Ypred2)
    Ypred_list.append(Ypred2)
    # Ypred2 = None
    # ---------- Generate Feature 3 (Model2) For Testing Data ------------- 
    all_SD_pairs = pickle.load(open(os.path.join(DATA_PATH,'all_SD_pairs.pkl'), 'rb'))
    train_head_ali = os.path.join(DATA_PATH,'heads.afa')
    train_tail_ali = os.path.join(DATA_PATH,'tails.afa')
    Xtest3 = get_sd_residue_feature(test_dic,all_SD_pairs,cv_out_path,THREADS,train_head_ali,train_tail_ali)
    Ypred3 = make_prediction(os.path.join(DATA_PATH,'model_2.pkl'),Xtest3)
    score_out_file = os.path.join(OUTPUT_PATH,'{}_DD_affinity_prediction.csv'.format('Model2'))
    write_dd_predicted_scores(score_out_file,test_ids,Ypred3)
    Ypred_list.append(Ypred3)
    # Ypred3 = None
    # -------------  stack predictions made by different features ---------------------
    Ypred_stacked = stack_predictions(Ypred_list,WEIGHT3)
    score_out_file = os.path.join(OUTPUT_PATH,'{}_DD_affinity_prediction.csv'.format('DDAP'))
    write_dd_predicted_scores(score_out_file,test_ids,Ypred_stacked)

    return Ypred2,Ypred3,Ypred_stacked,test_ids

def predict_pathway(Ypred2,Ypred3,Ypred_stacked,test_ids,gene_dic,OUTPUT_PATH,DATA_PATH):
    '''Output the predicted score for each possible pathway
    '''
    out_folder = os.path.join(OUTPUT_PATH,'Pathway_Prediction_Result')
    create_folder(out_folder)
    out_f = open(os.path.join(out_folder,'pathway_prediction.csv'),'w')
    DELIM = ','
    header = DELIM.join(['PathwayID','Pathway','Score','Score2','Score3','ProductSMILES','Polymer'])
    print(header,file=out_f)
    
    pred_dic = defaultdict(lambda:defaultdict())
    for i in range(len(test_ids)):
        pred_dic[test_ids[i]]['pred'] = Ypred_stacked[i]
        pred_dic[test_ids[i]]['pred2'] = Ypred2[i]
        pred_dic[test_ids[i]]['pred3'] = Ypred3[i]

    path_score_dic = defaultdict(int)
    
    gene_list = list(gene_dic.keys())
    pathway = list(range(len(gene_list)))   # pathway = [0,1,2 ....]
    for full_path in itertools.permutations(pathway):
        path_score,path_score2,path_score3 = get_path_score(full_path,gene_list,pred_dic)
        path = ';'.join([gene_list[x] for x in full_path]) 
        smiles_str,polymer_str = get_SMILES(path,gene_dic,DATA_PATH)
        path_score_dic[path] = [path_score,path_score2,path_score3,smiles_str,polymer_str]

    sorted_by_path_score = sorted(path_score_dic.keys(), key=lambda x: path_score_dic[x][0],reverse = True)
    path_id = 1
    for path in sorted_by_path_score:
        row = [path_id,path]+path_score_dic[path]
        print(DELIM.join([str(x) for x in row]),file=out_f )
        path_id+=1
    out_f.close()

    # --------- plot compound structure ----------
    max_n = 10 # plot structure for top n 
    top_n = min(max_n,len(sorted_by_path_score))

    plot_folder = os.path.join(out_folder,'structure_images')
    create_folder(plot_folder)
    for i,path in enumerate(sorted_by_path_score[0:top_n]):
        smiles_str = path_score_dic[path][3]
        path_id = i+1
        plot_smiles(path_id,smiles_str,plot_folder)
    return path_score_dic

def linker_pair_prediction_performance(Ytest, Ypred):
    '''This function defines the DD affinity prediction scoring metric
    '''
    score = metrics.roc_auc_score(Ytest, Ypred)
    return score 

def pathway_prediction_performance(Ypred_permu,bgc_dic,test_ids_permu,out_file):
    
    out_f = open(out_file,"w")
    print("BGC_ID","pathway","score",sep="\t",file=out_f)
    pred_dic = defaultdict(float)
    for a_id, pred in zip(test_ids_permu, Ypred_permu):
        pred_dic[a_id] = float(pred)

    bgc_score_dic = defaultdict(int)
    for bgc in bgc_dic.keys():
        rank_of_true_path = 1
        n_dd = int(bgc_dic[bgc])
        if n_dd < 4 or n_dd > 10:
            continue
        true_path = list(range(1,n_dd+1))
        true_path_score = get_path_score(true_path,bgc,pred_dic)
        path =  ';'.join([str(x) for x in true_path]) 
        print(bgc,path,true_path_score,sep='\t',file=out_f )

        mid_dd = list(range(2,n_dd))
        for permu in itertools.permutations(mid_dd):
            full_path = [1]+ list(permu) +[n_dd]
            if ';'.join([str(x) for x in full_path]) == ';'.join([str(x) for x in true_path]):
                continue
            path_score = get_path_score(full_path,bgc,pred_dic)
            if path_score > true_path_score:
                rank_of_true_path += 1       
            path = ';'.join([str(x) for x in full_path]) 
            print(bgc,path,path_score,sep='\t',file=out_f )
        bgc_score_dic[bgc] = int(rank_of_true_path)
    out_f.close()
    return bgc_score_dic
    # bgc_score_dic['BGC_2']['1234']

def get_path_score(full_path,gene_list,pred_dic):
    '''calculate score for a given pathway
    INPUT: 
    full_path = [1,4,3,2]
    '''

    score = []
    score2 = []
    score3 = []
    for i in range(len(full_path)-1): # iterate over all DD pairs 
        head = "{}_H".format(gene_list[full_path[i]])
        tail = "{}_T".format(gene_list[full_path[i+1]])
        pair_ID = head+":"+tail
        score.append(pred_dic[pair_ID]['pred'])
        score2.append(pred_dic[pair_ID]['pred2'])
        score3.append(pred_dic[pair_ID]['pred3'])
    return np.mean(score),np.mean(score2),np.mean(score3)


def get_seq_sim_feature(train_dic,test_dic,cv_out_path,TOP_N ,STRATEGY):
    '''get sequence similarity features
    '''
    test_ids = list(test_dic.keys())
    Xtest = [[] for x in test_ids]
    test_data_feature2 = generate_features_sequence_similarity('test',train_dic,
            test_dic,cv_out_path,TOP_N ,STRATEGY)
    for i,one_name in enumerate(test_ids,0):
        tmp=[]
        for jj in range(TOP_N):
            tmp += list(test_data_feature2[one_name][6*jj:6*jj+4])
        Xtest[i] = tmp
    return Xtest

def get_sd_residue_feature(test_dic,all_SD_pairs,cv_out_path,THREADS,train_head_ali,train_tail_ali):
    '''get cognate residue features
    '''
    test_ids = list(test_dic.keys())
    Xtest = [[] for x in test_ids]
    test_data_feature3 = \
        generate_sd_residue_feature_for_TEST(test_dic,cv_out_path,all_SD_pairs,train_head_ali,train_tail_ali,THREADS)

    for i,one_name in enumerate(test_ids,0):
        Xtest[i] = test_data_feature3[one_name]
    return Xtest

def print_log(string,log_f):
    print(string)
    print(string,file=log_f)

def write_dd_pair_prediction_result(test_ids,test_ids_permu,Ytest_permu,Ypred_permu,file_name):
    with open(file_name,"w") as out_f:
        print("pair_ID","true_label","predicted_label",sep="\t",file=out_f)
        for i, one_name in enumerate(test_ids_permu):
            if one_name in test_ids:
                print(one_name,Ytest_permu[i],Ypred_permu[i],sep="\t",file=out_f)
            else:
                print(one_name,"+0",Ypred_permu[i],sep="\t",file=out_f)

def write_dd_predicted_scores(out_file,test_ids,Ypred):
    '''Write the predicted dd affinity scores to a file
    '''
    sep=','
    with open(out_file,"w") as out_f:
        print('Head(C-terminal-DD)', 'Tail(N-terminal-DD)','PredictedScore', 'PredictedLabel',sep=sep,file=out_f)
        for i, one_name in enumerate(test_ids):
            head,tail=one_name.split(':')
            lab = 1 if Ypred[i] > 0.5 else 0
            print(head,tail,Ypred[i],lab,sep=sep,file=out_f)


def write_ranking_result_for_all_pathway(file_name,pathway_score_dic):
    out_f_path=open(file_name,"w")
    print("BGC_ID","rank_of_true_pathway",sep="\t",file=out_f_path)
    for k,v in pathway_score_dic.items():
        print(k,v,sep="\t",file=out_f_path)
    out_f_path.close()

def make_prediction(model_path,Xtest):
    '''The base learner
    '''
    est = pickle.load(open(model_path, 'rb'))
    Ypred = list(est.predict(Xtest))
    return Ypred

def stack_predictions(Ypred_permu_list, weight):
    weight=float(weight)
    if len(Ypred_permu_list) == 1:
        return Ypred_permu_list[0]
    else:
        Ypred_stacked = np.array(Ypred_permu_list[0])*(1-weight)+np.array(Ypred_permu_list[1])*weight
        return Ypred_stacked.tolist()

# ----------
def get_positive_negative_samples(head_seq_dic,tail_seq_dic):
    '''get all positive and negative samples from sequences dictionary
    '''
    # list of all negative and positive samples
    all_negative = []
    all_positive = []
    # ---------  get all positive and negative samples from the linker sequence dictionary ---------
    for name_head, name_tail in itertools.product(head_seq_dic.keys(),tail_seq_dic.keys()):
        # an example of name_head is: "BGC_18_LINKER_7_HC", "HC" means Head, C-terminal
        # an example of name_head is: "BGC_18_LINKER_7_TN", "TN" means Tail, N-terminal
        name_h = re.sub('_HC','',name_head)
        name_h_pathway = re.sub('_LINKER.*','',name_h)
    
        name_t = re.sub('_TN','',name_tail)
        name_t_pathway= re.sub('_LINKER.*','',name_t)
        if name_h == name_t:
            all_positive.append(name_h+"_HC:"+name_t+"_TN")
        elif name_h_pathway == name_t_pathway:
            all_negative.append(name_h+"_HC:"+name_t+"_TN")

    # ---------- random sample the same number of negative samples as positive samples -----------
    random.seed(123)
    all_negative = random.sample(all_negative, len(all_positive))

    return all_positive,all_negative

    
def get_test_dic_from_bgc_info(bgc_dic):
    test_dic_permu = defaultdict()
    all_heads = defaultdict()
    all_tails = defaultdict()
    for gene in bgc_dic.keys():
        all_heads[gene+'_H'] = bgc_dic[gene]['H']
        all_tails[gene+'_T'] = bgc_dic[gene]['T']
    for h,t in itertools.product(all_heads.keys(),all_tails.keys()):
        test_dic_permu['{}:{}'.format(h,t)] = [0,all_heads[h],all_tails[t]]
    return test_dic_permu

def load_smiles(file):
    """Load smiles from a dictionary mapping residues to SMILES string
    From antiSMASH source code
    """
    aa_smiles = {}
    smiles_monomer = open(file,'r')

    for line in smiles_monomer.readlines():
        line = line.strip()
        if not line or line.startswith('#') or line == "END":
            continue
        smiles = line.split()
        assert len(smiles) == 2, "Invalid smiles line {!r}".format(line)

        aa_smiles[smiles[0]] = smiles[1]

    smiles_monomer.close()
    return aa_smiles

def get_SMILES(path,gene_dic,DATA_PATH):
    smiles_file = os.path.join(DATA_PATH,'aaSMILES.txt')
    smiles_map = load_smiles(smiles_file)
    smiles_str =''
    genes_in_path = path.split(';')
    polymer_str = ''
    for gene in genes_in_path:
        polymer_str+='('
        if 'mono' not in gene_dic[gene].keys():
            continue
        mono_code = gene_dic[gene]['mono']
        polymer_str+='-'.join(mono_code)
        for mono in mono_code:
            if mono in smiles_map.keys():
                smiles_str += smiles_map[mono]
                
            else:
                print('unknown monomer:{}'.format(mono))
        polymer_str+=')'
    return smiles_str,polymer_str

def plot_smiles(path_id,smiles_str,out_folder):

    # size = (150,150) # size of the figure
    m = Chem.MolFromSmiles(smiles_str)
    try:
        fig_path = os.path.join(out_folder,'pathway_{}.png'.format(path_id))
    except:
        print(smiles_str)
    svg_code = moltosvg(m)
    svg2png(bytestring=svg_code,write_to=fig_path)


def moltosvg(mol,molSize=(450,450),kekulize=True):
    mc = Chem.Mol(mol.ToBinary())
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg.replace('svg:','')


if __name__ == "__main__":
    main()
