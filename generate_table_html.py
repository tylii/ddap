import re
import csv
import glob
import json 
import shutil
import os

def main():
    # write the fixed stuff
    name1 = '/Users/ran/Documents/UM/BGCstructurePrediction/restructured_project/dev/ddap/database_page/tmp.txt'      # input file
    name2 = '/Users/ran/Documents/UM/BGCstructurePrediction/restructured_project/dev/ddap/database_page/html/table1.html'  # output file
    
    fin = open(name1, "r")
    data2 = fin.read()
    fin.close()
    fout = open(name2, "w")
    fout.write(data2)
    print('',file=fout)

    #----- write column names ---------
    columns = ['BGC ID','BGC source','Product','Reference','Pathway Source','Gene Order in Pathway']
    for name in columns:
        new_name = re.sub(' ','&nbsp;',name)
        print('<th>'+new_name+'</th>',file = fout)
    
    #---- print fixed content ----
    print('</tr>\n</thead>\n<tbody>',file = fout)
    
    #----- print published BGC ------
    info_f = open('/Users/ran/Documents/UM/BGCstructurePrediction/restructured_project/data/raw_data/collect_t1pks/BGC_list/combined_bgc_list_full_information_REVISED.csv','r')
    next(info_f)
    list_of_products = []
    csv_info = csv.reader(info_f)
    count_bgc = 0 
    for row in csv_info:
        count_bgc +=1
        print("<tr data-status='published'>",file=fout)
        print("<td>{}</td>".format(row[0]),file=fout)               # BGC_ID

        if re.search('json$',row[1]):
            print("<td><a href='https://mibig.secondarymetabolites.org/repository/{}/'>{}</a></td>".format(row[1].replace('.json',''),row[1].replace('.json','')),file=fout)       # MIBIG ID
        else:
            print("<td>{}</td>".format(row[1]),file=fout) 
        print("<td>{}</td>".format(row[2].lower()),file=fout)       # product name
        list_of_products.append(row[2].lower())
        print("<td><a href='https://www.ncbi.nlm.nih.gov/pubmed/{}'>{}</a></td>".format(row[4],row[4]),file=fout)       # PMID
        print("<td><span class='label label-success'>Published</span></td>",file=fout)      # publish or predicted
        pathway = row[8]
        pathway = re.sub(';','&#8594',pathway)
        print("<td>{}</td>".format(pathway),file=fout)              # pathway
        print("</tr>",file=fout)
    info_f.close()
    #----- print predicted BGC ------
    result_files = glob.glob('/Users/ran/Documents/UM/BGCstructurePrediction/restructured_project/model_results/predict_pathway_for_PKS_in_MIBIG/ddap_prediction_collection/*')
    info_f = open('/Users/ran/Documents/UM/BGCstructurePrediction/restructured_project/data/raw_data/collect_t1pks/BGC_list/combined_bgc_list_full_information_REVISED.csv','r')
    
    
    for file in result_files:
        if os.stat(file).st_size == 0:
            continue
        mibig_id = shutil._basename(file).split('.')[0]
        json_name = mibig_id+'.json'
        product_name,pmid = get_bgc_info(json_name)

        # ---- remove replicated products -----------
        if product_name.lower() in list_of_products:
            print(product_name)
            continue
        elif product_name.lower().split(' ')[0] in [x.split(" ")[0] for x in list_of_products]:
            print(product_name)
            continue
        count_bgc+=1
        # ------ copy this file to desired dir -----------
        copy_to_dir = '/Users/ran/Documents/UM/BGCstructurePrediction/restructured_project/model_results/predict_pathway_for_PKS_in_MIBIG/ddap_prediction_collection_for_download'
        max_line = 5e4
        bgc_id = 'BGC_{}'.format(count_bgc)
        copy_file(file,copy_to_dir,max_line,bgc_id)
        # ----- actually print the pathway information ---------
        
        n_pathways = 10  # only print top 10 pathways
        result_f = open(file,'r')
        next(result_f)
        result_csv = csv.reader(result_f)

        print("<tr data-status='predicted'>",file=fout)
        print("<td>BGC_{}</td>".format(count_bgc),file=fout)           # BGC_ID
        print("<td><a href='https://mibig.secondarymetabolites.org/repository/{}/'>{}</a></td>".format(mibig_id,mibig_id),file=fout)       # MIBIG ID
        print("<td>{}</td>".format(product_name.lower()),file=fout)       # product name  
        print_pmid_list(pmid,fout)      # PMID
        print("<td><span class='label label-danger'>Predicted</span></td>",file=fout)      # publish or predicted

        count_path=1
        for row in result_csv:
            if count_path==1:
                pathway1 = row[1]
                pathway1 = re.sub(';','&#8594',pathway1)
                print('<td class="text-left">{}&nbsp;&nbsp;&nbsp;<button class="btn btn-sm btn-light" id="show_{}">+</button></td>'.format(pathway1,count_bgc),file=fout)              # pathway
                print("</tr>",file=fout)

                print('<tr data-status="predicted">\n<td colspan="4">\n<div id="extra_{}" style="display: none;">'.format(count_bgc),file=fout) #style="display: none;""
                print('<table class="sub-table table-hover">\n<thead class="thead-light">\n',file=fout)
                print('<tr>\n<th>#</th>\n<th>Pathway</th>\n<th>Score</th>\n<th>SMILES&nbsp;String</th>\n</tr>\n</thead>\n<tbody>'.format(count_bgc),file=fout)
                
            if count_path<=n_pathways:
                print('<tr>',file=fout)
                print("<td>{}</td>".format(count_path),file=fout)
                pathway = row[1]
                pathway = re.sub(';','&#8594',pathway)
                print("<td>{}</td>".format(pathway),file=fout)
                print("<td>{}</td>".format(row[2]),file=fout)
                print("<td>{}</td>".format(row[5]),file=fout)
                count_path+=1
            else:
                break

            print('</tr>',file=fout)
        
        print('</tbody>\n</table>\n</div>\n</td>\n</tr>',file=fout)
    
    print('</tbody>\n</table>\n</div>',file=fout)

    fout.close()

def get_bgc_info(json_file_name):
    f = open('/Users/ran/Documents/UM/BGCstructurePrediction/restructured_project/data/raw_data/MiBIG/mibig_json_1.4_only_t1pks/{}'.format(json_file_name),'r')
    data = json.load(f)
    prodcut_name = data['general_params']['compounds'][0]['compound']
    try:
        pmid = data['general_params']['publications']
        if not (type(pmid) is str):
            pmid=pmid[0]
    except:
        pmid = ""
    return prodcut_name, pmid 

def print_pmid_list(pmid,fout):
    print("<td>",file=fout)
    ids = []
    for a_id in pmid.replace(' ','').split(','):
        if re.search('doi',a_id):
            ids.append("<a href='https://{}'>{}</a>".format(a_id,a_id)) 
        else:
            ids.append("<a href='https://www.ncbi.nlm.nih.gov/pubmed/{}'>{}</a>".format(a_id,a_id)) 
    print(', '.join(ids),file=fout)
    print("</td>",file=fout)

def copy_file(file,copy_to_dir,max_line,bcg_id):
    '''copy file to new folder and change name
    '''
    in_file = open(file, 'r')
    out_file = open(os.path.join(copy_to_dir,bcg_id+"_pathway_prediction.csv"),'w')
    for i in range(1,int(max_line)):
        line=in_file.readline()
        print(line,end='',file=out_file)
        
if __name__ == '__main__':
    main()