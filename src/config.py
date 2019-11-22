
import os
home_path = os.path.abspath(os.path.join(os.getcwd(), ".."))
IEA_FLAG=False
GODIM=128  # the dimension of the feature vectors of GO terms
ONTOLOGY_FILE=home_path+"/datasets/gene_ontology_data/go_20180919_termdb.obo-xml"
ONTOLOGY_EDGELIST_FOLDER=home_path+"/term_encoding_output"

ECOLI_ASSOCIAION_FILE = home_path + '/datasets/gene_ontology_data/ecocyc.gaf'
HUMAN_ASSOCIAION_FILE = home_path + '/datasets/gene_ontology_data/goa_human_2018.gaf'
YEAST_ASSOCIAION_FILE = home_path + '/datasets/gene_ontology_data/yeast_2018.gaf'

# They are output files of the term encoding module,which save the feature vector representation of GO terms
BP_TERM_EMB_FILE_PATH = ONTOLOGY_EDGELIST_FOLDER + "/BP.emb"
CC_TERM_EMB_FILE_PATH = ONTOLOGY_EDGELIST_FOLDER + "/CC.emb"
MF_TERM_EMB_FILE_PATH = ONTOLOGY_EDGELIST_FOLDER + "/MF.emb"

not_found_go_terms=[]