
###########step 1: 根据注释文件，得到每个蛋白质（基因产物)的Term 集合#####################

import os
from src.helper import  *
import src.config as config


home_path = os.path.abspath(os.path.join(os.getcwd(), ".."))

def load_test_data(dataset_name,ontology_name,split_mode,iea_name):
    test_data_file=home_path+"/intermediate_outputs/"+split_mode+"/"+dataset_name+"/"+iea_name+"/"+ontology_name+"/test.npy"
    return load_data(test_data_file)

def load_train_data(dataset_name,ontology_name,split_mode,iea_name):
    train_data_file=home_path+"/intermediate_outputs/"+split_mode+"/"+dataset_name+"/"+iea_name+"/"+ontology_name+"/train.npy"
    return load_data(train_data_file)

def load_cross_data(dataset_name,ontology_name,split_mode,iea_name):
    cross_data_file=home_path+"/intermediate_outputs/"+split_mode+"/"+dataset_name+"/"+iea_name+"/"+ontology_name+"/cross.npy"
    return load_data(cross_data_file)

def load_test_data_for_fivefold(dataset_name,ontology_name,split_mode,iea_name,part_name):
    test_data_file=home_path+"/intermediate_outputs/"+split_mode+"/"+dataset_name+"/"+iea_name+"/"+ontology_name+"/test."+part_name+".npy"
    return load_data(test_data_file)

def load_train_data_for_fivefold(dataset_name,ontology_name,split_mode,iea_name,part_name):
    train_data_file=home_path+"/intermediate_outputs/"+split_mode+"/"+dataset_name+"/"+iea_name+"/"+ontology_name+"/train."+part_name+".npy"
    return load_data(train_data_file)

def load_cross_data_for_fivefold(dataset_name,ontology_name,split_mode,iea_name,part_name):
    cross_data_file=home_path+"/intermediate_outputs/"+split_mode+"/"+dataset_name+"/"+iea_name+"/"+ontology_name+"/cross."+part_name+".npy"
    return load_data(cross_data_file)

# if __name__ == '__main__':
#
#     # # step1: 获得注释数据
#     home_path = os.path.abspath(os.path.join(os.getcwd(), ".."))
#     dataset_names=["yeast","ecoli","human"]
#     for dataset_name in dataset_names:
#         if config.IEA_FLAG==True:
#             iea_name="iea+"
#         else:
#             iea_name="iea-"
#         if  dataset_name=="ecoli":
#             ASSOCIAION_FILE=config.ECOLI_ASSOCIAION_FILE
#         if dataset_name=="yeast":
#             ASSOCIAION_FILE=config.YEAST_ASSOCIAION_FILE
#         if dataset_name=="human":
#             ASSOCIAION_FILE=config.HUMAN_ASSOCIAION_FILE
#         geneAssociation = GeneAssociation()
#         geneAssociation.parse_association_file(ASSOCIAION_FILE,config.IEA_FLAG)
#         print("BP CC MF 注释数目")
#         print(len(geneAssociation.BPassociations))
#         print(len(geneAssociation.CCassociations))
#         print(len(geneAssociation.MFassociations))
#
#         src_folder = home_path + "/datasets/processed_ppi/connectivity/" + dataset_name + "/" + iea_name
#         proteinlistBP, proteinlistCC, proteinlistMF = get_all_proteins(src_folder)
#
#         proteinlistset = set(proteinlistBP + proteinlistCC + proteinlistMF)
#         print(dataset_name)
#         print(len(proteinlistset))
#
#         # 解析基因本体文件
#         ontology=Ontology()
#         ontology.parse_from_xml_file(config.ONTOLOGY_FILE)
#         # construct termVector
#         BPTermVectors=TermVectors()
#         BPTermVectors.parse_term_embedding_file(config.BP_TERM_EMB_FILE_PATH)
#         BPTermVectors.construct_term_vector_dict()
#         print(len(BPTermVectors.termVectorDict))
#
#         CCTermVectors = TermVectors()
#         CCTermVectors.parse_term_embedding_file(config.CC_TERM_EMB_FILE_PATH)
#         CCTermVectors.construct_term_vector_dict()
#         print(len(CCTermVectors.termVectorDict))
#
#         MFTermVectors = TermVectors()
#         MFTermVectors.parse_term_embedding_file(config.MF_TERM_EMB_FILE_PATH)
#         MFTermVectors.construct_term_vector_dict()
#         print(len(MFTermVectors.termVectorDict))
#
#         # BP_protein_vector_dict=present_protein_with_go_term_vector(proteinlistBP,geneAssociation.BPassociations,BPTermVectors,ontology.BP_terms,ontology.AbandonedBPGOTerms)
#         # CC_protein_vector_dict = present_protein_with_go_term_vector(proteinlistCC, geneAssociation.CCassociations, CCTermVectors,ontology.CC_terms,ontology.AbandonedCCGOTerms)
#         # MF_protein_vector_dict = present_protein_with_go_term_vector(proteinlistMF, geneAssociation.MFassociations, MFTermVectors,ontology.MF_terms,ontology.AbandonedMFGOTerms)
#
#         BP_protein_vector_dict = present_protein_with_go_term_vector(proteinlistBP, geneAssociation.BPassociations,
#                                                                      BPTermVectors, ontology.BP_terms,
#                                                                      )
#         CC_protein_vector_dict = present_protein_with_go_term_vector(proteinlistCC, geneAssociation.CCassociations,
#                                                                      CCTermVectors, ontology.CC_terms,
#                                                                      )
#         MF_protein_vector_dict = present_protein_with_go_term_vector(proteinlistMF, geneAssociation.MFassociations,
#                                                                      MFTermVectors, ontology.MF_terms,
#                                                                      )
#
#         #in oder to save runtime, save train_data cross_data, test_data to intermediate_outputs where protien is represented by term vecor
#         # src_folder=home_path+"/datasets/processed_ppi/fivefold/"+dataset_name+"/"+iea_name
#         # out_folder=home_path+"/intermediate_outputs/fivefold/"+dataset_name+"/"+iea_name
#         # save_all_protein_pairs_as_term_vector_for_fivefold(src_folder, out_folder, BP_protein_vector_dict,CC_protein_vector_dict,MF_protein_vector_dict)
#
#         src_folder = home_path + "/datasets/processed_ppi/connectivity/" + dataset_name + "/" + iea_name
#         out_folder = home_path + "/intermediate_outputs/connectivity/" + dataset_name + "/" + iea_name
#         save_all_protein_pairs_as_term_vector(src_folder, out_folder, BP_protein_vector_dict,
#                                                            CC_protein_vector_dict, MF_protein_vector_dict)
#
#         # load data
#
#
#
#
#
#
#
#



#
