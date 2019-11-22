


###########step 1: 根据注释文件，得到每个蛋白质（基因产物)的Term 集合#####################

import os
from src.helper import  *
import src.config as config

home_path = os.path.abspath(os.path.join(os.getcwd(), ".."))
 # load expression data

def load_expression_train_data(ontology_name,iea_name):
    train_data_file=home_path+"/intermediate_outputs/expression_data/train/"+iea_name+"/"+ontology_name+"/.npy"
    return  load_data(train_data_file)

def load_expression_test_data(ontology_name,iea_name):
    test_data_file = home_path + "/intermediate_outputs/expression_data/test/" + iea_name + "/" + ontology_name + "/.npy"
    return load_data(test_data_file)

if __name__ == '__main__':

    # # step1: 获得注释数据
    home_path = os.path.abspath(os.path.join(os.getcwd(), ".."))
    ASSOCIAION_FILE = config.YEAST_ASSOCIAION_FILE
    if config.IEA_FLAG == True:
        iea_name = "iea+"
    else:
        iea_name = "iea-"
    # 基因表达模型的训练集是整个Yeast数据集，测试集是Expression 数据集
    # 为Yeast 数据集生成intermediate_ouputs
    geneAssociation = GeneAssociation()
    geneAssociation.parse_association_file(ASSOCIAION_FILE,config.IEA_FLAG)

    ############### create intermediate_outputs for train expression data #######################
    # src_folder = home_path + "/datasets/raw_ppi/yeast/" + iea_name
    # proteinlistBP, proteinlistCC, proteinlistMF = get_all_proteins_for_expression(src_folder)
    # proteinlistset = set(proteinlistBP + proteinlistCC + proteinlistMF)
    # print("length of proteins in yeast datasets")
    # print(len(proteinlistset))
    # # 解析基因本体文件
    # ontology = Ontology()
    # ontology.parse_from_xml_file(config.ONTOLOGY_FILE)
    # # construct termVector
    # BPTermVectors = TermVectors()
    # BPTermVectors.parse_term_embedding_file(config.BP_TERM_EMB_FILE_PATH)
    # BPTermVectors.construct_term_vector_dict()
    # print(len(BPTermVectors.termVectorDict))
    #
    # CCTermVectors = TermVectors()
    # CCTermVectors.parse_term_embedding_file(config.CC_TERM_EMB_FILE_PATH)
    # CCTermVectors.construct_term_vector_dict()
    # print(len(CCTermVectors.termVectorDict))
    #
    # MFTermVectors = TermVectors()
    # MFTermVectors.parse_term_embedding_file(config.MF_TERM_EMB_FILE_PATH)
    # MFTermVectors.construct_term_vector_dict()
    # print(len(MFTermVectors.termVectorDict))
    #
    # # BP_protein_vector_dict=present_protein_with_go_term_vector(proteinlistBP,geneAssociation.BPassociations,BPTermVectors,ontology.BP_terms,ontology.AbandonedBPGOTerms)
    # # CC_protein_vector_dict = present_protein_with_go_term_vector(proteinlistCC, geneAssociation.CCassociations, CCTermVectors,ontology.CC_terms,ontology.AbandonedCCGOTerms)
    # # MF_protein_vector_dict = present_protein_with_go_term_vector(proteinlistMF, geneAssociation.MFassociations, MFTermVectors,ontology.MF_terms,ontology.AbandonedMFGOTerms)
    #
    # BP_protein_vector_dict = present_protein_with_go_term_vector(proteinlistBP, geneAssociation.BPassociations,
    #                                                              BPTermVectors, ontology.BP_terms,
    #                                                              )
    # CC_protein_vector_dict = present_protein_with_go_term_vector(proteinlistCC, geneAssociation.CCassociations,
    #                                                              CCTermVectors, ontology.CC_terms,
    #                                                              )
    # MF_protein_vector_dict = present_protein_with_go_term_vector(proteinlistMF, geneAssociation.MFassociations,
    #                                                              MFTermVectors, ontology.MF_terms,
    #                                                              )
    #
    # out_folder = home_path + "/intermediate_outputs/expression_data/train/"+ iea_name
    # save_as_term_vector_for_expression_train(src_folder, out_folder, BP_protein_vector_dict,
    #                                       CC_protein_vector_dict, MF_protein_vector_dict)

#

#     # ############### create intermediate_outputs for train expression data #######################
    src_folder = home_path + "/datasets/expression_data/"
    proteinlistBP, proteinlistCC, proteinlistMF = get_all_proteins_for_expression_test(src_folder)
    proteinlistset = set(proteinlistBP + proteinlistCC + proteinlistMF)
    print("length of proteins in yeast datasets")
    print(len(proteinlistset))
    # 解析基因本体文件
    ontology = Ontology()
    ontology.parse_from_xml_file(config.ONTOLOGY_FILE)
    # construct termVector
    BPTermVectors = TermVectors()
    BPTermVectors.parse_term_embedding_file(config.BP_TERM_EMB_FILE_PATH)
    BPTermVectors.construct_term_vector_dict()
    print(len(BPTermVectors.termVectorDict))

    CCTermVectors = TermVectors()
    CCTermVectors.parse_term_embedding_file(config.CC_TERM_EMB_FILE_PATH)
    CCTermVectors.construct_term_vector_dict()
    print(len(CCTermVectors.termVectorDict))

    MFTermVectors = TermVectors()
    MFTermVectors.parse_term_embedding_file(config.MF_TERM_EMB_FILE_PATH)
    MFTermVectors.construct_term_vector_dict()
    print(len(MFTermVectors.termVectorDict))

    # BP_protein_vector_dict=present_protein_with_go_term_vector(proteinlistBP,geneAssociation.BPassociations,BPTermVectors,ontology.BP_terms,ontology.AbandonedBPGOTerms)
    # CC_protein_vector_dict = present_protein_with_go_term_vector(proteinlistCC, geneAssociation.CCassociations, CCTermVectors,ontology.CC_terms,ontology.AbandonedCCGOTerms)
    # MF_protein_vector_dict = present_protein_with_go_term_vector(proteinlistMF, geneAssociation.MFassociations, MFTermVectors,ontology.MF_terms,ontology.AbandonedMFGOTerms)

    BP_protein_vector_dict = present_protein_with_go_term_vector(proteinlistBP, geneAssociation.BPassociations,
                                                                 BPTermVectors, ontology.BP_terms,
                                                                 )
    CC_protein_vector_dict = present_protein_with_go_term_vector(proteinlistCC, geneAssociation.CCassociations,
                                                                 CCTermVectors, ontology.CC_terms,
                                                                 )
    MF_protein_vector_dict = present_protein_with_go_term_vector(proteinlistMF, geneAssociation.MFassociations,
                                                                 MFTermVectors, ontology.MF_terms,
                                                                 )

    out_folder = home_path + "/intermediate_outputs/expression_data/test/"+ iea_name
    save_as_term_vector_for_expression_test(src_folder, out_folder, BP_protein_vector_dict,
                                          CC_protein_vector_dict, MF_protein_vector_dict)

