"""
The protein encoding module consists of two LSTM networks,
each encoding a protein into a vector (called protein vector)
via the feature vectors of the terms annotated with the protein.

The function of this file is to represent a protein as a sequence
of feature vectors of terms annotated to it, which is the input for
the protein encoding module
"""
from src.helper import  *

def generate_for_ppi():
    # Deal for PPI
    home_path = os.path.abspath(os.path.join(os.getcwd(), ".."))
    dataset_names = ["human"] #"yeast", "ecoli",
    split_modes = ["connectivity"]
    for dataset_name in dataset_names:
        for split_mode in split_modes:
            if config.IEA_FLAG == True:
                iea_name = "iea+"
            else:
                iea_name = "iea-"
            if dataset_name == "ecoli":
                ASSOCIAION_FILE = config.ECOLI_ASSOCIAION_FILE
            if dataset_name == "yeast":
                ASSOCIAION_FILE = config.YEAST_ASSOCIAION_FILE
            if dataset_name == "human":
                ASSOCIAION_FILE = config.HUMAN_ASSOCIAION_FILE

            # 解析基因本体文件
            ontology = Ontology()
            ontology.parse_from_xml_file(config.ONTOLOGY_FILE)
            geneAssociation = GeneAssociation()
            geneAssociation.parse_association_file(ASSOCIAION_FILE, config.IEA_FLAG,ontology,dataset_name)
            print("The number of annotation for BP CC MF")
            print(len(geneAssociation.BPassociations))
            print(len(geneAssociation.CCassociations))
            print(len(geneAssociation.MFassociations))

            src_folder = home_path + "/datasets/processed_ppi/" + split_mode + "/" + dataset_name + "/" + iea_name
            out_folder = home_path + "/intermediate_outputs/" + split_mode + "/" + dataset_name + "/" + iea_name
            print(src_folder)
            proteinlistBP, proteinlistCC, proteinlistMF = get_all_proteins(src_folder)

            proteinlistset = set(proteinlistBP + proteinlistCC + proteinlistMF)
            print(dataset_name)
            print(len(proteinlistset))


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
            # deal with PPI protein
            # in oder to save runtime, save train_data cross_data, test_data to intermediate_outputs
            # where protien is represented by vecor of GO term
            if split_mode == "connectivity":
                save_all_protein_pairs_as_term_vector(src_folder, out_folder, BP_protein_vector_dict,
                                                      CC_protein_vector_dict, MF_protein_vector_dict)
            if split_mode == "fivefold":
                save_all_protein_pairs_as_term_vector_for_fivefold(src_folder, out_folder, BP_protein_vector_dict,
                                                               CC_protein_vector_dict, MF_protein_vector_dict)

def generate_for_expression():
    # Deal for Expression
    home_path = os.path.abspath(os.path.join(os.getcwd(), ".."))
    dataset_name = "yeast"
    if config.IEA_FLAG == True:
        iea_name = "iea+"
    else:
        iea_name = "iea-"
    ASSOCIAION_FILE = config.YEAST_ASSOCIAION_FILE

    # 解析基因本体文件
    ontology = Ontology()
    ontology.parse_from_xml_file(config.ONTOLOGY_FILE)

    geneAssociation = GeneAssociation()
    geneAssociation.parse_association_file(ASSOCIAION_FILE, config.IEA_FLAG,ontology)

    # generate for train
    train_src_folder = home_path + "/datasets/raw_ppi/yeast/"+iea_name
    train_out_folder = home_path + "/intermediate_outputs/expression_data/train/" + iea_name

    proteinlistBP, proteinlistCC, proteinlistMF = get_all_proteins_for_expression(train_src_folder)
    proteinlistset = set(proteinlistBP + proteinlistCC + proteinlistMF)
    print(dataset_name)
    print(len(proteinlistset))


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

    BP_protein_vector_dict = present_protein_with_go_term_vector(proteinlistBP, geneAssociation.BPassociations,
                                                                 BPTermVectors, ontology.BP_terms,
                                                                 )
    CC_protein_vector_dict = present_protein_with_go_term_vector(proteinlistCC, geneAssociation.CCassociations,
                                                                 CCTermVectors, ontology.CC_terms,
                                                                 )
    MF_protein_vector_dict = present_protein_with_go_term_vector(proteinlistMF, geneAssociation.MFassociations,
                                                                 MFTermVectors, ontology.MF_terms)

    save_as_term_vector_for_expression_train(train_src_folder, train_out_folder, BP_protein_vector_dict, CC_protein_vector_dict, MF_protein_vector_dict)

    # generate for test
    test_src_folder = home_path + "/datasets/expression_data/"
    test_out_folder = home_path + "/intermediate_outputs/expression_data/test/" + iea_name
    proteinlistBP, proteinlistCC, proteinlistMF = get_all_proteins_for_expression_test(test_src_folder)
    proteinlistset = set(proteinlistBP + proteinlistCC + proteinlistMF)
    print(dataset_name)
    print(len(proteinlistset))
    BP_protein_vector_dict = present_protein_with_go_term_vector(proteinlistBP, geneAssociation.BPassociations,
                                                                 BPTermVectors, ontology.BP_terms,
                                                                 )
    CC_protein_vector_dict = present_protein_with_go_term_vector(proteinlistCC, geneAssociation.CCassociations,
                                                                 CCTermVectors, ontology.CC_terms,
                                                                 )
    MF_protein_vector_dict = present_protein_with_go_term_vector(proteinlistMF, geneAssociation.MFassociations,
                                                                 MFTermVectors, ontology.MF_terms)
    save_as_term_vector_for_expression_test(test_src_folder, test_out_folder, BP_protein_vector_dict, CC_protein_vector_dict, MF_protein_vector_dict)


def generate_for_sequence():
    # Deal for Sequence datasets
    home_path = os.path.abspath(os.path.join(os.getcwd(), ".."))
    dataset_name = "yeast"
    if config.IEA_FLAG == True:
        iea_name = "iea+"
    else:
        iea_name = "iea-"

    # 解析基因本体文件
    ontology = Ontology()
    ontology.parse_from_xml_file(config.ONTOLOGY_FILE)

    ASSOCIAION_FILE = config.YEAST_ASSOCIAION_FILE
    geneAssociation = GeneAssociation()
    geneAssociation.parse_association_file(ASSOCIAION_FILE, config.IEA_FLAG,ontology)

    # generate for train
    train_src_folder = home_path + "/datasets/raw_ppi/yeast/"+iea_name
    train_out_folder = home_path + "/intermediate_outputs/sequence_data/train/" + iea_name


    proteinlistBP, proteinlistCC, proteinlistMF = get_all_proteins_for_expression(train_src_folder)
    proteinlistset = set(proteinlistBP + proteinlistCC + proteinlistMF)
    print(dataset_name)
    print("the size of protein set: {}".format(len(proteinlistset)))


    # construct termVector
    BPTermVectors = TermVectors()
    BPTermVectors.parse_term_embedding_file(config.BP_TERM_EMB_FILE_PATH)
    BPTermVectors.construct_term_vector_dict()

    print("the size of go term set for BP: {}".format(len(BPTermVectors.termVectorDict)))

    CCTermVectors = TermVectors()
    CCTermVectors.parse_term_embedding_file(config.CC_TERM_EMB_FILE_PATH)
    CCTermVectors.construct_term_vector_dict()
    print("the size of go term set for CC: {}".format(len(CCTermVectors.termVectorDict)))

    MFTermVectors = TermVectors()
    MFTermVectors.parse_term_embedding_file(config.MF_TERM_EMB_FILE_PATH)
    MFTermVectors.construct_term_vector_dict()
    print("the size of go term set for MF: {}".format(len(MFTermVectors.termVectorDict)))


    print(train_src_folder)
    BP_protein_vector_dict = present_protein_with_go_term_vector(proteinlistBP, geneAssociation.BPassociations,
                                                                 BPTermVectors, ontology.BP_terms,
                                                                 )
    CC_protein_vector_dict = present_protein_with_go_term_vector(proteinlistCC, geneAssociation.CCassociations,
                                                                 CCTermVectors, ontology.CC_terms,
                                                                 )
    MF_protein_vector_dict = present_protein_with_go_term_vector(proteinlistMF, geneAssociation.MFassociations,
                                                                 MFTermVectors, ontology.MF_terms)

    # the train data for sequence data is same with the train data for expression data
    save_as_term_vector_for_expression_train(train_src_folder, train_out_folder, BP_protein_vector_dict, CC_protein_vector_dict, MF_protein_vector_dict)

    # generate for test
    test_src_folder = home_path + "/datasets/sequence_data/"
    test_out_folder = home_path + "/intermediate_outputs/sequence_data/test/" + iea_name
    proteinlistBP, proteinlistCC, proteinlistMF = get_all_proteins_for_sequence_test(test_src_folder,geneAssociation)
    proteinlistset = set(proteinlistBP + proteinlistCC + proteinlistMF)
    # construct termVector
    BPTermVectors = TermVectors()
    BPTermVectors.parse_term_embedding_file(config.BP_TERM_EMB_FILE_PATH)
    BPTermVectors.construct_term_vector_dict()

    print("the size of go term set for BP: {}".format(len(BPTermVectors.termVectorDict)))

    CCTermVectors = TermVectors()
    CCTermVectors.parse_term_embedding_file(config.CC_TERM_EMB_FILE_PATH)
    CCTermVectors.construct_term_vector_dict()
    print("the size of go term set for CC: {}".format(len(CCTermVectors.termVectorDict)))

    MFTermVectors = TermVectors()
    MFTermVectors.parse_term_embedding_file(config.MF_TERM_EMB_FILE_PATH)
    MFTermVectors.construct_term_vector_dict()
    print("the size of go term set for MF: {}".format(len(MFTermVectors.termVectorDict)))

    print(test_src_folder)
    print(len(proteinlistset))
    BP_protein_vector_dict = present_protein_with_go_term_vector(proteinlistBP, geneAssociation.BPassociations,
                                                                 BPTermVectors, ontology.BP_terms,
                                                                 )
    CC_protein_vector_dict = present_protein_with_go_term_vector(proteinlistCC, geneAssociation.CCassociations,
                                                                 CCTermVectors, ontology.CC_terms,
                                                                 )
    MF_protein_vector_dict = present_protein_with_go_term_vector(proteinlistMF, geneAssociation.MFassociations,
                                                                 MFTermVectors, ontology.MF_terms)
    save_as_term_vector_for_sequence_test(test_src_folder, test_out_folder, BP_protein_vector_dict, CC_protein_vector_dict, MF_protein_vector_dict)

if __name__ == '__main__':

    generate_for_ppi()
    #print("not fuound term size:{}".format(len(config.not_found_go_terms)))

    #generate_for_expression()

    #generate_for_sequence()











