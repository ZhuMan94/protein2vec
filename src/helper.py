import os
import numpy as np
import random
import xml.etree.ElementTree as ET

import src.config as config



# class Ontology:
#     def __init__(self):
#         self.all_terms=[]
#         self.BP_terms=[]
#         self.CC_terms=[]
#         self.MF_terms=[]
#
#     def parse_from_xml_file(self,ontology_file):
#         tree = ET.parse(ontology_file)
#         root = tree.getroot()
#         domain=""
#         term_id=""
#         alt_ids=[]
#         considers=[]
#         for term in root.iter('term'):
#             is_obsolete = False
#             for node in term:
#                 if node.tag == "namespace":
#                     domain = node.text
#                 if node.tag == "id":
#                     term_id = node.text
#                 if node.tag == "is_obsolete":
#                     is_obsolete=True
#             if is_obsolete==False:
#                 self.all_terms.append(term_id)
#                 if  domain=="biological_process":
#                     self.BP_terms.append(term_id)
#                 elif domain=="molecular_function":
#                     self.MF_terms.append(term_id)
#                 else:
#                     self.CC_terms.append(term_id)


class Ontology:
    def __init__(self):
        self.all_terms={}
        self.BP_terms={}
        self.CC_terms={}
        self.MF_terms={}

    def parse_from_xml_file(self,ontology_file):
        tree = ET.parse(ontology_file)
        root = tree.getroot()
        domain=""
        term_id=""
        alt_ids=[]
        considers=[]
        for term in root.iter('term'):
            is_obsolete = False
            for node in term:
                if node.tag == "namespace":
                    domain = node.text
                if node.tag == "id":
                    term_id = node.text
                if node.tag == "is_obsolete":
                    is_obsolete=True
            if is_obsolete==False:
                if  domain=="biological_process":
                    term=GOTerm(ID=term_id,domain='P')
                    self.BP_terms[term_id]=term
                elif domain=="molecular_function":
                    term = GOTerm(ID=term_id, domain='F')
                    self.MF_terms[term_id]=term
                else:
                    term = GOTerm(ID=term_id, domain='C')
                    self.CC_terms[term_id]=term
                self.all_terms[term_id] =term


class AbandonedGOTerm:
    def __init__(self,AbandonedID,domain,CurrentID):
        self.AbandonedID=AbandonedID
        self.domain=domain
        self.CurrentID=CurrentID

class GOTerm:
    def __init__(self,ID,domain):
        self.ID=ID
        self.domain=domain
        self.alt_ids=None
        self.considers=None

class GeneAssociation:
    def __init__(self):
        self.associations=dict()
        self.BPassociations=dict()
        self.MFassociations=dict()
        self.CCassociations=dict()
        self.BPGOtermNum=0
        self.CCGOtermNum =0
        self.MFGOtermNum =0
        self.geneID2protienID={}
        self.proteinID2geneID={}

    def get_geneID_proteinID_mapping(self,association_file):
        with open(association_file,'r') as f:
            lines=f.readlines()
            for line in lines:
                if line.startswith('!'):
                    continue
                if line.startswith('EcoGene'):
                    continue
                parts=line.split("\t")
                proteinID=parts[1]
                geneID=parts[2]
                if proteinID in self.proteinID2geneID:
                    if geneID not in self.proteinID2geneID[proteinID]:
                        self.proteinID2geneID[proteinID].append(geneID)
                else:
                    self.proteinID2geneID[proteinID]=[geneID]
                if geneID in self.geneID2protienID:
                    if proteinID not in self.geneID2protienID[geneID]:
                        self.geneID2protienID[geneID].append(proteinID)
                else:
                    self.geneID2protienID[geneID]=[proteinID]




    def parse_association_file(self,association_file,IEA_Flag,ontology,dataset_name):
        qualifiers=["NOT","not","colocalizes_with","Colocalizes_with","COLOCALIZES_WITH","contributes_to","Contributes_to","CONTRIBUTES_TO"]

        mf_proteins_with_contribute = {}
        with open(association_file,'r') as f:
            lines=f.readlines()
            for line in lines:
                if line.startswith('!'):
                    continue
                if line.startswith('EcoGene'):
                    continue
                old_parts = line.split("\t")
                parts=[]
                for part in old_parts:
                    if part!='':  # remove ''
                        parts.append(part)

                if dataset_name=="human":
                    proteinID = parts[2].strip() # dataset is provided with gene pair
                else:
                    proteinID = parts[1].strip()
                domain=None
                TermID=None
                evidence=None
                constraint =None
                for part in parts:
                    if part=='C' or part=='F' or part=='P':
                        if domain==None:
                            domain=part
                    elif part.startswith("GO:"):
                        if TermID==None:
                            TermID=part
                    elif part=='IEA':
                        evidence='IEA'
                    elif part.startswith("NOT") or part.startswith("not") or part in qualifiers:
                        constraint=part

                assert len(domain)==1,"domain: "+domain+",line:"+str(parts)
                assert len(TermID)==10,"TermID: "+TermID+",line:"+str(parts)
                assert evidence==None or evidence=="IEA","evidence"+evidence+"line:"+str(parts)
                assert TermID.startswith("GO:"),"Term Name:"+TermID #check TermID name

                if evidence=="IEA" and IEA_Flag==False: # without using IEA
                    continue
                if part.startswith("NOT") or part.startswith("not"):
                    continue
                if (constraint in qualifiers[2:5]):
                    continue
                if (constraint in qualifiers[5:8])and domain=='p':
                    continue
                if TermID not in ontology.all_terms:  # is_obsolete
                    continue
                if (constraint in qualifiers[5:8]) and domain == 'F':
                    if proteinID in mf_proteins_with_contribute:
                        mf_proteins_with_contribute[proteinID].append(TermID)
                    else:
                        mf_proteins_with_contribute[proteinID]=[TermID]
                    continue
                if proteinID in self.associations:
                    self.associations[proteinID].append(TermID)
                else:
                    self.associations[proteinID] = [TermID]
                if domain=="P":
                    #assert TermID in ontology.BP_terms, TermID+":"+domain+' not in BP'+"line:"+str(parts)
                    self.BPGOtermNum+=1
                    if proteinID in self.BPassociations:
                        self.BPassociations[proteinID].append(TermID)
                    else:
                        self.BPassociations[proteinID] = [TermID]
                if domain=="C":
                    #assert TermID in ontology.CC_terms,TermID+":"+domain+' not in CC'+"line:"+str(parts)
                    self.CCGOtermNum += 1
                    if proteinID in self.CCassociations:
                        self.CCassociations[proteinID].append(TermID)
                    else:
                        self.CCassociations[proteinID] = [TermID]
                if domain == "F":
                    #assert TermID in ontology.MF_terms, TermID+":"+domain+'ontology domain:'+ontology.all_terms[TermID].domain+"line:"+str(parts)
                    self.MFGOtermNum += 1
                    if proteinID in self.MFassociations:
                        self.MFassociations[proteinID].append(TermID)
                    else:
                        self.MFassociations[proteinID] = [TermID]

            for protien_id_with_contribute in mf_proteins_with_contribute:
                if protien_id_with_contribute  in self.CCassociations:
                    self.MFassociations[protien_id_with_contribute]=mf_proteins_with_contribute[protien_id_with_contribute]

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' + directory)

def get_all_proteins_for_expression_test(src_folder):
    proteinlistCC = []
    proteinlistBP = []
    proteinlistMF = []
    files=os.listdir(src_folder)
    for file in files:
        in_file_path=os.path.join(src_folder,file)
        with open(in_file_path, 'r') as f:
            lines = f.readlines()
            if file.find('.c.')>-1 or file.endswith(".c"):
                for line in lines:
                    if line.find(",")>-1:
                        protein1, protein2,v = line.strip('\n').strip(' ').split(',')
                    else:
                        protein1, protein2,v = line.strip('\n').strip(' ').split()
                    # print(protein1+' '+protein2)
                    # 去除protein2 后的换行符
                    if protein1 not in proteinlistCC:
                        proteinlistCC.append(protein1)
                    if protein2 not in proteinlistCC:
                        proteinlistCC.append(protein2)
            if file.find('.f.')>-1 or file.endswith(".f"):
                for line in lines:
                    if line.find(",")>-1:
                        protein1, protein2,v = line.strip('\n').strip(' ').split(',')
                    else:
                        protein1, protein2 ,v= line.strip('\n').strip(' ').split()
                    # print(protein1+' '+protein2)
                    # 去除protein2 后的换行符
                    if protein1 not in proteinlistMF:
                        proteinlistMF.append(protein1)
                    if protein2 not in proteinlistMF:
                        proteinlistMF.append(protein2)
            if file.find('.p.')>-1 or file.endswith(".p"):
                for line in lines:
                    if line.find(",") > -1:
                        protein1, protein2,v = line.strip('\n').strip(' ').split(',')
                    else:
                        protein1, protein2,v = line.strip('\n').strip(' ').split()
                    # print(protein1+' '+protein2)
                    # 去除protein2 后的换行符
                    if protein1 not in proteinlistBP:
                        proteinlistBP.append(protein1)
                    if protein2 not in proteinlistBP:
                        proteinlistBP.append(protein2)
    return proteinlistBP,proteinlistCC,proteinlistMF


def get_all_proteins_for_expression(src_folder):
    """
    :param dataset_name: "ecoli","human","yeast"
    :return:
    """
    proteinlistCC=[]
    proteinlistBP=[]
    proteinlistMF=[]
    files_negatives = os.listdir(src_folder+"/negatives/")
    files_positives = os.listdir(src_folder+"/positives/")
    files=files_negatives+files_positives
    for file in files:
        if not os.path.isdir(file):
            print(file)
            if file.startswith("negatives"):
                in_file_path = os.path.join(src_folder+"/negatives/", file)
            if file.startswith("positives"):
                in_file_path = os.path.join(src_folder + "/positives/", file)
            with open(in_file_path, 'r') as f:
                lines = f.readlines()
                if file.find('.c.')>-1 or file.endswith(".c"):
                    for line in lines:
                        if line.find(",")>-1:
                            protein1, protein2 = line.strip('\n').strip(' ').split(',')
                        else:
                            protein1, protein2 = line.strip('\n').strip(' ').split()
                        # print(protein1+' '+protein2)
                        # 去除protein2 后的换行符
                        if protein1 not in proteinlistCC:
                            proteinlistCC.append(protein1)
                        if protein2 not in proteinlistCC:
                            proteinlistCC.append(protein2)
                if file.find('.f.')>-1 or file.endswith(".f"):
                    for line in lines:
                        if line.find(",")>-1:
                            protein1, protein2 = line.strip('\n').strip(' ').split(',')
                        else:
                            protein1, protein2 = line.strip('\n').strip(' ').split()
                        # print(protein1+' '+protein2)
                        # 去除protein2 后的换行符
                        if protein1 not in proteinlistMF:
                            proteinlistMF.append(protein1)
                        if protein2 not in proteinlistMF:
                            proteinlistMF.append(protein2)
                if file.find('.p.')>-1 or file.endswith(".p"):
                    for line in lines:
                        if line.find(",") > -1:
                            protein1, protein2 = line.strip('\n').strip(' ').split(',')
                        else:
                            protein1, protein2 = line.strip('\n').strip(' ').split()
                        # print(protein1+' '+protein2)
                        # 去除protein2 后的换行符
                        if protein1 not in proteinlistBP:
                            proteinlistBP.append(protein1)
                        if protein2 not in proteinlistBP:
                            proteinlistBP.append(protein2)
    return proteinlistBP,proteinlistCC,proteinlistMF


def get_all_proteins(src_folder):
    """
    :param dataset_name: "ecoli","human","yeast"
    :return:
    """
    proteinlistCC=[]
    proteinlistBP=[]
    proteinlistMF=[]
    subfolders=["train","cross","test"]
    for subfolder in subfolders:
        files_negatives = os.listdir(src_folder+"/negatives/"+subfolder)
        files_positives = os.listdir(src_folder+"/positives/"+subfolder)
        files=files_negatives+files_positives
        for file in files:
            if not os.path.isdir(file):
                print(file)
                if file.startswith("negatives"):
                    in_file_path = os.path.join(src_folder+"/negatives/"+subfolder, file)
                if file.startswith("positives"):
                    in_file_path = os.path.join(src_folder + "/positives/"+subfolder, file)
                with open(in_file_path, 'r') as f:
                    lines = f.readlines()
                    if file.find('.c.')>-1 or file.endswith(".c"):
                        for line in lines:
                            if line.find(",")>-1:
                                print(line)
                                protein1, protein2 = line.strip('\n').strip(' ').split(',')
                            else:
                                protein1, protein2 = line.strip('\n').strip(' ').split()
                            # print(protein1+' '+protein2)
                            # 去除protein2 后的换行符
                            if protein1 not in proteinlistCC:
                                proteinlistCC.append(protein1)
                            if protein2 not in proteinlistCC:
                                proteinlistCC.append(protein2)
                    if file.find('.f.')>-1 or file.endswith(".f"):
                        for line in lines:
                            if line.find(",")>-1:
                                protein1, protein2 = line.strip('\n').strip(' ').split(',')
                            else:
                                protein1, protein2 = line.strip('\n').strip(' ').split()
                            # print(protein1+' '+protein2)
                            # 去除protein2 后的换行符
                            if protein1 not in proteinlistMF:
                                proteinlistMF.append(protein1)
                            if protein2 not in proteinlistMF:
                                proteinlistMF.append(protein2)
                    if file.find('.p.')>-1 or file.endswith(".p"):
                        for line in lines:
                            if line.find(",") > -1:
                                protein1, protein2 = line.strip('\n').strip(' ').split(',')
                            else:
                                protein1, protein2 = line.strip('\n').strip(' ').split()
                            # print(protein1+' '+protein2)
                            # 去除protein2 后的换行符
                            if protein1 not in proteinlistBP:
                                proteinlistBP.append(protein1)
                            if protein2 not in proteinlistBP:
                                proteinlistBP.append(protein2)
    return proteinlistBP,proteinlistCC,proteinlistMF


class TermVectors:
    def __init__(self):
        self.termVectorDict=dict()
        self.terms=None
        self.vectors=None

    def parse_term_embedding_file(self,file_path):
        with open(file_path) as f:
            lines=f.readlines()
            self.vectors=lines[1:] # 第一个是term的数量，不是termID
            terms = []
            for line in lines[1:]:
                term=line.split()[0]
                terms.append(term)
            self.terms = terms

    def str_to_vector(self,Str):
        Str = Str.strip('\n')
        nums = Str.split()
        vec = []
        for num in nums:
            vec.append(float(num))
        return np.array(vec)

    def construct_term_vector_dict(self):
        for term in self.terms:
            termindex = self.terms.index(term)
            line=self.vectors[termindex]
            s =line[len(term):].lstrip(' ')
            s_vec =self.str_to_vector(s)
            assert len(s_vec)==config.GODIM
            self.termVectorDict[term]=s_vec



def present_protein_with_go_term_vector(proteinlist,associations,termVectors,ontology_terms):
    """
    :param proteinlist:
    :param geneAssociation:
    :param termVectors:
    :return:
    """
    go_term_not_in_ontology=[]
    protein_vector_dict=dict()
    for protein in proteinlist:
        term_vector_array=[]
        if protein in associations:
            annoate_terms=associations[protein]
        else:
            continue
            print("{} protein not in association file".format(protein))
        for term in annoate_terms:
            term_ID=term[3:].lstrip('0')
            if term_ID in termVectors.termVectorDict:
                term_vector_array.append(termVectors.termVectorDict[term_ID])
            elif term in ontology_terms:
                print("{} in ontology".format(term))
            else:
                # GO term in annotation file can not be found or marked as is_obsolete in ontology file or
                #print("{} not in ontology".format(term))
                go_term_not_in_ontology.append(term_ID)
                continue
        protein_vector_dict[protein]=np.array(term_vector_array)
    print("the totoal num not in ontology:{}".format(len(set(go_term_not_in_ontology))))
    print(set(go_term_not_in_ontology))

    return protein_vector_dict


def save_protein_pair_as_term_vector_for_expression_test(file_path, protein_vector_dict,
                                                                 out_file_path):
    records = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            record = []
            if line.find(",")>-1:
                protein1, protein2,value = line.strip("\n").split(",")
            else:
                protein1, protein2,value = line.strip("\n").split("\t")
            if protein1 not in protein_vector_dict:
                print("{} protein key error".format(protein1))
                continue
            if protein2 not in protein_vector_dict:
                print("{} protein key error".format(protein2))
                continue
            protein1_vector = protein_vector_dict[protein1]
            protein2_vector = protein_vector_dict[protein2]

            if len(protein1_vector) == 0 or len(protein2_vector) == 0:
                print(line)
                continue
            record.append(protein1_vector)
            record.append(protein2_vector)
            record.append(float(value))
            record = np.array(record)
            records.append(record)
    records = np.array(records)
    np.save(out_file_path, records)

def save_protein_pair_as_term_vector_for_sequence_test(file_path, protein_vector_dict,
                                                                 out_file_path):
    records = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            record = []
            if line.find(",")>-1:
                protein1, protein2,LRBS,RRBS = line.strip("\n").split(",")
            else:
                protein1, protein2,LRBS,RRBS = line.strip("\n").split()
            if protein1 not in protein_vector_dict:
                print("{} protein key error".format(protein1))
                continue
            if protein2 not in protein_vector_dict:
                print("{} protein key error".format(protein2))
                continue
            protein1_vector = protein_vector_dict[protein1]
            protein2_vector = protein_vector_dict[protein2]

            if len(protein1_vector) == 0 or len(protein2_vector) == 0:
                print(line)
                continue
            record.append(protein1_vector)
            record.append(protein2_vector)
            record.append(float(LRBS))
            record.append(float(RRBS))
            record = np.array(record)
            records.append(record)
    records = np.array(records)
    np.save(out_file_path, records)


def save_protein_pair_as_term_vector(positive_file, negative_file,protein_vector_dict,out_file_path):
    records = []
    with open(positive_file,'r') as f:
       lines=f.readlines()
       for line in lines:
            record=[]
            if line.find(",")>-1:
                protein1, protein2 = line.strip("\n").split(",")
            else:
                protein1,protein2=line.split()
            if protein1 not in protein_vector_dict:
                print("{} protein key error".format(protein1))
                continue
            if protein2 not in protein_vector_dict:
                print("{} protein key error".format(protein2))
                continue
            protein1_vector=protein_vector_dict[protein1]
            protein2_vector=protein_vector_dict[protein2]

            if len(protein1_vector)==0 or len(protein2_vector)==0:
                print(line)
                continue
            record.append(protein1_vector)
            record.append(protein2_vector)
            record.append(1)
            record=np.array(record)
            records.append(record)
    with open(negative_file,'r') as f2:
       lines2=f2.readlines()
       for line in lines2:
            record_n=[]
            if line.find(",")>-1:
                protein1_n, protein2_n = line.strip("\n").split(",")
            else:
                protein1_n, protein2_n = line.split()
            if protein1_n not in protein_vector_dict:
                print("{} protein key error".format(protein1_n))
                continue
            if protein2_n not in protein_vector_dict:
                print("{} protein key error".format(protein2_n))
                continue
            protein1_vector_n=protein_vector_dict[protein1_n]
            protein2_vector_n=protein_vector_dict[protein2_n]
            if len(protein1_vector_n)==0 or len(protein2_vector_n)==0:
                print(line)
                continue
            record_n.append(protein1_vector_n)
            record_n.append(protein2_vector_n)
            record_n.append(0)
            record_n = np.array(record_n)
            records.append(record_n)
    records=np.array(records)
    np.save(out_file_path,records)

def save_all_protein_pairs_as_term_vector_for_fivefold(src_folder,out_folder,BP_protein_vector_dict,CC_protein_vector_dict,MF_protein_vector_dict):

    sub_folders=["train","cross","test"]
    for sub_folder in sub_folders:
        positives_folder = src_folder + "/positives/"+sub_folder
        negatives_folder = src_folder + "/negatives/"+sub_folder
        positives_files=os.listdir(positives_folder)
        for positive_file in positives_files:
            negative_file=positive_file.replace("positives","negatives")
            if positive_file.find(".f.")>-1:
                createFolder(out_folder + "/MolecularFunction/")
                out_file_path=out_folder+"/MolecularFunction/"+sub_folder+"."+positive_file.split(".")[-1]
                positive_file_path=os.path.join(positives_folder,positive_file)
                negative_file_path = os.path.join(negatives_folder, negative_file)
                save_protein_pair_as_term_vector(positive_file_path, negative_file_path,MF_protein_vector_dict,out_file_path)
            if positive_file.find(".p.")>-1:
                createFolder(out_folder+"/BiologicalProcess/")
                out_file_path=out_folder+"/BiologicalProcess/"+sub_folder+"."+positive_file.split(".")[-1]
                positive_file_path=os.path.join(positives_folder,positive_file)
                negative_file_path =os.path.join(negatives_folder, negative_file)
                save_protein_pair_as_term_vector(positive_file_path, negative_file_path,BP_protein_vector_dict,out_file_path)
            if positive_file.find(".c.")>-1:
                createFolder(out_folder + "/CellularComponent/")
                out_file_path=out_folder+"/CellularComponent/"+sub_folder+"."+positive_file.split(".")[-1]
                positive_file_path=os.path.join(positives_folder,positive_file)
                negative_file_path = os.path.join(negatives_folder, negative_file)
                save_protein_pair_as_term_vector(positive_file_path, negative_file_path,CC_protein_vector_dict,out_file_path)

def save_as_term_vector_for_expression_train(src_folder,out_folder,BP_protein_vector_dict,CC_protein_vector_dict,MF_protein_vector_dict):
    positives_folder = src_folder + "/positives/"
    negatives_folder = src_folder + "/negatives/"
    positives_files = os.listdir(positives_folder)
    for positive_file in positives_files:
        negative_file = positive_file.replace("positives", "negatives")
        if positive_file.endswith(".f"):
            createFolder(out_folder + "/MolecularFunction/")
            out_file_path = out_folder + "/MolecularFunction/"
            positive_file_path = os.path.join(positives_folder, positive_file)
            negative_file_path = os.path.join(negatives_folder, negative_file)
            save_protein_pair_as_term_vector(positive_file_path, negative_file_path, MF_protein_vector_dict,
                                             out_file_path)
        if positive_file.endswith(".p"):
            createFolder(out_folder + "/BiologicalProcess/")
            out_file_path = out_folder + "/BiologicalProcess/"
            positive_file_path = os.path.join(positives_folder, positive_file)
            negative_file_path = os.path.join(negatives_folder, negative_file)
            save_protein_pair_as_term_vector(positive_file_path, negative_file_path, BP_protein_vector_dict,
                                             out_file_path)
        if positive_file.endswith(".c"):
            createFolder(out_folder + "/CellularComponent/")
            out_file_path = out_folder + "/CellularComponent/"
            positive_file_path = os.path.join(positives_folder, positive_file)
            negative_file_path = os.path.join(negatives_folder, negative_file)
            save_protein_pair_as_term_vector(positive_file_path, negative_file_path, CC_protein_vector_dict,
                                             out_file_path)

def save_as_term_vector_for_expression_test(src_folder,out_folder,BP_protein_vector_dict,CC_protein_vector_dict,MF_protein_vector_dict):

    files = os.listdir(src_folder)
    for file in files:
        if file.endswith(".f"):
            createFolder(out_folder + "/MolecularFunction/")
            out_file_path = out_folder + "/MolecularFunction/"
            file_path = os.path.join(src_folder, file)
            save_protein_pair_as_term_vector_for_expression_test(file_path, MF_protein_vector_dict,
                                             out_file_path)
        if file.endswith(".p"):
            createFolder(out_folder + "/BiologicalProcess/")
            out_file_path = out_folder + "/BiologicalProcess/"
            file_path = os.path.join(src_folder, file)
            save_protein_pair_as_term_vector_for_expression_test(file_path, BP_protein_vector_dict,
                                                                 out_file_path)
        if file.endswith(".c"):
            createFolder(out_folder + "/CellularComponent/")
            out_file_path = out_folder + "/CellularComponent/"
            file_path = os.path.join(src_folder, file)
            save_protein_pair_as_term_vector_for_expression_test(file_path, CC_protein_vector_dict,
                                                                 out_file_path)


def save_all_protein_pairs_as_term_vector(src_folder,out_folder,BP_protein_vector_dict,CC_protein_vector_dict,MF_protein_vector_dict):

    sub_folders=["train","cross","test"]
    for sub_folder in sub_folders:
        positives_folder = src_folder + "/positives/"+sub_folder
        negatives_folder = src_folder + "/negatives/"+sub_folder
        positives_files=os.listdir(positives_folder)
        for positive_file in positives_files:
            negative_file=positive_file.replace("positives","negatives")
            if positive_file.endswith(".f"):
                createFolder(out_folder + "/MolecularFunction/")
                out_file_path=out_folder+"/MolecularFunction/"+sub_folder
                positive_file_path=os.path.join(positives_folder,positive_file)
                negative_file_path = os.path.join(negatives_folder, negative_file)
                save_protein_pair_as_term_vector(positive_file_path, negative_file_path,MF_protein_vector_dict,out_file_path)
            if positive_file.endswith(".p"):
                createFolder(out_folder+"/BiologicalProcess/")
                out_file_path=out_folder+"/BiologicalProcess/"+sub_folder
                positive_file_path=os.path.join(positives_folder,positive_file)
                negative_file_path =os.path.join(negatives_folder, negative_file)
                save_protein_pair_as_term_vector(positive_file_path, negative_file_path,BP_protein_vector_dict,out_file_path)
            if positive_file.endswith(".c"):
                createFolder(out_folder + "/CellularComponent/")
                out_file_path=out_folder+"/CellularComponent/"+sub_folder
                positive_file_path=os.path.join(positives_folder,positive_file)
                negative_file_path = os.path.join(negatives_folder, negative_file)
                save_protein_pair_as_term_vector(positive_file_path, negative_file_path,CC_protein_vector_dict,out_file_path)


def load_data(file_path):
    records=np.load(file_path)
    random_idexes=[i for i in range(len(records))]
    random.shuffle(random_idexes)
    protein1_array = []
    protein2_array = []
    label_array = []
    for index in random_idexes:
        protein1 = records[index][0]
        protein2 = records[index][1]
        label = records[index][2]
        protein1_array.append(protein1)
        protein2_array.append(protein2)
        label_array.append(label)
    return protein1_array, protein2_array, label_array

def get_all_proteins_for_sequence_test(src_folder,geneAssociation):
    sequence_file=src_folder+"/sequence_homology"
    proteinlistBP=[]
    proteinlistCC= []
    proteinlistMF=[]
    with open(sequence_file, 'r') as f:
        lines = f.readlines()
        for line in lines[1:]:
            if line.find(",") > -1:
                protein1, protein2, LRBS,RRBS = line.strip('\n').strip(' ').split(',')
            else:
                protein1, protein2, LRBS,RRBS = line.strip('\n').strip(' ').split()
            # print(protein1+' '+protein2)
            # 去除protein2 后的换行符
            if protein1 in geneAssociation.BPassociations and protein2 in geneAssociation.BPassociations:
                if protein1 not in proteinlistBP:
                    proteinlistBP.append(protein1)
                if protein2 not in proteinlistBP:
                    proteinlistBP.append(protein2)
            if protein1 in geneAssociation.CCassociations and protein2 in geneAssociation.CCassociations:
                if protein1 not in proteinlistCC:
                    proteinlistCC.append(protein1)
                if protein2 not in proteinlistCC:
                    proteinlistCC.append(protein2)
            if protein1 in geneAssociation.MFassociations and protein2 in geneAssociation.MFassociations:
                if protein1 not in proteinlistMF:
                    proteinlistMF.append(protein1)
                if protein2 not in proteinlistMF:
                    proteinlistMF.append(protein2)
        return proteinlistBP, proteinlistCC, proteinlistMF


def save_as_term_vector_for_sequence_test(src_folder, out_folder, BP_protein_vector_dict,
                                          CC_protein_vector_dict, MF_protein_vector_dict):
    files = os.listdir(src_folder)
    for file in files:
        file_path=src_folder+"/"+file
        createFolder(out_folder + "/MolecularFunction/")
        out_file_path = out_folder + "/MolecularFunction/"
        file_path = os.path.join(src_folder, file)
        save_protein_pair_as_term_vector_for_sequence_test(file_path, MF_protein_vector_dict,
                                                             out_file_path)

        createFolder(out_folder + "/BiologicalProcess/")
        out_file_path = out_folder + "/BiologicalProcess/"
        file_path = os.path.join(src_folder, file)
        save_protein_pair_as_term_vector_for_sequence_test(file_path, BP_protein_vector_dict,
                                                             out_file_path)

        createFolder(out_folder + "/CellularComponent/")
        out_file_path = out_folder + "/CellularComponent/"
        file_path = os.path.join(src_folder, file)
        save_protein_pair_as_term_vector_for_sequence_test(file_path, CC_protein_vector_dict,
                                                             out_file_path)

























