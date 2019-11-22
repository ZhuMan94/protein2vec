import src.config as config
import os
import xml.etree.ElementTree as ET
#from src.run_node2vec import run_node2vec


def parseFromXml(file):
    """
    parse ontology file, and save the edges of ontology file to a edgelist file
    :param file: ontology file
    """
    tree = ET.parse(file)
    root = tree.getroot()

    if not os.path.exists(config.ONTOLOGY_EDGELIST_FOLDER):
        os.makedirs(config.ONTOLOGY_EDGELIST_FOLDER)
    CC_edgelist_file=config.ONTOLOGY_EDGELIST_FOLDER+"/CC.edgelist"
    BP_edgelist_file = config.ONTOLOGY_EDGELIST_FOLDER + "/BP.edgelist"
    MF_edgelist_file = config.ONTOLOGY_EDGELIST_FOLDER + "/MF.edgelist"
    of_CC = open(CC_edgelist_file, 'w')
    of_BP = open(BP_edgelist_file, 'w')
    of_MF = open(MF_edgelist_file, 'w')
    for term in root.iter('term'):
        is_a = []
        flag = 0
        part_of = []
        is_obsolete_falg = 0
        consider = []
        regulates = []
        for node in term:
            # print(node.tag)
            # print(node.text)
            if node.tag=="namespace":
                domain=node.text
            elif node.tag == "id":
                go_id = node.text
            elif node.tag == "is_a":
                is_a.append(node.text)
            elif node.tag == "is_obsolete":
                if node.text == "1":
                    is_obsolete_falg = node.text
            elif node.tag == "consider":
                consider.append(node.text)
            elif node.tag == "relationship":
                for relation in node:
                    if relation.tag == "type":
                        if relation.text == "part_of":
                            flag = 1
                        if relation.text == "regulates":
                            flag = 2
                    if relation.tag == "to":
                        if flag == 1:
                            part_of.append(relation.text)
                            flag = 0
                        if flag == 2:
                            regulates.append(relation.text)
                            flag = 0

        if is_obsolete_falg == 0:
            for parent in is_a:
                if domain=="biological_process":
                    of_BP.write(parent[3:] + " " + go_id[3:] + '\n')
                elif domain=="molecular_function":
                    of_MF.write(parent[3:] + " " + go_id[3:] + '\n')
                else:
                    of_CC.write(parent[3:] + " " + go_id[3:] + '\n')
            for parent in part_of:
                if domain == "biological_process":
                    of_BP.write(parent[3:] + " " + go_id[3:] + '\n')
                elif domain == "molecular_function":
                    of_MF.write(parent[3:] + " " + go_id[3:] + '\n')
                else:
                    of_CC.write(parent[3:] + " " + go_id[3:] + '\n')
        else:
           continue
    of_CC.close()
    of_BP.close()
    of_MF.close()



if __name__ == '__main__':
    #Step1: parse gene ontology file and generate edgelist files for CC,BP,MF
    parseFromXml(config.ONTOLOGY_FILE)

    # The term encoding module employs a node2vec model to generate dense feature vectors for GO terms
    # save the feature vector representation of GO terms in *.emb file
   # run_node2vec()

