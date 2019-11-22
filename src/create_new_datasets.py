from src.helper import  *

def deal_expression_data(src_folder,out_folder,geneAssociation):
    files = os.listdir(src_folder)
    for file in files:
        in_file_path = os.path.join(src_folder, file)
        out_file_path=os.path.join(out_folder,file)
        of=open(out_file_path,'w')
        with open(in_file_path, 'r') as f:
            lines = f.readlines()
            if file.find('.c.') > -1 or file.endswith(".c"):
                for line in lines:
                    if line.find(",") > -1:
                        protein1, protein2, v = line.strip('\n').strip(' ').split(',')
                    else:
                        protein1, protein2, v = line.strip('\n').strip(' ').split()
                    # print(protein1+' '+protein2)
                    # 去除protein2 后的换行符
                    if protein1 in geneAssociation.CCassociations and protein2  in geneAssociation.CCassociations:
                       of.write(line)
                       of.flush()
            if file.find('.f.') > -1 or file.endswith(".f"):
                for line in lines:
                    if line.find(",") > -1:
                        protein1, protein2, v = line.strip('\n').strip(' ').split(',')
                    else:
                        protein1, protein2, v = line.strip('\n').strip(' ').split()
                    if protein1 in geneAssociation.MFassociations and protein2  in geneAssociation.MFassociations:
                       of.write(line)
                       of.flush()
            if file.find('.p.') > -1 or file.endswith(".p"):
                for line in lines:
                    if line.find(",") > -1:
                        protein1, protein2, v = line.strip('\n').strip(' ').split(',')
                    else:
                        protein1, protein2, v = line.strip('\n').strip(' ').split()
                    # print(protein1+' '+protein2)
                    # 去除protein2 后的换行符
                    if protein1 in geneAssociation.BPassociations and protein2 in geneAssociation.BPassociations:
                       of.write(line)
                       of.flush()

def deal_sequence_data(src_folder,out_folder,geneAssociation):
    files = os.listdir(src_folder)
    for file in files:
        in_file_path = os.path.join(src_folder, file)
        out_file_path=os.path.join(out_folder,file)
        of=open(out_file_path,'w')
        with open(in_file_path, 'r') as f:
            lines = f.readlines()
            of.write(lines[0])
            for line in lines[1:]:
                if line.find(",") > -1:
                    protein1, protein2, lrbs,rrbs = line.strip('\n').strip(' ').split(',')
                else:
                    protein1, protein2, lrbs,rrbs = line.strip('\n').strip(' ').split()
                # print(protein1+' '+protein2)
                # 去除protein2 后的换行符
                if protein1 in geneAssociation.associations and protein2  in geneAssociation.associations:
                   of.write(line)
                   of.flush()

def get_all_ppis_from_split_data(src_folder,geneAssociation):
    """
    :param dataset_name: "ecoli","human","yeast"
    :return:
    """
    ppisCC=[]
    ppisBP=[]
    ppisMF=[]
    subfolders=["train","cross","test"]
    for subfolder in subfolders:
        files_negatives = os.listdir(src_folder+"/negatives/"+subfolder)
        files_positives = os.listdir(src_folder+"/positives/"+subfolder)
        files=files_negatives+files_positives
        for file in files:
            if not os.path.isdir(file):
                #print(file)
                if file.startswith("negatives"):
                    in_file_path = os.path.join(src_folder+"/negatives/"+subfolder, file)
                if file.startswith("positives"):
                    in_file_path = os.path.join(src_folder + "/positives/"+subfolder, file)
                with open(in_file_path, 'r') as f:
                    lines = f.readlines()
                    if file.find('.c.')>-1 or file.endswith(".c"):
                        for line in lines:
                            if line.find(",")>-1:
                                protein1, protein2 = line.strip('\n').strip(' ').split(',')
                            else:
                                protein1, protein2 = line.strip('\n').strip(' ').split()
                            # if [protein1,protein2] not in ppisCC:
                            #     ppisCC.append([protein1,protein2])
                            if protein1 in geneAssociation.CCassociations and protein2 in geneAssociation.CCassociations:
                                if protein1+protein2 not in ppisCC:
                                    ppisCC.append(protein1+protein2)
                            else:
                                print("prototein pair is removed:{}".format(line))

                    if file.find('.f.')>-1 or file.endswith(".f"):
                        for line in lines:
                            if line.find(",")>-1:
                                protein1, protein2 = line.strip('\n').strip(' ').split(',')
                            else:
                                protein1, protein2 = line.strip('\n').strip(' ').split()
                            # if [protein1, protein2] not in ppisMF:
                            #     ppisMF.append([protein1, protein2])
                            if protein1 in geneAssociation.MFassociations and protein2 in geneAssociation.MFassociations:
                                if protein1+protein2 not in ppisMF:
                                    ppisMF.append(protein1+protein2)
                            else:
                                print("prototein pair is removed:{}".format(line))
                    if file.find('.p.')>-1 or file.endswith(".p"):
                        for line in lines:
                            if line.find(",") > -1:
                                protein1, protein2 = line.strip('\n').strip(' ').split(',')
                            else:
                                protein1, protein2 = line.strip('\n').strip(' ').split()
                            if protein1 in geneAssociation.BPassociations and protein2 in geneAssociation.BPassociations:
                                if protein1+protein2 not in ppisBP:
                                    ppisBP.append(protein1+protein2)
                            else:
                                print("prototein pair is removed:{}".format(line))
    return ppisBP,ppisCC,ppisMF

def stat_ppi_num():
    # deal ppi data
    home_path = os.path.abspath(os.path.join(os.getcwd(), ".."))
    all_bp_ppis = []
    all_cc_ppis = []
    all_mf_ppis = []
    dataset_names = ["human"]  # "yeast", "ecoli",
    split_modes = ["connectivity", "fivefold"]
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
            geneAssociation.parse_association_file(ASSOCIAION_FILE, config.IEA_FLAG, ontology)
            print("The number of annotation for BP CC MF")
            print(len(geneAssociation.BPassociations))
            print(len(geneAssociation.CCassociations))
            print(len(geneAssociation.MFassociations))

            src_folder = home_path + "/datasets/processed_ppi/" + split_mode + "/" + dataset_name + "/" + iea_name
            print(src_folder)
            ppisBP, ppisCC, ppisMF = get_all_ppis_from_split_data(src_folder,geneAssociation)

            all_bp_ppis = all_bp_ppis + ppisBP
            all_cc_ppis = all_cc_ppis + ppisCC
            all_mf_ppis = all_mf_ppis + ppisMF

            print("dataset name:{};iea name:{},split_mode:{},CC size:{}".format(dataset_name, iea_name, split_mode,
                                                                                len(ppisCC)))
            print("dataset name:{};iea name:{},split_mode:{},BP size:{}".format(dataset_name, iea_name, split_mode,
                                                                                len(ppisBP)))
            print("dataset name:{};iea name:{},split_mode:{},MF size:{}".format(dataset_name, iea_name, split_mode,
                                                                                len(ppisMF)))

        print("ALL: dataset name:{};iea name:{},split_mode:{},CC size:{}".format(dataset_name, iea_name, split_mode,
                                                                                 len(set(all_cc_ppis))))
        print("ALL: dataset name:{};iea name:{},split_mode:{},BP size:{}".format(dataset_name, iea_name, split_mode,
                                                                                 len(set(all_bp_ppis))))
        print("ALL: dataset name:{};iea name:{},split_mode:{},MF size:{}".format(dataset_name, iea_name, split_mode,
                                                                                 len(set(all_mf_ppis))))


def deal_ppi(src_folder,out_folder,geneAssociation):
    """
    :param dataset_name: "ecoli","human","yeast"
    :return:
    """
    ppisCC=[]
    ppisBP=[]
    ppisMF=[]
    subfolders=["train","cross","test"]
    for subfolder in subfolders:
        files_negatives = os.listdir(src_folder+"/negatives/"+subfolder)
        files_positives = os.listdir(src_folder+"/positives/"+subfolder)
        files=files_negatives+files_positives
        for file in files:
            if not os.path.isdir(file):
                #print(file)
                if file.startswith("negatives"):
                    in_file_path = os.path.join(src_folder+"/negatives/"+subfolder, file)
                    if not os.path.exists(out_folder+"/negatives/"+subfolder):
                        os.makedirs(out_folder+"/negatives/"+subfolder)
                    out_file_path=os.path.join(out_folder+"/negatives/"+subfolder, file)
                if file.startswith("positives"):
                    in_file_path = os.path.join(src_folder + "/positives/"+subfolder, file)
                    if not os.path.exists(out_folder + "/positives/" + subfolder):
                        os.makedirs(out_folder + "/positives/" + subfolder)
                    out_file_path = os.path.join(out_folder + "/positives/" + subfolder, file)
                of=open(out_file_path,'w')
                with open(in_file_path, 'r') as f:
                    lines = f.readlines()
                    if file.find('.c.')>-1 or file.endswith(".c"):
                        for line in lines:
                            if line.find(",")>-1:
                                protein1, protein2 = line.strip('\n').strip(' ').split(',')
                            else:
                                protein1, protein2 = line.strip('\n').strip(' ').split()

                            if protein1 in geneAssociation.CCassociations and protein2 in geneAssociation.CCassociations:
                                 of.write(line)
                                 of.flush()
                                 if [protein1,protein2] not in ppisCC:
                                      ppisCC.append([protein1,protein2])
                            else:
                                print("prototein pair is removed:{}".format(line))

                    if file.find('.f.')>-1 or file.endswith(".f"):
                        for line in lines:
                            if line.find(",")>-1:
                                protein1, protein2 = line.strip('\n').strip(' ').split(',')
                            else:
                                protein1, protein2 = line.strip('\n').strip(' ').split()
                            # if [protein1, protein2] not in ppisMF:
                            #     ppisMF.append([protein1, protein2])
                            if protein1 in geneAssociation.MFassociations and protein2 in geneAssociation.MFassociations:
                                of.write(line)
                                of.flush()
                                if [protein1, protein2] not in ppisMF:
                                    ppisMF.append([protein1, protein2])
                            else:
                                print("prototein pair is removed:{}".format(line))
                    if file.find('.p.')>-1 or file.endswith(".p"):
                        for line in lines:
                            if line.find(",") > -1:
                                protein1, protein2 = line.strip('\n').strip(' ').split(',')
                            else:
                                protein1, protein2 = line.strip('\n').strip(' ').split()
                            if protein1 in geneAssociation.BPassociations and protein2 in geneAssociation.BPassociations:
                                of.write(line)
                                of.flush()
                                if [protein1, protein2] not in ppisBP:
                                    ppisBP.append([protein1, protein2])
                            else:
                                print("prototein pair is removed:{}".format(line))
    return ppisBP,ppisCC,ppisMF

def deal_ppi_data():
    # deal ppi data
    home_path = os.path.abspath(os.path.join(os.getcwd(), ".."))
    dataset_names = ["human","yeast",'ecoli']  # "yeast", "ecoli",
    split_modes = ["connectivity", "fivefold"]
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
            geneAssociation.parse_association_file(ASSOCIAION_FILE, config.IEA_FLAG, ontology)
            print("The number of annotation for BP CC MF")
            print(len(geneAssociation.BPassociations))
            print(len(geneAssociation.CCassociations))
            print(len(geneAssociation.MFassociations))

            src_folder = home_path + "/datasets/processed_ppi/" + split_mode + "/" + dataset_name + "/" + iea_name
            out_folder = home_path + "/new_datasets/processed_ppi/" + split_mode + "/" + dataset_name + "/" + iea_name
            print(src_folder)
            ppisBP, ppisCC, ppisMF = deal_ppi(src_folder,out_folder,geneAssociation)
            print("dataset name:{};iea name:{},split_mode:{},CC size:{}".format(dataset_name, iea_name, split_mode,
                                                                                len(ppisCC)))
            print("dataset name:{};iea name:{},split_mode:{},BP size:{}".format(dataset_name, iea_name, split_mode,
                                                                                len(ppisBP)))
            print("dataset name:{};iea name:{},split_mode:{},MF size:{}".format(dataset_name, iea_name, split_mode,
                                                                                len(ppisMF)))

def deal_raw_ppis(src_folder,out_folder,geneAssociation):
    """
    :param dataset_name: "ecoli","human","yeast"
    :return:
    """
    ppisCC=[]
    ppisBP=[]
    ppisMF=[]
    files_negatives = os.listdir(src_folder+"/negatives/")
    files_positives = os.listdir(src_folder+"/positives/")
    files=files_negatives+files_positives
    for file in files:
        if not os.path.isdir(file):
            #print(file)
            if file.startswith("negatives"):
                in_file_path = os.path.join(src_folder+"/negatives/", file)
                if not os.path.exists(out_folder + "/negatives/"):
                    os.makedirs(out_folder + "/negatives/")
                out_file_path = os.path.join(out_folder + "/negatives/", file)
            if file.startswith("positives"):
                in_file_path = os.path.join(src_folder + "/positives/", file)
                if not os.path.exists(out_folder + "/positives/"):
                    os.makedirs(out_folder + "/positives/")
                out_file_path = os.path.join(out_folder + "/positives/", file)
            of=open(out_file_path,'w')
            with open(in_file_path, 'r') as f:
                lines = f.readlines()
                if file.find('.c.')>-1 or file.endswith(".c"):
                    for line in lines:
                        if line.find(",")>-1:
                            protein1, protein2 = line.strip('\n').strip(' ').split(',')
                        else:
                            protein1, protein2 = line.strip('\n').strip(' ').split()

                        if protein1 in geneAssociation.CCassociations and protein2 in geneAssociation.CCassociations:
                            if [protein1,protein2] not in ppisCC:
                                ppisCC.append([protein1,protein2])
                            of.write(line)
                            of.flush()
                        else:
                            print("prototein pair is removed:{}".format(line))

                if file.find('.f.')>-1 or file.endswith(".f"):
                    for line in lines:
                        if line.find(",")>-1:
                            protein1, protein2 = line.strip('\n').strip(' ').split(',')
                        else:
                            protein1, protein2 = line.strip('\n').strip(' ').split()
                        # if [protein1, protein2] not in ppisMF:
                        #     ppisMF.append([protein1, protein2])
                        if protein1 in geneAssociation.MFassociations and protein2 in geneAssociation.MFassociations:
                            if [protein1, protein2] not in ppisMF:
                                ppisMF.append([protein1, protein2])
                            of.write(line)
                            of.flush()
                        else:
                            print("prototein pair is removed:{}".format(line))
                if file.find('.p.')>-1 or file.endswith(".p"):
                    for line in lines:
                        if line.find(",") > -1:
                            protein1, protein2 = line.strip('\n').strip(' ').split(',')
                        else:
                            protein1, protein2 = line.strip('\n').strip(' ').split()
                        if protein1 in geneAssociation.BPassociations and protein2 in geneAssociation.BPassociations:
                            if [protein1, protein2] not in ppisBP:
                                ppisBP.append([protein1, protein2])
                            of.write(line)
                            of.flush()
                        else:
                            print("prototein pair is removed:{}".format(line))
    return ppisBP,ppisCC,ppisMF



def deal_raw_ppi():
    # deal ppi data
    home_path = os.path.abspath(os.path.join(os.getcwd(), ".."))
    dataset_names = ["human", "yeast", 'ecoli']  # "yeast", "ecoli",
    for dataset_name in dataset_names:

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
        geneAssociation.parse_association_file(ASSOCIAION_FILE, config.IEA_FLAG, ontology)
        print("The number of annotation for BP CC MF")
        print(len(geneAssociation.BPassociations))
        print(len(geneAssociation.CCassociations))
        print(len(geneAssociation.MFassociations))

        src_folder = home_path + "/datasets/raw_ppi/" +dataset_name + "/" + iea_name
        out_folder = home_path + "/new_datasets/raw_ppi/"+ dataset_name + "/" + iea_name
        ppisBP,ppisCC,ppisMF=deal_raw_ppis(src_folder, out_folder, geneAssociation)


if __name__ == '__main__':
    home_path = os.path.abspath(os.path.join(os.getcwd(), ".."))
    dataset_name = "yeast"
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

    ASSOCIAION_FILE = config.YEAST_ASSOCIAION_FILE
    geneAssociation = GeneAssociation()
    geneAssociation.parse_association_file(ASSOCIAION_FILE, config.IEA_FLAG, ontology)


    # deal expression data
    # src_folder=home_path+"/datasets/expression_data"
    # out_folder = home_path + "/new_datasets/expression_data"
    # if not os.path.exists(out_folder):
    #     os.makedirs(out_folder)
    # deal_expression_data(src_folder, out_folder, geneAssociation)
    #
    # # deal sequence data
    # src_folder = home_path + "/datasets/sequence_data"
    # out_folder = home_path + "/new_datasets/sequence_data"
    # if not os.path.exists(out_folder):
    #     os.makedirs(out_folder)
    # deal_sequence_data(src_folder, out_folder, geneAssociation)

    #deal_ppi_data()

    deal_raw_ppi()





