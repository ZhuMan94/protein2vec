from src.helper import  *

def stat_num(geneAssociation,file_name,type):
    """
     统计原始数据集的数目
    :param geneAssociation:
    :param file_name:
    :return:
    """
    protein_set=[]
    ppis=[]
    with open(file_name,'r') as f:
        lines=f.readlines()
        if type == "sequence_data":
            for line in lines:
                p1,p2,LRBS,RRBS=line.strip("\n").split()
                if p1 in geneAssociation and p2 in geneAssociation:
                    ppis.append([p1,p2])
                    if p1 not in protein_set:
                        protein_set.append(p1)
                    if p2 not in protein_set:
                        protein_set.append(p2)
                else:
                    print("Not in geneAssociation:{}".format(line))
        if type == "ppi_data":
            for line in lines:
                if line.find(","):
                    p1,p2=line.strip("\n").split(",")
                else:
                    p1, p2 = line.strip("\n").split()
                if p1 in geneAssociation and p2 in geneAssociation:
                    ppis.append([p1,p2])
                    if p1 not in protein_set:
                        protein_set.append(p1)
                    if p2 not in protein_set:
                        protein_set.append(p2)
                else:
                    print("Not in geneAssociation:{}".format(line))
        if type == "expression_data":
            for line in lines:
                p1,p2,v=line.strip("\n").split()
                if p1 in geneAssociation and p2 in geneAssociation:
                    ppis.append([p1,p2])
                else:
                    print("Not in geneAssociation:{}".format(line))
    return ppis,protein_set


def check_repeats(train_file,test_file):
    train_ppis=[]
    with open(train_file,'r') as f:
       lines=f.readlines()
       for line in lines:
           ppis=line.split()
           train_ppis.append((ppis[0],ppis[1]))

    with open(test_file,'r') as f:
        lines=f.readlines()
        cnt=-1
        for line in lines:
            ppis = line.split()
            cnt=cnt+1
            if (ppis[0],ppis[1]) in train_ppis or (ppis[1],ppis[0]) in train_ppis:
                print(cnt)
                print(line)
            else:
                continue

    return False

def check_connectivity_train_test():
    folder=home_path+"/datasets/processed_ppi/connectivity"
    dataset_names=['human','ecoli','yeast']
    iea_names=['iea+','iea-']
    subfoders=['negatives','positives']
    suffixs=['p','c','f']
    for dataset_name in dataset_names:
        for iea_name in iea_names:
            for subfoder in subfoders:
                for suffix in suffixs:
                     train_file=folder+"/"+dataset_name+"/"+iea_name+"/"+subfoder+"/train/"+subfoder+"."+dataset_name+"."+suffix
                     test_file =folder+"/"+dataset_name+"/"+iea_name+"/"+subfoder+"/test/"+subfoder+"."+dataset_name+"."+suffix
                     print(train_file)
                     assert check_repeats(train_file,test_file)==False

def check_connectivity_train_test():
    folder=home_path+"/datasets/processed_ppi/connectivity"
    dataset_names=['human','ecoli','yeast']
    iea_names=['iea+','iea-']
    subfoders=['negatives','positives']
    suffixs=['p','c','f']
    for dataset_name in dataset_names:
        for iea_name in iea_names:
            for subfoder in subfoders:
                for suffix in suffixs:
                     train_file=folder+"/"+dataset_name+"/"+iea_name+"/"+subfoder+"/train/"+subfoder+"."+dataset_name+"."+suffix
                     test_file =folder+"/"+dataset_name+"/"+iea_name+"/"+subfoder+"/test/"+subfoder+"."+dataset_name+"."+suffix
                     print(train_file)
                     assert check_repeats(train_file,test_file)==False

def check_single_repeats(file):
    ppis =[]
    with open(file) as f:
        lines=f.readlines()
        for line in lines:
            if line.find(',')>-1:
                p1,p2=line.split(",")
            else:
                p1,p2=line.split()
            if (p1,p2) in ppis or (p2,p1) in ppis:
                continue
            else:
                ppis.append((p1,p2))
    return len(ppis),len(lines)



def check_raw_ppi_repeats():
    folder = home_path + "/datasets/raw_ppi"
    dataset_names = ['human', 'ecoli', 'yeast']
    iea_names = ['iea+', 'iea-']
    subfoders = ['negatives', 'positives']
    suffixs = ['p', 'c', 'f']
    for dataset_name in dataset_names:
        for iea_name in iea_names:
            for subfoder in subfoders:
                    folder_path = folder + "/" + dataset_name + "/" + iea_name + "/" + subfoder
                    files=os.listdir(folder_path)
                    for file in files:
                        file_path=folder_path+"/"+file
                        ppi_len,file_len=check_single_repeats(file_path)
                        if ppi_len!=file_len:
                            print(file_path)
                            print("ppi len:{}".format(ppi_len))
                            print("file len{}".format(file_len))


if __name__ == '__main__':
    home_path = os.path.abspath(os.path.join(os.getcwd(), ".."))
    #check_connectivity_train_test()
    #check_raw_ppi_repeats()

    # step1: 获得注释数据
    home_path = os.path.abspath(os.path.join(os.getcwd(), ".."))
    dataset_names=["yeast",'human','ecoli']
    for dataset_name in dataset_names:
        if config.IEA_FLAG==True:
            iea_name="iea+"
        else:
            iea_name="iea-"
        if dataset_name=="ecoli":
            ASSOCIAION_FILE=config.ECOLI_ASSOCIAION_FILE
        if dataset_name=="yeast":
            ASSOCIAION_FILE=config.YEAST_ASSOCIAION_FILE
        if dataset_name=="human":
            ASSOCIAION_FILE=config.HUMAN_ASSOCIAION_FILE

        # 解析基因本体文件
        ontology = Ontology()
        ontology.parse_from_xml_file(config.ONTOLOGY_FILE)

        geneAssociation = GeneAssociation()
        geneAssociation.parse_association_file(ASSOCIAION_FILE,config.IEA_FLAG,ontology,dataset_name)

        if dataset_name=='yeast':
            file_name_prefix = home_path + "/datasets/raw_ppi/" + dataset_name + "/" + iea_name + "/positives/" + "positives.sgd.iea"
        else:
            file_name_prefix = home_path + "/datasets/raw_ppi/" + dataset_name + "/" + iea_name + "/positives/" + "positives." + dataset_name

        bp_file_name = file_name_prefix + ".p"
        bp_ppis,bp_protein_set=stat_num(geneAssociation.BPassociations,bp_file_name,"ppi_data")

        cc_file_name = file_name_prefix + ".c"
        cc_ppis,cc_protein_set=stat_num(geneAssociation.CCassociations, cc_file_name,"ppi_data")

        mf_file_name = file_name_prefix + ".f"
        mf_ppis,mf_protein_set=stat_num(geneAssociation.MFassociations, mf_file_name,"ppi_data")
        print("iea :{}".format(iea_name))
        print("data set name:{}".format(dataset_name))

        print("CC ppi num is:{}".format(len(cc_ppis)))
        print("BP ppi num is:{}".format(len(bp_ppis)))
        print("MF ppi num is:{}".format(len(mf_ppis)))

        newppis=bp_ppis+cc_ppis+mf_ppis
        totoal_ppis=[]
        for ppi in newppis:
            if [ppi[0],ppi[1]] not in totoal_ppis and  [ppi[1],ppi[0]] not in totoal_ppis:
                totoal_ppis.append(ppi)
        print("all ppi num is:{}".format(len(totoal_ppis)))

        # print("BP")
        # bp_expresssin_file = home_path + "/datasets/expression_data/expression.p"
        # bp_ppis = stat_num(geneAssociation.BPassociations, bp_expresssin_file,"expression_data")
        # print("CC")
        # cc_expresssin_file = home_path + "/datasets/expression_data/expression.c"
        # cc_ppis = stat_num(geneAssociation.CCassociations, cc_expresssin_file,"expression_data")
        # print("MF")
        # mf_expresssin_file = home_path + "/datasets/expression_data/expression.f"
        # mf_ppis = stat_num(geneAssociation.MFassociations, mf_expresssin_file,"expression_data")
        #
        # print("BP Expression num is:{}".format(len(bp_ppis)))
        # print("CC Expression num is:{}".format(len(cc_ppis)))
        # print("MF Expression num is:{}".format(len(mf_ppis)))
        #
        # newppis = bp_ppis + cc_ppis + mf_ppis
        # totoal_ppis = []
        # for ppi in newppis:
        #     if [ppi[0], ppi[1]] not in totoal_ppis and [ppi[1], ppi[0]] not in totoal_ppis:
        #         totoal_ppis.append(ppi)
        # print("total Expression num is:{}".format(len(totoal_ppis)))
        #
        #

