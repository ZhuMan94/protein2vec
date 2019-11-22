
from src.train_model_cpu import train_model_for_ecoli_and_yeast
from src.train_model_for_human import train_model_for_human

dataset_names=["human"]
split_modes=["connectivity"]
ontology_names=["MolecularFunction"]
iea_names=["iea-"]
for dataset_name in dataset_names:
    for split_mode in split_modes:
        for ontology_name in ontology_names:
            for iea_name in iea_names:
                    if split_mode=="fivefold":
                        part_names=["part0","part1","part2","part3","part4"]
                        for part_name in part_names:
                            if dataset_name=="human":
                                train_model_for_human(dataset_name, split_mode, ontology_name, iea_name,
                                                                part_name)
                            else:
                                train_model_for_ecoli_and_yeast(dataset_name,split_mode,ontology_name,iea_name,part_name)
                    else:
                        if dataset_name=="human":
                            train_model_for_human(dataset_name, split_mode, ontology_name, iea_name,"")
                        else:
                            train_model_for_ecoli_and_yeast(dataset_name, split_mode, ontology_name, iea_name,"")

