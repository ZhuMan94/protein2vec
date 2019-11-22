
from src.train_expression_model import train_expression_model
from src.train_sequence_model import train_sequence_model
ontology_names=["CellularComponent","MolecularFunction","BiologicalProcess"]
iea_names=["iea-"]
for ontology_name in ontology_names:
    for iea_name in iea_names:
        train_expression_model(ontology_name,iea_name)
        #train_sequence_model(ontology_name, iea_name)

