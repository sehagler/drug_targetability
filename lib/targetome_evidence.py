# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 10:21:26 2023

@author: haglers
"""

#
import csv

#
def _create_data_list_list(file):
    with open(file, 'r') as f:
        csv_reader = csv.reader(f, delimiter="\t")
        data_list_list = []
        next(csv_reader)
        for row in csv_reader:
            data_list_list.append(row)
    data_list_list = _unique(data_list_list)
    return data_list_list

#
def _is_within_bounds(row, max_assay_value):
    assay_type = row[7]
    assay_relation = row[8]
    assay_value = row[9]
    evidence_level = row[10]
    ret_val = False
    if evidence_level == 'III' and assay_relation == '=' and assay_value.isnumeric():
        if max_assay_value == None:
            ret_val = True
        elif float(assay_value) <= max_assay_value:
            ret_val = True
    return ret_val

#
def _unique(data_list_list):
    data_list_list = [list(x) for x in set(tuple(x) for x in data_list_list)]
    return data_list_list

#
class Targetome_evidence_object(object):
    
    #
    def __init__(self):
        self.max_assay_value = None
        
    #
    def _get_targetome_drug_list(self):
        targetome_drug_list = []
        for item in self.data_list_list:
            targetome_drug_list.append(item[0].lower())
        targetome_drug_list = list(set(targetome_drug_list))
        return targetome_drug_list
    
    #
    def get_drugs_in_targetome(self, drug_list):
        targetome_drug_list = self._get_targetome_drug_list()
        drug_in_targetome_list = []
        for drug in drug_list:
            if drug.lower() in targetome_drug_list:
                drug_in_targetome_list.append(drug)
            else:
                print('Targetome query failed to find drug name ' + drug)
        return drug_in_targetome_list
    
    #
    def map_drug_to_uniprot_accession(self, drug_name_list):
        uniprot_accession_list = []
        for item in self.data_list_list:
            if item[0] in drug_name_list and item[2] == 'Protein':
                if _is_within_bounds(item, self.max_assay_value):
                    uniprot_accession_list.append(item[3])
        uniprot_accession_list = list(set(uniprot_accession_list))
        return uniprot_accession_list
    
    #
    def map_uniprot_accession_to_drug(self, uniprot_accession):
        if type(uniprot_accession) == list:
            uniprot_accession_list = uniprot_accession
        else:
            uniprot_accession_list = [ uniprot_accession ]
        drug_list = []
        for item in self.data_list_list:
            if item[3] in uniprot_accession_list:
                if _is_within_bounds(item, self.max_assay_value):
                    drug_list.append(item[0])
        drug_list = list(set(drug_list))
        return sorted(drug_list)

    #
    def read_data(self, file):
        self.data_list_list = _create_data_list_list(file)
        
    #
    def set_max_assay_value(self, max_assay_value):
        self.max_assay_value = max_assay_value