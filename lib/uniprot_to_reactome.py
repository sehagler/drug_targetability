# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 12:10:58 2023

@author: haglers
"""

#
import csv
from utilities import _unique

#
def _create_data_list_list(file):
    with open(file, 'r', encoding='utf-8') as f:
        csv_reader = csv.reader(f, delimiter="\t")
        data_list_list = []
        for row in csv_reader:
            if 'HSA' in row[1]:
                data_list_list.append(row)
    data_list_list = _unique(data_list_list)
    return data_list_list

#
def _map_uniprot_accession_list_to_pathways(data_list_list,
                                            uniprot_accession_list):
    pathway_stable_identifier_list = []
    for item0 in uniprot_accession_list:
        for item1 in data_list_list:
            if item1[0] == item0:
                pathway_stable_identifier_list.append(item1[1])
    pathway_stable_identifier_list = list(set(pathway_stable_identifier_list))
    return pathway_stable_identifier_list

#
class Uniprot_to_reactome_object(object):
    
    #
    def get_uniprot_accession_list(self):
        uniprot_accession_list = []
        for item in self.all_levels_pathways_data_list_list:
            uniprot_accession_list.append(item[0])
        uniprot_accession_list = list(set(uniprot_accession_list))
        return uniprot_accession_list
    
    #
    def map_pathway_stable_identifier_to_uniprot_accession(self,
                                                           pathway_stable_identifier):
        if type(pathway_stable_identifier) == list:
            pathway_stable_identifier_list = pathway_stable_identifier
        else:
            pathway_stable_identifier_list = [ pathway_stable_identifier ]
        uniprot_accession_list = []
        for item0 in pathway_stable_identifier_list:
            for item1 in self.all_levels_pathways_data_list_list:
                if item1[1] == item0:
                    uniprot_accession_list.append(item1[0])
        uniprot_accession_list = list(set(uniprot_accession_list))
        return uniprot_accession_list
    
    #
    def read_data(self, all_levels_pathways_file):
        self.all_levels_pathways_data_list_list = \
            _create_data_list_list(all_levels_pathways_file)