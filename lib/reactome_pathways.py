# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 10:21:26 2023

@author: haglers
"""

#
import csv
from utilities import _unique

#
def _create_data_list_list(file):
    with open(file, 'r') as f:
        csv_reader = csv.reader(f, delimiter="\t")
        data_list_list = []
        for row in csv_reader:
            data_list_list.append(row)
    data_list_list = _unique(data_list_list)
    return data_list_list

#
class Reactome_pathways_object(object):
    
    #
    def map_pathway_name_to_pathway_stable_identifier(self, pathway_name):
        if type(pathway_name) == list:
            pathway_name_list = pathway_name
        else:
            pathway_name_list = [ pathway_name ]
        pathway_stable_identifier_list = []
        for item0 in pathway_name_list:
            for item1 in self.data_list_list:
                if item1[1] == item0:
                    pathway_stable_identifier_list.append(item1[0])
        pathway_stable_identifier_list = \
            list(set(pathway_stable_identifier_list))
        return pathway_stable_identifier_list
    
    #
    def map_pathway_stable_identifier_to_pathway_name(self, 
                                                      pathway_stable_identifier):
        if type(pathway_stable_identifier) == list:
            pathway_stable_identifier_list = pathway_stable_identifier
        else:
            pathway_stable_identifier_list = [ pathway_stable_identifier ]
        pathway_name_list = []
        for item0 in pathway_stable_identifier_list:
            for item1 in self.data_list_list:
                if item1[0] == item0:
                    pathway_name_list.append(item1[1])
        pathway_name_list = list(set(pathway_name_list))
        return pathway_name_list
    
    #
    def name(self, key):
        pathway = key
        for item in self.data_list_list:
            if item[0] == key:
                pathway = item[1]
        return pathway

    #
    def read_data(self, file):
        self.data_list_list = _create_data_list_list(file)