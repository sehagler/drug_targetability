# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 16:11:27 2024

@author: haglers
"""

#
import os
import sys

#
path = os.path.dirname(__file__)
sys.path.insert(0, path)

#
from metrics import Metrics_object
from pathway_lattice import Pathway_lattice_object
from reactome_pathways import Reactome_pathways_object
from targetome_evidence import Targetome_evidence_object
from uniprot_to_reactome import Uniprot_to_reactome_object
from uniprot_request import Uniprot_request_object

#
class Targetability_object(object):
        
    #
    def _initialize_metrics_object(self):
        self.metrics_object = Metrics_object()

    #
    def _initialize_pathway_lattice_object(self):
        self.pathway_lattice_object = Pathway_lattice_object()
        self.pathway_lattice_object.read_reactome_pathways_relation_file(
            os.path.join(self.reactome_flat_files_path,
                         self.pathways_relation_filename))

    #
    def _initialize_reactome_pathways_object(self):
        self.reactome_pathways_object = Reactome_pathways_object()
        self.reactome_pathways_object.read_data(
            os.path.join(self.reactome_flat_files_path,
                         self.reactome_pathways_filename))

    #
    def _initialize_targetome_evidence_object(self):
        self.targetome_evidence_object = Targetome_evidence_object()
        self.targetome_evidence_object.read_data(
            os.path.join(self.targetome_flat_files_path,
                         self.targetome_evidence_filename))

    #
    def _initialize_uniprot_request_object(self):
        self.uniprot_request_object = Uniprot_request_object()

    #
    def _initialize_uniprot_to_reactome_object(self):
        self.uniprot_to_reactome_object = Uniprot_to_reactome_object()
        self.uniprot_to_reactome_object.read_data(
            os.path.join(self.reactome_flat_files_path,
                         self.all_levels_pathways_filename))
        
    #
    def _load_pathway_lattice_object(self):
        self.pathway_lattice_object.load_reactome_pathways_object(
            self.reactome_pathways_object)
        self.pathway_lattice_object.load_targetome_evidence_object(
            self.targetome_evidence_object)
        self.pathway_lattice_object.load_uniprot_to_reactome_object(
            self.uniprot_to_reactome_object)

    #
    def _load_pkl_file(self, pkl_path, pkl_filename):
        self.pathway_lattice_object.load_data_dict(
            os.path.join(pkl_path, pkl_filename))

    #
    def _load_metrics_object(self):
        self.metrics_object.load_pathway_lattice_object(
            self.pathway_lattice_object)
        self.metrics_object.load_reactome_pathways_object(
            self.reactome_pathways_object)
        self.metrics_object.load_targetome_evidence_object(
            self.targetome_evidence_object)
        
    #
    def get_drug_data_dict(self, drug_list_in, max_assay_value):
        self.targetome_evidence_object.set_max_assay_value(max_assay_value)
        drug_list = \
            self.targetome_evidence_object.get_drugs_in_targetome(drug_list_in)
        drug_uniprot_accession_list = \
            self.targetome_evidence_object.map_drug_to_uniprot_accession(drug_list)
        targeting_uniprot_accession_list = \
            drug_uniprot_accession_list
        targeting_list, no_result_list = \
            self.uniprot_request_object.map_uniprot_accession_to_gene_name(targeting_uniprot_accession_list)
        augmented_lattice_complexity_dict = \
            self.metrics_object.augmented_lattice_complexities()
        self.pathway_lattice_object.create_protein_lattice(drug_uniprot_accession_list)
        drug_lattice_complexity_dict = \
            self.metrics_object.protein_lattice_complexities()
        drug_coverage_dict = \
            self.metrics_object.drug_coverage(drug_uniprot_accession_list,
                                              None)
        protein_fraction_dict = \
            self.metrics_object.root_protein_fraction(drug_uniprot_accession_list)
        data_dict = {}
        data_dict['drug_list_in'] = drug_list_in
        data_dict['max_assay_value'] = [str(max_assay_value)]
        data_dict['drug_list'] = drug_list
        data_dict['targeting_list'] = targeting_list
        data_dict['targeting_uniprot_accession_list'] = \
            targeting_uniprot_accession_list
        data_dict['augmented_lattice_complexity_dict'] = \
            augmented_lattice_complexity_dict
        data_dict['drug_lattice_complexity_dict'] = \
            drug_lattice_complexity_dict
        data_dict['drug_coverage_dict'] = drug_coverage_dict
        data_dict['protein_fraction_dict'] = protein_fraction_dict
        return data_dict
    
    #
    def get_drug_target_data_dict(self, drug_list_in, target_list_in, max_assay_value):
        self.targetome_evidence_object.set_max_assay_value(max_assay_value)
        drug_list = \
            self.targetome_evidence_object.get_drugs_in_targetome(drug_list_in)
        drug_uniprot_accession_list = \
            self.targetome_evidence_object.map_drug_to_uniprot_accession(drug_list)
        target_uniprot_accession_list = \
            self.uniprot_request_object.map_gene_name_to_uniprot_accession(target_list_in)
        target_list, no_result_list = \
            self.uniprot_request_object.map_uniprot_accession_to_gene_name(target_uniprot_accession_list)
        targeting_uniprot_accession_list = \
            list(set(drug_uniprot_accession_list) & \
                 set(target_uniprot_accession_list))
        overtargeting_uniprot_accession_list = \
            list(set(drug_uniprot_accession_list) - \
                 set(targeting_uniprot_accession_list))
        undertargeting_uniprot_accession_list = \
            list(set(target_uniprot_accession_list) - \
                 set(targeting_uniprot_accession_list))
        targeting_list, no_result_list = \
            self.uniprot_request_object.map_uniprot_accession_to_gene_name(targeting_uniprot_accession_list)
        overtargeting_list, no_result_list = \
            self.uniprot_request_object.map_uniprot_accession_to_gene_name(overtargeting_uniprot_accession_list)
        undertargeting_list, no_result_list = \
            self.uniprot_request_object.map_uniprot_accession_to_gene_name(undertargeting_uniprot_accession_list)
        drugging_list = []
        for drug in drug_list:
            tmp_uniprot_accession_list = \
                self.targetome_evidence_object.map_drug_to_uniprot_accession(drug)
            tmp_targets = list(set(tmp_uniprot_accession_list) & \
                               set(targeting_uniprot_accession_list))
            if len(tmp_targets) > 0:
                drugging_list.append(drug)
        drug_uniprot_accession_list = \
            self.targetome_evidence_object.map_drug_to_uniprot_accession(drugging_list)
        augmented_lattice_complexity_dict = \
            self.metrics_object.augmented_lattice_complexities()
        self.pathway_lattice_object.create_protein_lattice(drug_uniprot_accession_list)
        drug_lattice_complexity_dict = \
            self.metrics_object.protein_lattice_complexities()
        drug_coverage_dict = \
            self.metrics_object.drug_coverage(drug_uniprot_accession_list,
                                              targeting_uniprot_accession_list)
        protein_fraction_dict = \
            self.metrics_object.root_protein_fraction(drug_uniprot_accession_list)
        data_dict = {}
        data_dict['drug_list_in'] = drug_list_in
        data_dict['target_list_in'] = target_list_in
        data_dict['max_assay_value'] = [str(max_assay_value)]
        data_dict['drug_list'] = drug_list
        data_dict['target_list'] = target_list
        data_dict['targeting_list'] = targeting_list
        data_dict['overtargeting_list'] = overtargeting_list
        data_dict['undertargeting_list'] = undertargeting_list
        data_dict['drugging_list'] = drugging_list
        data_dict['targeting_uniprot_accession_list'] = \
            targeting_uniprot_accession_list
        data_dict['overtargeting_uniprot_accession_list'] = \
            overtargeting_uniprot_accession_list
        data_dict['undertargeting_uniprot_accession_list'] = \
            undertargeting_uniprot_accession_list
        data_dict['augmented_lattice_complexity_dict'] = \
            augmented_lattice_complexity_dict
        data_dict['drug_lattice_complexity_dict'] = \
            drug_lattice_complexity_dict
        data_dict['drug_coverage_dict'] = drug_coverage_dict
        data_dict['protein_fraction_dict'] = protein_fraction_dict
        return data_dict
    
    #
    def get_target_data_dict(self, target_list, max_assay_value):
        self.targetome_evidence_object.set_max_assay_value(max_assay_value)
        targeting_uniprot_accession_list = \
            self.uniprot_request_object.map_gene_name_to_uniprot_accession(target_list)
        targeting_list, no_result_list = \
            self.uniprot_request_object.map_uniprot_accession_to_gene_name(targeting_uniprot_accession_list)
        augmented_lattice_complexity_dict = \
            self.metrics_object.augmented_lattice_complexities()
        self.pathway_lattice_object.create_protein_lattice(targeting_uniprot_accession_list)
        protein_lattice_complexity_dict = \
            self.metrics_object.protein_lattice_complexities()
        protein_fraction_dict = \
            self.metrics_object.root_protein_fraction(targeting_uniprot_accession_list)
        data_dict = {}
        data_dict['target_list'] = target_list
        data_dict['max_assay_value'] = [str(max_assay_value)]
        data_dict['targeting_list'] = targeting_list
        data_dict['targeting_uniprot_accession_list'] = \
            targeting_uniprot_accession_list
        data_dict['augmented_lattice_complexity_dict'] = \
            augmented_lattice_complexity_dict
        data_dict['target_lattice_complexity_dict'] = \
            protein_lattice_complexity_dict
        data_dict['protein_fraction_dict'] = protein_fraction_dict
        return data_dict
    
    #
    def initialize(self):
        self._initialize_metrics_object()
        self._initialize_pathway_lattice_object()
        self._initialize_reactome_pathways_object()
        self._initialize_targetome_evidence_object()
        self._initialize_uniprot_request_object()
        self._initialize_uniprot_to_reactome_object()
        self._load_metrics_object()
        self._load_pathway_lattice_object()
        
    #
    def load_pkl_file(self, pkl_filename):
        self._load_pkl_file(self.pkl_path, pkl_filename)

    # Generate the data set for Paper_0
    def paper_0(self):
        self.targetome_evidence_object.set_max_assay_value(100)
        
        #
        root_pathways = [ 'R-HSA-1643685', 'R-HSA-162582' ]
        drug_list = [ 'Sunitinib Malate', 'Bosutinib', 'Imatinib Mesylate' ]
        lattice = self.pathway_lattice_object._pull_lattice_by_root('augmented_lattice',
                                                                    root_pathways[0])
        uniprot_accessions_0 = lattice[root_pathways[0]]['UNIPROT_ACCESSIONS']
        lattice = self.pathway_lattice_object._pull_lattice_by_root('augmented_lattice',
                                                                    root_pathways[1])
        uniprot_accessions_1 = lattice[root_pathways[1]]['UNIPROT_ACCESSIONS']
        uniprot_accessions_2 = \
            self.targetome_evidence_object.map_drug_to_uniprot_accession([drug_list[0]])
        uniprot_accessions_3 = \
            self.targetome_evidence_object.map_drug_to_uniprot_accession([drug_list[1]])
        uniprot_accessions_4 = \
            self.targetome_evidence_object.map_drug_to_uniprot_accession([drug_list[2]])
        uniprot_accessions = \
            list(set(uniprot_accessions_0) & set(uniprot_accessions_1) & \
                 set(uniprot_accessions_2) & set(uniprot_accessions_3) & \
                 set(uniprot_accessions_4))
        uniprot_accessions = sorted(uniprot_accessions)
        protein_uniprot_accession = uniprot_accessions[1]
        
        #
        protein_list = [ protein_uniprot_accession ]
        parameters_dict = {}
        parameters_dict['drug_name_0'] = 'Sunitinib Malate'
        parameters_dict['drug_name_1'] = 'Bosutinib'
        parameters_dict['drug_name_2'] = 'Imatinib Mesylate'
        parameters_dict['gene_name'], no_results_list = \
            self.uniprot_request_object.map_uniprot_accession_to_gene_name([protein_uniprot_accession])
        parameters_dict['gene_name'] = parameters_dict['gene_name'][0]
        parameters_dict['protein_name'], no_results_list = \
            self.uniprot_request_object.map_uniprot_accession_to_protein_fullname([protein_uniprot_accession])
        parameters_dict['protein_name'] = parameters_dict['protein_name'][0]
        parameters_dict['root_pathway_0_name'] = \
            self.reactome_pathways_object.map_pathway_stable_identifier_to_pathway_name(root_pathways[0])[0]
        parameters_dict['root_pathway_1_name'] = \
            self.reactome_pathways_object.map_pathway_stable_identifier_to_pathway_name(root_pathways[1])[0]
        
        #
        drug_list = [ parameters_dict['drug_name_0'],
                      parameters_dict['drug_name_1'],
                      parameters_dict['drug_name_2'],
                      parameters_dict['drug_name_2'] ]
        
        #
        data_dict = {}
        data_dict['parameters_dict'] = parameters_dict
        
        # full lattice
        full_lattice_complexity_dict = \
            self.metrics_object.full_lattice_complexities()
        data_dict['full_lattice_complecity_dict'] = \
            full_lattice_complexity_dict
        
        #augented lattice
        augmented_lattice_complexity_dict = \
            self.metrics_object.augmented_lattice_complexities()
        data_dict['augmented_lattice_complecity_dict'] = \
            augmented_lattice_complexity_dict
        
        #crosstalk lattice
        crosstalk_lattice_complexity_dict = \
            self.metrics_object.crosstalk_lattice_complexities()
        data_dict['crosstalk_lattice_complecity_dict'] = \
            crosstalk_lattice_complexity_dict
        
        #protein lattice
        self.pathway_lattice_object.create_protein_lattice(protein_list)
        protein_lattice_complexity_dict = \
            self.metrics_object.protein_lattice_complexities()
        data_dict['protein_lattice_complecity_dict'] = \
            protein_lattice_complexity_dict
        
        #example paths
        example_tree_path_list = {}
        for root_pathway in root_pathways:
            leaves = self.pathway_lattice_object.get_leaves('protein_lattice', root_pathway)
            pathway_name_list = \
                self.reactome_pathways_object.map_pathway_stable_identifier_to_pathway_name(root_pathway)
            leaves.sort()
            example_tree_path_list[pathway_name_list[0]] = \
                self.pathway_lattice_object.get_paths_to_vertex('trees_lattice', leaves[-1])
        data_dict['example_tree_path_list'] = example_tree_path_list
        example_crosstalk_lattice_path_list = \
            self.pathway_lattice_object.get_paths_to_vertex('crosstalk_lattice', 'R-HSA-629587')
        data_dict['example_crosstalk_lattice_path_list'] = \
            example_crosstalk_lattice_path_list
        
        #drug lattices
        drug_lattice_complexity_dict_dict = {}
        drug_coverage_dict_dict = {}
        protein_fraction_dict_dict = {}
        target_uniprot_accession_list, gene_name_found_list = \
            self.uniprot_request_object.map_gene_name_to_uniprot_accession([parameters_dict['gene_name']])
        for i in range(len(drug_list)):
            if i == len(drug_list) - 1:
                self.targetome_evidence_object.set_max_assay_value(1000)
            drug_name = drug_list[i]
            drug_uniprot_accession_list = \
                self.targetome_evidence_object.map_drug_to_uniprot_accession([drug_name])
            self.pathway_lattice_object.create_protein_lattice(drug_uniprot_accession_list)
            drug_lattice_complexity_dict = \
                self.metrics_object.protein_lattice_complexities()
            drug_lattice_complexity_dict_dict[i] = drug_lattice_complexity_dict
            drug_coverage_dict = \
                self.metrics_object.drug_coverage(drug_uniprot_accession_list,
                                                  target_uniprot_accession_list)
            drug_coverage_dict_dict[i] = drug_coverage_dict
            protein_fraction_dict = \
                self.metrics_object.root_protein_fraction(
                    drug_uniprot_accession_list)
            protein_fraction_dict_dict[i] = protein_fraction_dict
        data_dict['drug_lattice_complexity_dict'] = \
            drug_lattice_complexity_dict_dict
        data_dict['drug_coverage_dict_dict'] = drug_coverage_dict_dict
        data_dict['protein_fraction_dict_dict'] = protein_fraction_dict_dict
        
        #
        return data_dict
    
    #
    def push_all_levels_pathways_filename(self, all_levels_pathways_filename):
            self.all_levels_pathways_filename = all_levels_pathways_filename
    
    #
    def push_pathways_relation_filename(self, pathways_relation_filename):
        self.pathways_relation_filename = pathways_relation_filename
        
    #
    def push_pkl_path(self, pkl_path):
        self.pkl_path = pkl_path
    
    #
    def push_reactome_flat_files_path(self, reactome_flat_files_path):
        self.reactome_flat_files_path = reactome_flat_files_path
        
    #
    def push_reactome_pathways_filename(self, reactome_pathways_filename):
        self.reactome_pathways_filename = reactome_pathways_filename
        
    #
    def push_targetome_evidence_filename(self, targetome_evidence_filename):
        self.targetome_evidence_filename = targetome_evidence_filename
        
    #
    def push_targetome_flat_files_path(self, targetome_flat_files_path):
        self.targetome_flat_files_path = targetome_flat_files_path
    
    #
    def update(self, pkl_filename):
        self.pathway_lattice_object.create_lattices()
        self.pathway_lattice_object.save_data_dict(
            os.path.join(self.pkl_path, pkl_filename))