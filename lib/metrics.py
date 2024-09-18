# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 12:06:21 2023

@author: haglers
"""

#
from math import log10, floor

#
def _round_sig(x, sig_figs):
    if x != 0:
        y = round(x, sig_figs-int(floor(log10(abs(x))))-1)
    else:
        y = 0
    y_str = str(y)
    z = y_str.split('.')
    if len(z[0]) >= sig_figs:
        y_str = str(int(z[0]))
    elif z[0] == '0':
        num_zeros = 0
        stop_flg = False
        for i in range(len(z[1])):
            if not stop_flg:
                if z[1][i] == '0':
                    num_zeros += 1
                else:
                    stop_flg = True
        num_sig_figs = len(z[1]) - num_zeros
        if num_sig_figs < sig_figs:
            num_zeros = sig_figs - num_sig_figs
            for _ in range(num_zeros):
                y_str += '0'
    return y_str

#
class Metrics_object(object):
    
    #
    def __init__(self):
        self.none_str = ' '
        self.sig_figs = 2
    
    #
    def _get_complexity_dict(self, max_path_length, num_lineages, num_paths,
                             num_proteins, num_vertices,
                             sum_vertices_over_paths):
        undefined_str = 'U'
        complexity_dict = {}
        if num_lineages > 0:
            complexity_dict['|L|'] = num_lineages
        else:
            complexity_dict['|L|'] = self.none_str
        if num_proteins > 0:
            complexity_dict['P'] = num_proteins
        else:
            complexity_dict['P'] = self.none_str
        if num_paths > 0:
            complexity_dict['|Pi|'] = num_paths
        else:
            complexity_dict['|Pi|'] = self.none_str
        if num_vertices > 0:
            complexity_dict['|V|'] = num_vertices
        else:
            complexity_dict['|V|'] = self.none_str
        if max_path_length > 0:
            complexity_dict['M'] = max_path_length
        else:
            complexity_dict['M'] = self.none_str
        if sum_vertices_over_paths > 0:
            complexity_dict['D'] = sum_vertices_over_paths - num_paths
        else:
            complexity_dict['D'] = self.none_str
        if num_paths > 0 and sum_vertices_over_paths > 0:
            complexity_dict['C'] = sum_vertices_over_paths - num_paths + 1
        else:
            complexity_dict['C'] = self.none_str
        if num_paths > 0 and max_path_length > 0:
            complexity_dict['F'] = max_path_length + num_paths - 1
        else:
            complexity_dict['F'] = self.none_str
        if num_paths > 0 and num_vertices > 0 and sum_vertices_over_paths > 0 and max_path_length > 0:
            C = sum_vertices_over_paths - num_paths + 1
            F = max_path_length + num_paths - 1
            try:
                x = (num_vertices - F) / (C - F)
                complexity_dict['eta'] = _round_sig(x, self.sig_figs)
            except:
                complexity_dict['eta'] = undefined_str
        else:
            complexity_dict['eta'] = self.none_str
        return complexity_dict
    
    # 
    def _get_coverage_dict(self, lattice_flg, drug_uniprot_accession_list,
                           target_uniprot_accession_list):
        coverage_dict = {}
        root_pathways = self.pathway_lattice_object.get_root_pathways()
        for root_pathway in root_pathways:
            pathway_name_list = \
                self.reactome_pathways_object.map_pathway_stable_identifier_to_pathway_name(root_pathway)
            augmented_lattice = self.pathway_lattice_object.pull_augmented_lattice()
            augmented_lattice_name = 'augmented_lattice'
            paths = \
                self.pathway_lattice_object.get_paths(augmented_lattice_name, root_pathway)
            drug_path_ctr = 0
            gene_path_ctr = 0
            num_path_ctr = 0
            seen_key_list = []
            for path in paths:
                if lattice_flg == 'LEAF':
                    key = path[-1]
                elif lattice_flg == 'LINEAGE':
                    key = path[1]
                if key not in seen_key_list:
                    seen_key_list.append(key)
                    if set(drug_uniprot_accession_list) & set(augmented_lattice[key]['UNIPROT_ACCESSIONS']):
                        drug_path_ctr += 1
                    if target_uniprot_accession_list is not None:
                        if set(target_uniprot_accession_list) & set(augmented_lattice[key]['UNIPROT_ACCESSIONS']):
                            gene_path_ctr += 1
                    num_path_ctr += 1
            coverage_dict[pathway_name_list[0]] = {}
            if drug_path_ctr > 0:
                coverage = _round_sig(drug_path_ctr / num_path_ctr, self.sig_figs)
                coverage_dict[pathway_name_list[0]]['K'] = coverage
            else:
                coverage_dict[pathway_name_list[0]]['K'] = self.none_str
            if gene_path_ctr > 0:
                coverage = _round_sig(drug_path_ctr / gene_path_ctr, self.sig_figs)
                coverage_dict[pathway_name_list[0]]['k'] = coverage
            else:
                coverage_dict[pathway_name_list[0]]['k'] = self.none_str
        return coverage_dict
    
    #
    def augmented_lattice_complexities(self):
        return self.tree_complexities('augmented_lattice')
    
    #
    def crosstalk_lattice_complexities(self):
        max_path_length = \
            self.pathway_lattice_object.get_max_path_length('crosstalk_lattice',
                                                            'TOP')
        num_lineages = \
            self.pathway_lattice_object.get_num_lineages('crosstalk_lattice',
                                                         'TOP')
        num_paths = \
            self.pathway_lattice_object.get_num_paths('crosstalk_lattice',
                                                      'TOP')
        num_proteins = \
            self.pathway_lattice_object.get_num_proteins('crosstalk_lattice',
                                                         'TOP')
        num_vertices = \
            self.pathway_lattice_object.get_num_vertices('crosstalk_lattice',
                                                         'TOP')
        sum_vertices_over_paths = \
            self.pathway_lattice_object.get_sum_vertices_over_paths('crosstalk_lattice',
                                                                    'TOP')
        complexity_dict = \
            self._get_complexity_dict(max_path_length, num_lineages, num_paths,
                                      num_proteins, num_vertices, 
                                      sum_vertices_over_paths)
        return complexity_dict
    
    #
    def drug_coverage(self, drug_uniprot_accession_list,
                      target_uniprot_accession_list):
        coverage_dict = {}
        coverage_dict['LINEAGE'] = \
            self._get_coverage_dict('LINEAGE', drug_uniprot_accession_list,
                                    target_uniprot_accession_list)
        for key in coverage_dict['LINEAGE'].keys():
            if coverage_dict['LINEAGE'][key] == 0:
                coverage_dict['LINEAGE'][key] = ' '
        coverage_dict['LEAF'] = \
            self._get_coverage_dict('LEAF', drug_uniprot_accession_list,
                                    target_uniprot_accession_list)
        for key in coverage_dict['LEAF'].keys():
            if coverage_dict['LEAF'][key] == 0:
                coverage_dict['LEAF'][key] = ' '
        return coverage_dict
    
    #
    def full_lattice_complexities(self):
        return self.tree_complexities('full_lattice')
    
    #
    def load_pathway_lattice_object(self, pathway_lattice_object):
        self.pathway_lattice_object = pathway_lattice_object
        
    #
    def load_reactome_pathways_object(self, reactome_pathways_object):
        self.reactome_pathways_object = reactome_pathways_object
        
    #
    def load_targetome_evidence_object(self, targetome_evidence_object):
        self.targetome_evidence_object = targetome_evidence_object
    
    #
    def protein_lattice_complexities(self):
        return self.tree_complexities('protein_lattice')
    
    #
    def root_protein_fraction(self, uniprot_accession_list):
        lattice_name = 'augmented_lattice'
        protein_fraction_dict = {}
        root_pathways = self.pathway_lattice_object.get_root_pathways()
        for root_pathway in root_pathways:
            pathway_name_list = \
                self.reactome_pathways_object.map_pathway_stable_identifier_to_pathway_name(root_pathway)
            x = \
                self.pathway_lattice_object.get_root_protein_fraction(lattice_name,
                                                                      root_pathway,
                                                                      uniprot_accession_list)
            if x > 0:
                protein_fraction_dict[pathway_name_list[0]] = \
                    _round_sig(x, self.sig_figs)
            else:
                protein_fraction_dict[pathway_name_list[0]] = self.none_str
        return protein_fraction_dict
    
    #
    def set_sig_figs(self, sig_figs):
        self.sig_figs = sig_figs
        
    #
    def tree_complexities(self, lattice_name):
        complexity_dict = {}
        root_pathways = self.pathway_lattice_object.get_root_pathways()
        for root_pathway in root_pathways:
            pathway_name_list = \
                self.reactome_pathways_object.map_pathway_stable_identifier_to_pathway_name(root_pathway)
            max_path_length = \
                self.pathway_lattice_object.get_max_path_length(lattice_name,
                                                                root_pathway)
            num_lineages = \
                self.pathway_lattice_object.get_num_lineages(lattice_name,
                                                             root_pathway)
            num_paths = \
                self.pathway_lattice_object.get_num_paths(lattice_name,
                                                          root_pathway)
            num_proteins = \
                self.pathway_lattice_object.get_num_proteins(lattice_name,
                                                             root_pathway)
            num_vertices = \
                self.pathway_lattice_object.get_num_vertices(lattice_name,
                                                             root_pathway)
            sum_vertices_over_paths = \
                self.pathway_lattice_object.get_sum_vertices_over_paths(lattice_name,
                                                                        root_pathway)
            complexity_dict[pathway_name_list[0]] = \
                self._get_complexity_dict(max_path_length, num_lineages,
                                          num_paths, num_proteins, 
                                          num_vertices,sum_vertices_over_paths)
        return complexity_dict