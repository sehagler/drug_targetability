# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 16:24:13 2023

@author: haglers
"""

#
import copy
import csv
import pickle

#
def _children(pathway, lattice):
    pathway_set = lattice[pathway]['SET']
    possible_children = []
    for key in lattice.keys():
        if key != pathway and \
           lattice[key]['SET'].issubset(pathway_set):
            possible_children.append(key)
    children = []
    for possible_child_0 in possible_children:
        is_child = True
        for possible_child_1 in possible_children:
            if possible_child_0 != possible_child_1 and \
               lattice[possible_child_0]['SET'].issubset(lattice[possible_child_1]['SET']):
                is_child = False
        if is_child:
            children.append(possible_child_0)
    return children

#
def _create_lattice_dict(reactome_pathways_relation):
    all_pathways = []
    parents = []
    children = []
    for item in reactome_pathways_relation:
        all_pathways.extend(item)
        parents.append(item[0])
        children.append(item[1])
    all_pathways = list(set(all_pathways))
    singletons = list(set(children) - set(parents))
    lattice_dict = {}
    for item in reactome_pathways_relation:
        if item[0] in lattice_dict.keys():
            lattice_dict[item[0]].add(item[1])
        else:
            lattice_dict[item[0]] = set(item)
    for parent in parents:
        children = lattice_dict[parent]
        for key in lattice_dict.keys():
            if parent in lattice_dict[key]:
                lattice_dict[key] |= children
    for singleton in singletons:
        lattice_dict[singleton] = set([singleton])
    lattice_dict['TOP'] = set(all_pathways)
    lattice_dict['BOTTOM'] = set()
    return lattice_dict

#
def _get_crosstalk_list(reactome_pathways_relation,
                        reactome_lattice):
    children = []
    for item in reactome_pathways_relation:
        children.append(item[1])
    unique_children = list(set(children))
    crosstalk_list = []
    for child in unique_children:
        count = children.count(child)
        if count > 1:
            crosstalk_list.append(child)
    pathway_list = copy.copy(crosstalk_list)
    while len(pathway_list) > 0:
        updated_pathway_list = []
        for pathway in pathway_list:
            children = _children(pathway, reactome_lattice)
            crosstalk_list.extend(children)
            updated_pathway_list.extend(children)
        pathway_list = copy.copy(updated_pathway_list)
    crosstalk_list = list(set(crosstalk_list))
    return crosstalk_list

#
def _get_lowest_pathways(pathways, lattice):
    lowest_pathways = []
    for pathway_0 in pathways:
        if pathway_0 != 'BOTTOM':
            is_lowest = True
            for pathway_1 in pathways:
                if pathway_1 != 'BOTTOM' and pathway_0 != pathway_1 and \
                   lattice[pathway_1]['SET'].issubset(lattice[pathway_0]['SET']):
                       is_lowest = False
            if is_lowest:
                lowest_pathways.append(pathway_0)
    return lowest_pathways

#
def _get_paths(lattice_name, pathways, lattice):
    lowest_pathways = _get_lowest_pathways(pathways, lattice)
    path_list = []
    for pathway in lowest_pathways:
        if pathway != 'BOTTOM' and \
           len(list(lattice[pathway]['SET'])) == 1:
            path_list.append(['BOTTOM', pathway])
        else:
            path_list.append([pathway])  
    run_flg = True
    while run_flg:
        updated_path_list = []
        for path in path_list:
            pathway = path[-1]
            parents = _parents(pathway, lattice)
            if len(parents) > 0:
                for parent in parents:
                    path_tmp = copy.copy(path)
                    path_tmp.append(parent)
                    updated_path_list.append(path_tmp)
            else:
                path_tmp = copy.copy(path)
                updated_path_list.append(path_tmp)
        if set(map(tuple,path_list)) == set(map(tuple,updated_path_list)):
            run_flg = False
        else:
            path_list = copy.copy(updated_path_list)
    for i in range(len(path_list)):
        path_list[i].reverse()
    if lattice_name != 'crosstalk_lattice':
        for i in range(len(path_list)):
            if 'TOP' in path_list[i]:
                path_list[i].remove('TOP')
            if 'BOTTOM' in path_list[i]:
                path_list[i].remove('BOTTOM')
    return path_list

#
def _get_paths_to_vertex(pathway_stable_identifier, lattice):
    pathway_set = lattice[pathway_stable_identifier]['SET']
    pathways = []
    for key in lattice.keys():
        if pathway_set.issubset(lattice[key]['SET']):
            pathways.append(key)
    path_list = \
        _get_paths('full_lattice', pathways, lattice)
    return path_list

#
def _parents(pathway, lattice):
    pathway_set = lattice[pathway]['SET']
    possible_parents = []
    for key in lattice.keys():
        if key != pathway and \
           pathway_set.issubset(lattice[key]['SET']):
            possible_parents.append(key)
    parents = []
    for possible_parent_0 in possible_parents:
        is_parent = True
        for possible_parent_1 in possible_parents:
            if possible_parent_0 != possible_parent_1 and \
               lattice[possible_parent_1]['SET'].issubset(lattice[possible_parent_0]['SET']):
                is_parent = False
        if is_parent:
            parents.append(possible_parent_0)
    return parents

#
def _read_reactome_pathways_relation_file(file):
    with open(file, 'r') as f:
        csv_reader = csv.reader(f, delimiter="\t")
        reactome_pthways_relation = []
        for row in csv_reader:
            if 'HSA' in row[0]:
                reactome_pthways_relation.append(row)
    reactome_pthways_relation = \
        [list(x) for x in set(tuple(x) for x in reactome_pthways_relation)]
    return reactome_pthways_relation

#
class Pathway_lattice_object(object):
    
    #
    def _create_augmented_lattice(self, lattice):
        augmented_lattice = copy.deepcopy(lattice)
        augment_idx = 0
        for key0 in lattice.keys():
            children = list(set(_children(key0, lattice)).intersection(set(lattice.keys())))
            if len(children) > 0:
                parent_uniprot_accessions = lattice[key0]['UNIPROT_ACCESSIONS']
                children_uniprot_accessions = []
                for child in children:
                    children_uniprot_accessions.extend(lattice[child]['UNIPROT_ACCESSIONS'])
                orphaned_uniprot_accessions = \
                    list(set(parent_uniprot_accessions) - set(children_uniprot_accessions))
                if len(orphaned_uniprot_accessions) > 0:
                    augment_key = 'AUG-' + str(augment_idx)
                    augment_idx += 1
                    augmented_lattice[augment_key] = {}
                    augmented_lattice[augment_key]['SET'] = { augment_key }
                    augmented_lattice[augment_key]['UNIPROT_ACCESSIONS'] = \
                        orphaned_uniprot_accessions
                    for key1 in lattice.keys():
                        if key0 in augmented_lattice[key1]['SET']:
                            augmented_lattice[key1]['SET'].add(augment_key)
        return augmented_lattice
    
    #
    def _create_reactome_lattice(self):
        lattice_dict = \
            _create_lattice_dict(self.reactome_pathways_relation)
        self.reactome_lattice = self._create_lattice(lattice_dict)
    
    #
    def _create_trees_lattice(self):
        lattice_dict = \
            _create_lattice_dict(self.reactome_pathways_relation)
        trees_lattice_dict = lattice_dict
        crosstalk_list = _get_crosstalk_list(self.reactome_pathways_relation,
                                             self.reactome_lattice)
        for crosstalk in crosstalk_list:
            if crosstalk in trees_lattice_dict.keys():
                del trees_lattice_dict[crosstalk]
        lattice = self._create_lattice(trees_lattice_dict)
        for key in lattice.keys():
            lattice[key]['SET'] = \
                lattice[key]['SET'] - set(crosstalk_list)
        self.trees_lattice = lattice
    
    #
    def _create_lattice(self, lattice_dict):
        lattice = {}
        for key in lattice_dict.keys():
            lattice[key] = {}
            if key != 'BOTTOM' and key != 'TOP':
                lattice[key]['PATHWAY_NAME'] = \
                    self.reactome_pathways_object.name(key)
                lattice[key]['UNIPROT_ACCESSIONS'] = \
                    self.uniprot_to_reactome_object.map_pathway_stable_identifier_to_uniprot_accession(key)
            else:
                lattice[key]['PATHWAY_NAME'] = None
                lattice[key]['UNIPROT_ACCESSIONS'] = []
            lattice[key]['SET'] = lattice_dict[key]
        return lattice
    
    #
    def _get_sublattice(self, lattice, pathway):
        pathway_set = lattice[pathway]['SET']
        sublattice = {}
        for key in lattice.keys():
            if lattice[key]['SET'].issubset(pathway_set):
                sublattice[key] = lattice[key]
        return sublattice
        
    #
    def _pull_lattice_by_root(self, lattice_name, root_pathway):
        if lattice_name == 'augmented_lattice':
            lattice = self.pull_augmented_full_lattice(root_pathway)
        elif lattice_name == 'crosstalk_lattice':
            lattice = self.pull_crosstalk_lattice(root_pathway)
        elif lattice_name == 'protein_lattice':
            lattice = self.pull_protein_lattice(root_pathway)
        elif lattice_name == 'full_lattice':
            lattice = self.pull_full_lattice(root_pathway)
        elif lattice_name == 'reactome_pathways_relation':
            lattice = self.reactome_pathways_relation(root_pathway)
        return lattice
    
    #
    def create_lattices(self):
        self._create_reactome_lattice()
        self._create_trees_lattice()
        
        # crosstalk lattice
        lattice_dict = \
            _create_lattice_dict(self.reactome_pathways_relation)
        crosstalk_lattice_dict = copy.deepcopy(lattice_dict)
        crosstalk_list = _get_crosstalk_list(self.reactome_pathways_relation,
                                             self.reactome_lattice)
        del_keys = []
        for key in crosstalk_lattice_dict:
            if not ( set(crosstalk_lattice_dict[key]) & set(crosstalk_list)):
                del_keys.append(key)
        for key in del_keys:
            del crosstalk_lattice_dict[key]
        self.crosstalk_lattice = self._create_lattice(crosstalk_lattice_dict)
        
        # augmented lattice
        self.augmented_full_lattice = \
            self._create_augmented_lattice(self.trees_lattice)
            
        # lineaage lattice
        self.lineage_pathway_lattice = \
            copy.deepcopy(self.trees_lattice)
        keep_pathways = [ 'TOP' ]
        root_pathways = self.get_root_pathways()
        keep_pathways.extend(root_pathways)
        for root_pathway in root_pathways:
            keep_pathways.extend(_children(root_pathway, self.trees_lattice))
        delete_pathways = list(set(self.trees_lattice.keys()) - set(keep_pathways))
        for pathway in delete_pathways:
            del self.lineage_pathway_lattice[pathway]
        for key in self.lineage_pathway_lattice.keys():
            self.lineage_pathway_lattice[key]['SET'] -= set(delete_pathways)
        
        # augmented lineage lattice
        self.augmented_lineage_lattice = \
            self._create_augmented_lattice(self.lineage_pathway_lattice)
            
   #
    def create_protein_lattice(self, uniprot_accession_list):
        uniprot_accession_set = set(uniprot_accession_list)
        self.protein_lattice = {}
        for key in self.augmented_full_lattice.keys():
            if uniprot_accession_set & set(self.augmented_full_lattice[key]['UNIPROT_ACCESSIONS']):
                self.protein_lattice[key] = self.augmented_full_lattice[key]
        
    #
    def get_drugs(self, key):
        uniprot_accessions = \
            self.reactome_lattice[key]['UNIPROT_ACCESSIONS']
        drug_list = \
            self.targetome_evidence_object.map_uniprot_accession_to_drug(uniprot_accessions)
        return drug_list
    
    #
    def get_leaves(self, lattice_name, root_pathway):
        lattice = self._pull_lattice_by_root(lattice_name, root_pathway)
        paths = _get_paths(lattice_name, list(lattice.keys()), lattice)
        leaves = []
        for path in paths:
            leaves.append(path[-1])
        return leaves
    
    #
    def get_max_path_length(self, lattice_name, root_pathway):
        paths = self.get_paths(lattice_name, root_pathway)
        max_path_length = 0
        for path in paths:
            if len(path) > max_path_length:
                max_path_length = len(path)
        return max_path_length
    
    #
    def get_num_lineages(self, lattice_name, root_pathway):
        lattice = self._pull_lattice_by_root(lattice_name, root_pathway)
        if root_pathway in lattice.keys():
            children = _children(root_pathway, lattice)
        else:
            children = []
        return len(children)
    
    #
    def get_num_paths(self, lattice_name, root_pathway):
        lattice = self._pull_lattice_by_root(lattice_name, root_pathway)
        paths = _get_paths(lattice_name, list(lattice.keys()), lattice)
        return len(paths)
    
    #
    def get_num_proteins(self, lattice_name, root_pathway):
        lattice = self._pull_lattice_by_root(lattice_name, root_pathway)
        try:
            root_uniprot_accession_list = \
                lattice[root_pathway]['UNIPROT_ACCESSIONS']
        except:
            root_uniprot_accession_list = []
        return len(root_uniprot_accession_list)
    
    #
    def get_num_vertices(self, lattice_name, root_pathway):
        lattice = self._pull_lattice_by_root(lattice_name, root_pathway)
        paths = _get_paths(lattice_name, list(lattice.keys()), lattice)
        for i in range(len(paths)):
            if paths[i][0] == 'TOP':
                paths[i] = paths[i][1:]
            if paths[i][-1] == 'BOTTOM':
                paths[i] = paths[i][:-1]
        vertices = []
        for path in paths:
            vertices.extend(path)
        vertices = list(set(vertices))
        return len(vertices)
    
    #
    def get_paths(self, lattice_name, root_pathway):
        lattice = self._pull_lattice_by_root(lattice_name, root_pathway)
        paths = _get_paths(lattice_name, list(lattice.keys()), lattice)
        return paths
    
    #
    def get_paths_from_vertex(self, pathway_stable_identifier):
        pathway_set = self.trees_lattice[pathway_stable_identifier]['SET']
        pathways = []
        for key in self.trees_lattice.keys():
            if len(self.trees_lattice[key]['SET']) < len(pathway_set) and \
               len(self.trees_lattice[key]['SET']) > 0 and \
               self.trees_lattice[key]['SET'].issubset(pathway_set):
                pathways.append(key)
        path_list = \
            _get_paths('full_lattice', pathways, self.trees_lattice)
        return path_list
    
    #
    def get_paths_to_vertex(self, lattice_name, pathway_stable_identifier):
        if lattice_name == 'crosstalk_lattice':
            path_list_tmp = \
                _get_paths_to_vertex(pathway_stable_identifier, self.crosstalk_lattice)
        elif lattice_name == 'trees_lattice':
            path_list_tmp = \
                _get_paths_to_vertex(pathway_stable_identifier, self.trees_lattice)
        path_list = []
        for i in range(len(path_list_tmp)):
            path = []
            for j in range(len(path_list_tmp[i])):
                pathway_name = \
                    self.get_pathway_name(path_list_tmp[i][j])
                path.append([ path_list_tmp[i][j], pathway_name ])
            path_list.append(path)
        return path_list
    
    #
    def get_pathway_name(self, key):
        if key[:3] != 'AUG':
            return self.reactome_lattice[key]['PATHWAY_NAME']
        else:
            return key.upper()
    
    #
    def get_protein_paths(self, uniprot_accession):
        pathways = []
        for key in self.trees_lattice.keys():
            if uniprot_accession in self.trees_lattice[key]['UNIPROT_ACCESSIONS']:
                pathways.append(key)
        path_list = \
            _get_paths('full_lattice', pathways, self.trees_lattice)
        return path_list
    
    #
    def get_root_pathways(self):
        root_pathways = _children('TOP', self.trees_lattice)
        return root_pathways
    
    #
    def get_root_protein_fraction(self, lattice_name, root_pathway, 
                                  uniprot_accession_list):
        lattice = self._pull_lattice_by_root(lattice_name, root_pathway)
        root_uniprot_accession_list = \
            lattice[root_pathway]['UNIPROT_ACCESSIONS']
        common_uniprot_accession_list = \
            list(set(root_uniprot_accession_list) & set(uniprot_accession_list))
        fraction = len(common_uniprot_accession_list)
        fraction = fraction / len(root_uniprot_accession_list)
        return fraction
    
    #
    def get_sum_vertices_over_paths(self, lattice_name, root_pathway):
        lattice = self._pull_lattice_by_root(lattice_name, root_pathway)
        paths = _get_paths(lattice_name, list(lattice.keys()), lattice)
        for i in range(len(paths)):
            if paths[i][0] == 'TOP':
                paths[i] = paths[i][1:]
            if paths[i][-1] == 'BOTTOM':
                paths[i] = paths[i][:-1]
        vertices = []
        for path in paths:
            vertices.extend(path)
        return len(vertices)
    
    #
    def load_data_dict(self, file):
        with open(file, 'rb') as f:
            data_dict = pickle.load(f)
        self.augmented_full_lattice = data_dict['augmented_full_lattice']
        self.augmented_lineage_lattice = data_dict['augmented_lineage_lattice']
        self.crosstalk_lattice = data_dict['crosstalk_lattice']
        self.reactome_lattice = data_dict['reactome_lattice']
        self.trees_lattice = data_dict['trees_lattice']

    #
    def load_reactome_pathways_object(self, reactome_pathways_object):
        self.reactome_pathways_object = reactome_pathways_object
        
    #
    def load_targetome_evidence_object(self, targetome_evidence_object):
        self.targetome_evidence_object = targetome_evidence_object
        
    #
    def load_uniprot_to_reactome_object(self, uniprot_to_reactome_object):
        self.uniprot_to_reactome_object = uniprot_to_reactome_object
    
    #
    def pull_augmented_full_lattice(self, pathway):
        lattice = self.augmented_full_lattice
        sublattice = {}
        if pathway in lattice.keys():
            pathway_set = lattice[pathway]['SET']
            for key in lattice.keys():
                if lattice[key]['SET'].issubset(pathway_set):
                    sublattice[key] = lattice[key]
        return sublattice
    
    #
    def pull_augmented_lattice(self):
        lattice = self.augmented_full_lattice
        return lattice
    
    #
    def pull_crosstalk_lattice(self, pathway):
        lattice = self.crosstalk_lattice
        sublattice = {}
        if pathway in lattice.keys():
            pathway_set = lattice[pathway]['SET']
            for key in lattice.keys():
                if lattice[key]['SET'].issubset(pathway_set):
                    sublattice[key] = lattice[key]
        return sublattice

    #
    def pull_full_lattice(self, pathway):
        lattice = self.trees_lattice
        sublattice = {}
        if pathway in lattice.keys():
            pathway_set = lattice[pathway]['SET']
            for key in lattice.keys():
                if lattice[key]['SET'].issubset(pathway_set):
                    sublattice[key] = lattice[key]
        return sublattice
    
    #
    def pull_protein_lattice(self, pathway):
        lattice = self.protein_lattice
        sublattice = {}
        if pathway in lattice.keys():
            pathway_set = lattice[pathway]['SET']
            for key in lattice.keys():
                if lattice[key]['SET'].issubset(pathway_set):
                    sublattice[key] = lattice[key]
        return sublattice
    
    #
    def read_reactome_pathways_relation_file(self, file):
        self.reactome_pathways_relation = \
            _read_reactome_pathways_relation_file(file)
            
    #
    def save_data_dict(self, file):
        data_dict = {}
        data_dict['augmented_full_lattice'] = self.augmented_full_lattice
        data_dict['augmented_lineage_lattice'] = \
            self.augmented_lineage_lattice
        data_dict['crosstalk_lattice'] = self.crosstalk_lattice
        data_dict['reactome_lattice'] = self.reactome_lattice
        data_dict['trees_lattice'] = self.trees_lattice
        with open(file, 'wb') as f:
            pickle.dump(data_dict, f)