# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 09:52:39 2023

@author: haglers
"""

#
from math import log10, floor

#
def _conjointness_percentage(eta0, eta1, sig):
    if eta1 == ' ':
        K = ' '
    elif eta1 == 'U':
        K = 'U'
    elif eta1 > eta0:
        eta0 = float(eta0)
        eta1 = float(eta1)
        K = (eta1 - eta0) / (1 - eta0)
        K = str(_round_sig(K, sig))
    elif eta1 == eta0:
        K = '0'
    elif eta1 < eta0:
        eta0 = float(eta0)
        eta1 = float(eta1)
        K = (eta1 - eta0) / eta0
        K = str(_round_sig(K, sig))
    return K

#
def _round_sig(x, sig):
    if x != 0:
        y = round(x, sig-int(floor(log10(abs(x))))-1)
    else:
        y = 0
    y_str = str(y)
    z = y_str.split('.')
    if len(z[0]) >= sig:
        y = int(z[0])
    return y

#
class Paper_0_output_object(object):
    
    #
    def __append_crosstalk_lattice_data(self, output_text, lattice_dict, data_prefix):
        output_text += data_prefix + '_key=Crosstalk Lattice*\n'
        output_text += data_prefix + '_|L|=' + str(lattice_dict['|L|']) + '*\n'
        output_text += data_prefix + '_|Pi|=' + str(lattice_dict['|Pi|']) + '*\n'
        output_text += data_prefix + '_|V|=' + str(lattice_dict['|V|']) + '*\n'
        output_text += data_prefix + '_C=' + str(lattice_dict['C']) + '*\n'
        output_text += data_prefix + '_D=' + str(lattice_dict['D']) + '*\n'
        output_text += data_prefix + '_F=' + str(lattice_dict['F']) + '*\n'
        output_text += data_prefix + '_M=' + str(lattice_dict['M']) + '*\n'
        output_text += data_prefix + '_P=' + str(lattice_dict['P']) + '*\n'
        output_text += data_prefix + '_eta=' + str(lattice_dict['eta']) + '*\n'
        return output_text
    
    #
    def __append_lattice_data(self, output_text, lattice_dict, data_prefix):
        keys = lattice_dict.keys()
        keys = sorted(keys)
        n = 0
        for i in range(len(keys)):
            output_text += data_prefix + '_key' + str(n) + '=' + keys[i] + '\n'
            output_text += data_prefix + '_|L|' + str(n) + '=' + \
                           str(lattice_dict[keys[i]]['|L|']) + '\n'
            output_text += data_prefix + '_|Pi|' + str(n) + '=' + \
                           str(lattice_dict[keys[i]]['|Pi|']) + '\n'
            output_text += data_prefix + '_|V|' + str(n) + '=' + \
                           str(lattice_dict[keys[i]]['|V|']) + '\n'
            output_text += data_prefix + '_C' + str(n) + '=' + \
                           str(lattice_dict[keys[i]]['C']) + '\n'
            output_text += data_prefix + '_D' + str(n) + '=' + \
                           str(lattice_dict[keys[i]]['D']) + '\n'
            output_text += data_prefix + '_F' + str(n) + '=' + \
                           str(lattice_dict[keys[i]]['F']) + '\n'
            output_text += data_prefix + '_M' + str(n) + '=' + \
                           str(lattice_dict[keys[i]]['M']) + '\n'
            output_text += data_prefix + '_P' + str(n) + '=' + \
                           str(lattice_dict[keys[i]]['P']) + '\n'
            output_text += data_prefix + '_eta' + str(n) + '=' + \
                           str(lattice_dict[keys[i]]['eta']) + '\n'
            n += 1
        return output_text
    
    #
    def _append_augmented_lattice_complexity_metrics_data(self, output_text):
        lattice_dict = self.augmented_lattice_complexity_dict
        data_prefix = 'augmented_lattice_complexity_metrics'
        output_text = \
            self.__append_lattice_data(output_text, lattice_dict, data_prefix)
        return output_text
    
    #
    def _append_calculations(self, output_text):
        keys = list(self.augmented_lattice_complexity_dict.keys())
        keys = sorted(keys)
        for n in [ 10, 26 ]:
            key = keys[n]
            if n == 10:
                data_prefix = 'calculation_disease'
            elif n == 26:
                data_prefix = 'calculation_signal_transduction'
            x = self.protein_lattice_complexity_dict[key]['|L|']
            x = x / self.augmented_lattice_complexity_dict[key]['|L|']
            x = _round_sig(100*x, 2)
            output_text += data_prefix + '_augmented_tree_protein_lineage_percentage' + \
                '=' + str(x) + '\n'
            for i in range(3):
                x = self.drug_lattice_complexity_dict_dict[i][key]['|L|']
                x = x / self.augmented_lattice_complexity_dict[key]['|L|']
                x = _round_sig(x, 2)
                output_text += data_prefix + '_augmented_tree_drug_' + str(i) + '_lineage_percentage' + \
                    '=' + str(x) + '\n'
                x = self.drug_lattice_complexity_dict_dict[i][key]['|L|']
                x = x / self.protein_lattice_complexity_dict[key]['|L|']
                x = _round_sig(x, 2)
                output_text += data_prefix + '_protein_drug_' + str(i) + '_lineage_percentage' + \
                    '=' + str(x) + '\n'
            x = self.protein_lattice_complexity_dict[key]['|Pi|']
            x = x / self.augmented_lattice_complexity_dict[key]['|Pi|']
            x = _round_sig(100*x, 2)
            output_text += data_prefix + '_augmented_tree_protein_path_percentage' + \
                '=' + str(x) + '\n'
            for i in range(3):
                x = self.drug_lattice_complexity_dict_dict[i][key]['|Pi|']
                x = x / self.augmented_lattice_complexity_dict[key]['|Pi|']
                x = _round_sig(x, 2)
                output_text += data_prefix + '_augmented_tree_drug_' + str(i) + '_path_percentage' + \
                    '=' + str(x) + '\n'
                x = self.drug_lattice_complexity_dict_dict[i][key]['|Pi|']
                x = x / self.protein_lattice_complexity_dict[key]['|Pi|']
                x = _round_sig(x, 2)
                output_text += data_prefix + '_protein_drug_' + str(i) + '_path_percentage' + \
                    '=' + str(x) + '\n'
                x = self.drug_lattice_complexity_dict_dict[i][key]['|V|']
                x = x / self.drug_lattice_complexity_dict_dict[i][key]['|Pi|']
                x = _round_sig(x, 2)
                output_text += data_prefix + '_drug_' + str(i) + '_vertices_per_path' + \
                    '=' + str(x) + '\n'
                x = self.drug_lattice_complexity_dict_dict[i][key]['D']
                x = x / self.drug_lattice_complexity_dict_dict[i][key]['|Pi|']
                x = _round_sig(x, 2)
                output_text += data_prefix + '_drug_' + str(i) + '_avg_path_depth' + \
                    '=' + str(x) + '\n'
                x = self.drug_lattice_complexity_dict_dict[i][key]['C']
                x = x / self.drug_lattice_complexity_dict_dict[i][key]['|Pi|']
                x = _round_sig(x, 2)
                output_text += data_prefix + '_drug_' + str(i) + '_ceiling_vertices_per_path' + \
                    '=' + str(x) + '\n'
                x = self.drug_lattice_complexity_dict_dict[i][key]['F']
                x = x / self.drug_lattice_complexity_dict_dict[i][key]['|Pi|']
                x = _round_sig(x, 2)
                output_text += data_prefix + '_drug_' + str(i) + '_floor_vertices_per_path' + \
                    '=' + str(x) + '\n'
                x = self.drug_lattice_complexity_dict_dict[i][key]['F']
                x = x - self.drug_lattice_complexity_dict_dict[i][key]['|Pi|'] + 1
                output_text += data_prefix + '_drug_' + str(i) + '_max_path_depth' + \
                    '=' + str(x) + '\n'
        return output_text
            
    #
    def _append_crosstalk_lattice_complexity_metrics_data(self, output_text):
        lattice_dict = self.crosstalk_lattice_complexity_dict
        data_prefix = 'crosstalk_lattice_complexity_metrics'
        output_text = \
            self.__append_crosstalk_lattice_data(output_text, lattice_dict, data_prefix)
        return output_text
    
    #
    def _append_drug_coverage_metrics_data(self, output_text):
        keys_i = list(self.drug_coverage_dict_dict.keys())
        for i in range(len(keys_i)):
            drug_coverage_dict = self.drug_coverage_dict_dict[keys_i[i]]
            keys_j = list(drug_coverage_dict['LINEAGE'].keys())
            keys_j = sorted(keys_j)
            n = 0
            for j in range(len(keys_j)):
                output_text += 'drug_' + str(keys_i[i]) + '_coverage_metrics_key' + str(n) + '=' + keys_j[j] + '\n'
                output_text += 'drug_' + str(keys_i[i]) + '_coverage_lineage_metrics_K' + str(n) + '=' + \
                               str(drug_coverage_dict['LINEAGE'][keys_j[j]]['K']) + '\n'
                output_text += 'drug_' + str(keys_i[i]) + '_coverage_lineage_metrics_k' + str(n) + '=' + \
                               str(drug_coverage_dict['LINEAGE'][keys_j[j]]['k']) + '\n'
                output_text += 'drug_' + str(keys_i[i]) + '_coverage_leaf_metrics_K' + str(n) + '=' + \
                               str(drug_coverage_dict['LEAF'][keys_j[j]]['K']) + '\n'
                output_text += 'drug_' + str(keys_i[i]) + '_coverage_leaf_metrics_k' + str(n) + '=' + \
                               str(drug_coverage_dict['LEAF'][keys_j[j]]['k']) + '\n'
                eta0 = self.augmented_lattice_complexity_dict[keys_j[j]]['eta']
                eta1 = self.drug_lattice_complexity_dict_dict[keys_i[i]][keys_j[j]]['eta']
                K = _conjointness_percentage(eta0, eta1, 2)
                output_text += 'drug_' + str(keys_i[i]) + '_coverage_metrics_K_eta' + str(n) + '=' + K + '\n'
                n += 1
        return output_text
    
    #
    def _append_drug_lattice_complexity_metrics_data(self, output_text):
        keys = list(self.drug_lattice_complexity_dict_dict.keys())
        for i in range(len(keys)):
            lattice_dict = self.drug_lattice_complexity_dict_dict[keys[i]]
            data_prefix = 'drug_' + str(keys[i]) + '_lattice_complexity_metrics'
            output_text = \
                self.__append_lattice_data(output_text, lattice_dict, data_prefix)
        return output_text
    
    #
    def _append_example_crosstalk_lattice_paths_data(self, output_text):
        n = 0
        for i in range(len(self.example_crosstalk_lattice_path_list[0])):
            output_text += 'example_crosstalk_lattice_path_id' + str(n) + '=' + \
                            self.example_crosstalk_lattice_path_list[0][i][0] + '\n'
            output_text += 'example_crosstalk_lattice_path_name' + str(n) + '=' + \
                            self.example_crosstalk_lattice_path_list[0][i][1] + '\n'
            n += 1
        return output_text
    
    #
    def _append_example_tree_paths_data(self, output_text):
        keys = self.example_tree_path_list.keys()
        keys = sorted(keys)
        for i in range(len(keys)):
            key = keys[i].replace(' ', '_')
            key = key.lower()
            n = 0
            for j in range(len(self.example_tree_path_list[keys[i]][0])):
                output_text += 'example_' + key + '_tree_path_id' + str(n) + '=' + \
                                self.example_tree_path_list[keys[i]][0][j][0] + '\n'
                output_text += 'example_' + key + '_tree_path_name' + str(n) + '=' + \
                                self.example_tree_path_list[keys[i]][0][j][1] + '\n'
                n += 1
        return output_text
    
    #
    def _append_gene_lattice_complexity_metrics_data(self, output_text):
        lattice_dict = self.gene_lattice_complexity_dict
        data_prefix = 'gene_lattice_complexity_metrics'
        output_text = \
            self.__append_lattice_data(output_text, lattice_dict, data_prefix)
        return output_text
    
    #
    def _append_lattice_complexity_metrics_data(self, output_text):
        lattice_dict = self.lattice_complexity_dict
        data_prefix = 'lattice_complexity_metrics'
        output_text = \
            self.__append_lattice_data(output_text, lattice_dict, data_prefix)
        return output_text
    
    #
    def _append_parameters(self, output_text):
        output_text += 'drug_0_name=' + self.parameters_dict['drug_name_0'] + '\n'
        output_text += 'drug_1_name=' + self.parameters_dict['drug_name_1'] + '\n'
        output_text += 'gene_name=' + self.parameters_dict['gene_name'] + '\n'
        output_text += 'protein_name=' + self.parameters_dict['protein_name'] + '\n'
        output_text += 'root_pathway_0_name=' + self.parameters_dict['root_pathway_0_name'] + '\n'
        output_text += 'root_pathway_1_name=' + self.parameters_dict['root_pathway_1_name'] + '\n'
        return output_text
    
    #
    def _append_protein_fraction_data(self, output_text):
        keys_i = list(self.protein_fraction_dict_dict.keys())
        for i in range(len(keys_i)):
            protein_fraction_dict = self.protein_fraction_dict_dict[keys_i[i]]
            keys_j = protein_fraction_dict.keys()
            keys_j = sorted(keys_j)
            n = 0
            for j in range(len(keys_j)):
                output_text += 'protein_' + str(keys_i[i]) + '_fraction_key' + str(n) + '=' + keys_j[j] + '\n'
                output_text += 'protein_' + str(keys_i[i]) + '_fraction_val' + str(n) + '=' + \
                               str(protein_fraction_dict[keys_j[j]]) + '\n'
                n += 1
        return output_text
    
    #
    def _append_protein_lattice_complexity_metrics_data(self, output_text):
        lattice_dict = self.protein_lattice_complexity_dict
        data_prefix = 'protein_lattice_complexity_metrics'
        output_text = \
            self.__append_lattice_data(output_text, lattice_dict, data_prefix)
        return output_text
    
    #
    def _append_targetome_lattice_complexity_metrics_data(self, output_text):
        lattice_dict = self.targetome_lattice_complexity_dict
        data_prefix = 'targetome_lattice_complexity_metrics'
        output_text = \
            self.__append_lattice_data(output_text, lattice_dict, data_prefix)
        return output_text
    
    #
    def _append_targetome_coverage_metrics_data(self, output_text):
        keys = self.targetome_coverage_dict['LINEAGE'].keys()
        keys = sorted(keys)
        n = 0
        for i in range(len(keys)):
            output_text += 'targetome_coverage_metrics_key' + str(n) + '=' + keys[i] + '\n'
            output_text += 'targetome_coverage_lineage_metrics_K' + str(n) + '=' + \
                           str(self.targetome_coverage_dict['LINEAGE'][keys[i]]) + '\n'
            output_text += 'targetome_coverage_leaf_metrics_K' + str(n) + '=' + \
                           str(self.targetome_coverage_dict['LEAF'][keys[i]]) + '\n'
            n += 1
        return output_text
    
    #
    def create_python_data_file(self, filename):
        output_text = ''
        output_text = self._append_parameters(output_text)
        output_text = \
            self._append_example_crosstalk_lattice_paths_data(output_text)
        output_text = self._append_example_tree_paths_data(output_text)
        output_text = \
            self._append_lattice_complexity_metrics_data(output_text)
        output_text = \
            self._append_crosstalk_lattice_complexity_metrics_data(output_text)
        output_text = \
            self._append_augmented_lattice_complexity_metrics_data(output_text)
        output_text = \
            self._append_drug_lattice_complexity_metrics_data(output_text)
        output_text = \
            self. _append_protein_lattice_complexity_metrics_data(output_text)
        output_text = \
            self._append_protein_fraction_data(output_text)
        output_text = \
            self._append_drug_coverage_metrics_data(output_text)
        output_text = \
            self._append_calculations(output_text)
        with open(filename, 'w') as f:
            f.write(output_text)
            
    #
    def push_data_dict(self, data_dict):
        self.augmented_lattice_complexity_dict = \
            data_dict['augmented_lattice_complecity_dict']
        self.crosstalk_lattice_complexity_dict = \
            data_dict['crosstalk_lattice_complecity_dict']
        self.drug_coverage_dict_dict = \
            data_dict['drug_coverage_dict_dict']
        self.drug_lattice_complexity_dict_dict = \
            data_dict['drug_lattice_complexity_dict']
        self.example_crosstalk_lattice_path_list = \
            data_dict['example_crosstalk_lattice_path_list']
        self.example_tree_path_list = \
            data_dict['example_tree_path_list']
        self.lattice_complexity_dict = \
            data_dict['full_lattice_complecity_dict']
        self.parameters_dict = data_dict['parameters_dict']
        self.protein_fraction_dict_dict = \
            data_dict['protein_fraction_dict_dict']
        self.protein_lattice_complexity_dict = \
            data_dict['protein_lattice_complecity_dict']