# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 10:21:24 2023

@author: haglers
"""

#
import re

#
class Output_object(object):
    
    #
    def _generate_drug_coverage_csv(self):
        output_text = ''
        output_text += 'INPUT_DRUGS,'
        for i in range(len(self.data_dict['drug_list_in'])):
            output_text += self.data_dict['drug_list_in'][i]
        output_text += '\n'
        output_text += 'INPUT_MAX_ASSAY_VALUE,'
        for i in range(len(self.data_dict['max_assay_value'])):
            output_text += self.data_dict['max_assay_value'][i]
        output_text += '\n'
        output_text += 'INPUT_DRUGS_IN_TARGETOME,' 
        for i in range(len(self.data_dict['drug_list'])):
            output_text += self.data_dict['drug_list'][i]
        output_text += '\n'
        output_text += 'TARGETS_OF_INPUT_DRUGS_IN_TARGETOME,'
        for i in range(len(self.data_dict['targeting_list'])):
            output_text += self.data_dict['targeting_list'][i]
        output_text += '\n'
        output_text += 'TARGETS_OF_INPUT_DRUGS_IN_TARGETOME_UNIPROT_ACCESSION_IDS,'
        for i in range(len(self.data_dict['targeting_uniprot_accession_list'])):
            output_text += self.data_dict['targeting_uniprot_accession_list'][i]
        output_text += '\n'
        output_text += \
            self._write_coverage_csv(
                'DRUG_COVERAGE',
                self.data_dict['augmented_lattice_complexity_dict'],
                self.data_dict['drug_lattice_complexity_dict'],
                self.data_dict['drug_coverage_dict'],
                self.data_dict['protein_fraction_dict'])
        return output_text
    
    #
    def _generate_drug_coverage_report(self):
        output_text = ''
        output_text += 'INPUT_DRUGS\n\n' 
        output_text += self._write_text_list(self.data_dict['drug_list_in'])
        output_text += '\n'
        output_text += 'INPUT_MAX_ASSAY_VALUE\n\n' 
        output_text += self._write_text_list(self.data_dict['max_assay_value'])
        output_text += '\n'
        output_text += 'INPUT_DRUGS_IN_TARGETOME\n\n' 
        output_text += self._write_text_list(self.data_dict['drug_list'])
        output_text += '\n'
        output_text += 'TARGETS_OF_INPUT_DRUGS_IN_TARGETOME\n\n' 
        output_text += self._write_text_list(self.data_dict['targeting_list'])
        output_text += '\n'
        output_text += 'TARGETS_OF_INPUT_DRUGS_IN_TARGETOME_UNIPROT_ACCESSION_IDS\n\n' 
        output_text += \
            self._write_text_list(self.data_dict['targeting_uniprot_accession_list'])
        output_text += '\n'
        output_text += 'DRUG_COVERAGE\n\n' 
        output_text += \
            self._write_coverage_table(
                self.data_dict['augmented_lattice_complexity_dict'],
                self.data_dict['drug_lattice_complexity_dict'],
                self.data_dict['drug_coverage_dict'],
                self.data_dict['protein_fraction_dict'])
        return output_text
    
    #
    def _generate_drug_target_coverage_csv(self):
        output_text = ''
        output_text += 'INPUT_DRUGS,' 
        for i in range(len(self.data_dict['drug_list_in'])):
            output_text += self.data_dict['drug_list_in'][i]
        output_text += '\n'
        output_text += 'INPUT_TARGETS,' 
        for i in range(len(self.data_dict['target_list_in'])):
            output_text += self.data_dict['target_list_in'][i]
        output_text += '\n'
        output_text += 'INPUT_MAX_ASSAY_VALUE,' 
        for i in range(len(self.data_dict['max_assay_value'])):
            output_text += self.data_dict['max_assay_value'][i]
        output_text += '\n'
        output_text += 'INPUT_DRUGS_IN_TARGETOME,'
        for i in range(len(self.data_dict['drug_list'])):
            output_text += self.data_dict['drug_list'][i]
        output_text += '\n'
        output_text += 'INPUT_TARGETS_IN_UNIPROT,' 
        for i in range(len(self.data_dict['target_list'])):
            output_text += self.data_dict['target_list'][i]
        output_text += '\n'
        output_text += 'INPUT_DRUGS_TARGETING_INPUT_TARGETS,'
        for i in range(len(self.data_dict['drugging_list'])):
            output_text += self.data_dict['drugging_list'][i]
        output_text += '\n'
        output_text += 'INPUT_TARGETS_TARGETED_BY_INPUT_DRUGS,'
        for i in range(len(self.data_dict['targeting_list'])):
            output_text += self.data_dict['targeting_list'][i]
        output_text += '\n'
        output_text += 'ADDITIONAL_TARGETS_NOT_IN_INPUT_TARGETS,'
        for i in range(len(self.data_dict['overtargeting_list'])):
            output_text += self.data_dict['overtargeting_list'][i]
        output_text += '\n'
        output_text += 'INPUT_TARGETS_TARGETED_BY_INPUT_DRUGS_UNIPROT_ACCESSION_IDS,'
        for i in range(len(self.data_dict['targeting_uniprot_accession_list'])):
            output_text += self.data_dict['targeting_uniprot_accession_list'][i]
        output_text += '\n'
        output_text += 'ADDITIONAL_TARGETS_NOT_IN_INPUT_TARGETS_UNIPROT_ACCESSION_IDS,' 
        for i in range(len(self.data_dict['overtargeting_uniprot_accession_list'])):
            output_text += self.data_dict['overtargeting_uniprot_accession_list'][i]
        output_text += '\n'
        output_text += \
            self._write_coverage_csv(
                'DRUG_TARGET_COVERAGE',
                self.data_dict['augmented_lattice_complexity_dict'],
                self.data_dict['drug_lattice_complexity_dict'],
                self.data_dict['drug_coverage_dict'],
                self.data_dict['protein_fraction_dict'])
        return output_text
    
    #
    def _generate_drug_target_coverage_report(self):
        output_text = ''
        output_text += 'INPUT_DRUGS\n\n' 
        output_text += self._write_text_list(self.data_dict['drug_list_in'])
        output_text += '\n'
        output_text += 'INPUT_TARGETS\n\n' 
        output_text += self._write_text_list(self.data_dict['target_list_in'])
        output_text += '\n'
        output_text += 'INPUT_MAX_ASSAY_VALUE\n\n' 
        output_text += self._write_text_list(self.data_dict['max_assay_value'])
        output_text += '\n'
        output_text += 'INPUT_DRUGS_IN_TARGETOME\n\n' 
        output_text += self._write_text_list(self.data_dict['drug_list'])
        output_text += '\n'
        output_text += 'INPUT_TARGETS_IN_UNIPROT\n\n' 
        output_text += self._write_text_list(self.data_dict['target_list'])
        output_text += '\n'
        output_text += 'INPUT_DRUGS_TARGETING_INPUT_TARGETS\n\n' 
        output_text += self._write_text_list(self.data_dict['drugging_list'])
        output_text += '\n'
        output_text += 'INPUT_TARGETS_TARGETED_BY_INPUT_DRUGS\n\n' 
        output_text += self._write_text_list(self.data_dict['targeting_list'])
        output_text += '\n'
        output_text += 'ADDITIONAL_TARGETS_NOT_IN_INPUT_TARGETS\n\n' 
        output_text += \
            self._write_text_list(self.data_dict['overtargeting_list'])
        output_text += '\n'
        output_text += 'INPUT_TARGETS_TARGETED_BY_INPUT_DRUGS_UNIPROT_ACCESSION_IDS\n\n' 
        output_text += \
            self._write_text_list(self.data_dict['targeting_uniprot_accession_list'])
        output_text += '\n'
        output_text += 'ADDITIONAL_TARGETS_NOT_IN_INPUT_TARGETS_UNIPROT_ACCESSION_IDS\n\n' 
        output_text += \
            self._write_text_list(self.data_dict['overtargeting_uniprot_accession_list'])
        output_text += '\n'
        output_text += 'DRUG_TARGET_COVERAGE\n\n' 
        output_text += \
            self._write_coverage_table(
                self.data_dict['augmented_lattice_complexity_dict'],
                self.data_dict['drug_lattice_complexity_dict'],
                self.data_dict['drug_coverage_dict'],
                self.data_dict['protein_fraction_dict'])
        return output_text
    
    #
    def _generate_drug_tree_csv(self):
        output_text = ''
        output_text += 'INPUT_DRUGS,' 
        for i in range(len(self.data_dict['drug_list_in'])):
            output_text += self.data_dict['drug_list_in'][i]
        output_text += '\n'
        output_text += 'INPUT_MAX_ASSAY_VALUE,' 
        for i in range(len(self.data_dict['max_assay_value'])):
            output_text += self.data_dict['max_assay_value'][i]
        output_text += '\n'
        output_text += 'INPUT_DRUGS_IN_TARGETOME,'
        for i in range(len(self.data_dict['drug_list'])):
            output_text += self.data_dict['drug_list'][i]
        output_text += '\n'
        output_text += 'TARGETS_OF_INPUT_DRUGS_IN_TARGETOME,' 
        for i in range(len(self.data_dict['targeting_list'])):
            output_text += self.data_dict['targeting_list'][i]
        output_text += '\n'
        output_text += 'TARGETS_OF_INPUT_DRUGS_IN_TARGETOME_UNIPROT_ACCESSION_IDS,' 
        for i in range(len(self.data_dict['targeting_uniprot_accession_list'])):
            output_text += self.data_dict['targeting_uniprot_accession_list'][i]
        output_text += '\n'
        output_text += \
            self._write_drug_tree_csv(
                self.data_dict['augmented_lattice_complexity_dict'],
                self.data_dict['drug_lattice_complexity_dict'])
        return output_text
    
    #
    def _generate_drug_tree_report(self):
        output_text = ''
        output_text += 'INPUT_DRUGS\n\n' 
        output_text += self._write_text_list(self.data_dict['drug_list_in'])
        output_text += '\n'
        output_text += 'INPUT_MAX_ASSAY_VALUE\n\n' 
        output_text += self._write_text_list(self.data_dict['max_assay_value'])
        output_text += '\n'
        output_text += 'INPUT_DRUGS_IN_TARGETOME\n\n' 
        output_text += self._write_text_list(self.data_dict['drug_list'])
        output_text += '\n'
        output_text += 'TARGETS_OF_INPUT_DRUGS_IN_TARGETOME\n\n' 
        output_text += self._write_text_list(self.data_dict['targeting_list'])
        output_text += '\n'
        output_text += 'TARGETS_OF_INPUT_DRUGS_IN_TARGETOME_UNIPROT_ACCESSION_IDS\n\n' 
        output_text += \
            self._write_text_list(self.data_dict['targeting_uniprot_accession_list'])
        output_text += '\n'
        output_text += 'DRUG_TREE\n\n' 
        output_text += \
            self._write_drug_tree_table(
                self.data_dict['augmented_lattice_complexity_dict'],
                self.data_dict['drug_lattice_complexity_dict'])
        return output_text
    
    #
    def _generate_target_tree_csv(self):
        output_text = ''
        output_text += 'INPUT_TARGETS,' 
        for i in range(len(self.data_dict['target_list'])):
            output_text += self.data_dict['target_list'][i]
        output_text += '\n'
        output_text += 'INPUT_MAX_ASSAY_VALUE,'
        for i in range(len(self.data_dict['max_assay_value'])):
            output_text += self.data_dict['max_assay_value'][i]
        output_text += '\n'
        output_text += 'INPUT_TARGETS_IN_UNIPROT,'
        for i in range(len(self.data_dict['targeting_list'])):
            output_text += self.data_dict['targeting_list'][i]
        output_text += '\n'
        output_text += 'INPUT_TARGETS_IN_UNIPROT_UNIPROT_ACCESSION_IDS,'
        for i in range(len(self.data_dict['targeting_uniprot_accession_list'])):
            output_text += self.data_dict['targeting_uniprot_accession_list'][i]
        output_text += '\n'
        output_text += \
            self._write_target_tree_csv(
                self.data_dict['augmented_lattice_complexity_dict'],
                self.data_dict['target_lattice_complexity_dict'])
        return output_text
    
    #
    def _generate_target_tree_report(self):
        output_text = ''
        output_text += 'INPUT_TARGETS\n\n' 
        output_text += self._write_text_list(self.data_dict['target_list'])
        output_text += '\n'
        output_text += 'INPUT_MAX_ASSAY_VALUE\n\n' 
        output_text += self._write_text_list(self.data_dict['max_assay_value'])
        output_text += '\n'
        output_text += 'INPUT_TARGETS_IN_UNIPROT\n\n' 
        output_text += \
            self._write_text_list(self.data_dict['targeting_list'])
        output_text += '\n'
        output_text += 'INPUT_TARGETS_IN_UNIPROT_UNIPROT_ACCESSION_IDS\n\n' 
        output_text += \
            self._write_text_list(self.data_dict['targeting_uniprot_accession_list'])
        output_text += '\n'
        output_text += 'TARGET_TREE\n\n' 
        output_text += \
            self._write_target_tree_table(
                self.data_dict['augmented_lattice_complexity_dict'],
                self.data_dict['target_lattice_complexity_dict'])
        return output_text
    
    #
    def _get_max_length_complexity(self, tree_dict, key):
        len_list = [len(key)]
        keys = tree_dict.keys()
        for k in keys:
            len_list.append(len(str(tree_dict[k][key])))
        return max(len_list)
    
    #
    def _get_max_length_coverage(self, tree_dict, key1, key2):
        len_list = [len(key1)]
        keys = tree_dict.keys()
        for k in keys:
            len_list.append(len(str(tree_dict[k][key2])))
        return max(len_list)
    
    #
    def _get_max_length_fraction(self, tree_dict):
        len_list = [len('KP')]
        keys = tree_dict.keys()
        for k in keys:
            len_list.append(len(str(tree_dict[k])))
        return max(len_list)
    
    #
    def _write_coverage_column_labels(self, augmented_complexity_dict,
                                      drug_complexity_dict,
                                      coverage_dict, fraction_dict):
        output_text = ''
        keys = augmented_complexity_dict.keys()
        maxlen = len(max(keys, key=len))
        gap = ' '*(maxlen + 3)
        output_text += gap + 'P'
        maxlen = self._get_max_length_complexity(augmented_complexity_dict, 'P')
        gap = ' '*(maxlen - len('P') + 1)
        output_text += gap + '|L|'
        maxlen = self._get_max_length_complexity(augmented_complexity_dict, '|L|')
        gap = ' '*(maxlen - len('|L|') + 1)
        output_text += gap + '|Pi|'
        maxlen = self._get_max_length_complexity(augmented_complexity_dict, '|Pi|')
        gap = ' '*(maxlen - len('|Pi|') + 1)
        output_text += gap + 'eta'
        maxlen = self._get_max_length_complexity(augmented_complexity_dict, 'eta')
        gap = ' '*(maxlen - len('eta') + 1)
        output_text += gap + '| KP'
        maxlen = self._get_max_length_fraction(fraction_dict)
        gap = ' '*(maxlen - len('KP') + 1)
        output_text += gap + 'K|L|'
        maxlen = self._get_max_length_coverage(coverage_dict['LINEAGE'], 'K|L|', 'K')
        gap = ' '*(maxlen - len('k|L|') + 1)
        output_text += gap + 'k|L|'
        maxlen = self._get_max_length_coverage(coverage_dict['LEAF'], 'k|L|', 'k')
        gap = ' '*(maxlen - len('k|L|') + 1)
        output_text += gap + 'K|Pi|'
        maxlen = self._get_max_length_coverage(coverage_dict['LEAF'], 'K|Pi|', 'K')
        gap = ' '*(maxlen - len('K|Pi|') + 1)
        output_text += gap + 'k|Pi|'
        maxlen = self._get_max_length_coverage(coverage_dict['LEAF'], 'k|Pi|', 'k')
        gap = ' '*(maxlen - len('k|Pi|') + 1)
        output_text += gap + 'eta'
        output_text += '\n'
        return output_text
    
    #
    def _write_coverage_csv(self, label, augmented_complexity_dict, 
                            drug_complexity_dict, coverage_dict, 
                            fraction_dict):
        output_text = ''
        output_text_tmp = \
            self._write_coverage_column_labels(augmented_complexity_dict,
                                               drug_complexity_dict,
                                               coverage_dict,
                                               fraction_dict)
        output_text_tmp = re.sub(' +', ',', output_text_tmp)
        output_text_tmp = re.sub(',\|,', ',', output_text_tmp)
        output_text += label + ',TREE_NAME' + output_text_tmp
        keys = augmented_complexity_dict.keys()
        keys = sorted(keys)
        for i in range(len(keys)):
            output_text += label + ',' + keys[i] + ','
            output_text += str(augmented_complexity_dict[keys[i]]['P']) + ','
            output_text += str(augmented_complexity_dict[keys[i]]['|L|']) + ','
            output_text += str(augmented_complexity_dict[keys[i]]['|Pi|']) + ','
            output_text += str(augmented_complexity_dict[keys[i]]['eta']) + ','
            output_text += str(fraction_dict[keys[i]]) + ','
            output_text += str(coverage_dict['LINEAGE'][keys[i]]['K']) + ','
            output_text += str(coverage_dict['LINEAGE'][keys[i]]['k']) + ','
            output_text += str(coverage_dict['LEAF'][keys[i]]['K']) + ','
            output_text += str(coverage_dict['LEAF'][keys[i]]['k']) + ','
            output_text += str(drug_complexity_dict[keys[i]]['eta'])
            output_text += '\n'
            output_text = re.sub(', +', ',', output_text)
        return output_text
    
    #
    def _write_coverage_table(self, augmented_complexity_dict, 
                              drug_complexity_dict, coverage_dict, 
                              fraction_dict):
        output_text = ''
        output_text += \
            self._write_coverage_column_labels(augmented_complexity_dict,
                                               drug_complexity_dict,
                                               coverage_dict,
                                               fraction_dict)
        keys = augmented_complexity_dict.keys()
        keys = sorted(keys)
        maxlen = len(max(keys, key=len))
        for i in range(len(keys)):
            gap = '.' * (maxlen - len(keys[i]) + 2)
            output_text += ' ' + keys[i] + gap
            output_text += self._write_drug_element_complexity(
                augmented_complexity_dict, keys[i], 'P')
            output_text += self._write_drug_element_complexity(
                augmented_complexity_dict, keys[i], '|L|')
            output_text += self._write_drug_element_complexity(
                augmented_complexity_dict, keys[i], '|Pi|')
            output_text += self._write_drug_element_complexity(
                augmented_complexity_dict, keys[i], 'eta')
            output_text += '| '
            output_text += self._write_drug_element_fraction(
                fraction_dict, keys[i], 'KP')
            output_text += self._write_drug_element_coverage(
                coverage_dict['LINEAGE'], keys[i], 'K|L|', 'K')
            output_text += self._write_drug_element_coverage(
                coverage_dict['LINEAGE'], keys[i], 'k|L|', 'k')
            output_text += self._write_drug_element_coverage(
                coverage_dict['LEAF'], keys[i], 'K|L|', 'K')
            output_text += self._write_drug_element_coverage(
                coverage_dict['LEAF'], keys[i], 'k|L|', 'k')
            output_text += self._write_drug_element_coverage_eta(
                drug_complexity_dict, keys[i], 'eta')
            output_text += '\n'
        return output_text
    
    #
    def _write_drug_element_complexity(self, tree_dict, key0, key1):
        output_text = ''
        maxlen = self._get_max_length_complexity(tree_dict, key1)
        x = str(tree_dict[key0][key1])
        gap = ' '*(maxlen - len(x) + 1)
        output_text += x + gap
        return output_text
    
    #
    def _write_drug_element_coverage(self, tree_dict, key0, key1, key2):
        output_text = ''
        maxlen = self._get_max_length_coverage(tree_dict, key1, key2)
        x = str(tree_dict[key0][key2])
        gap = ' '*(maxlen - len(x) + 1)
        output_text += x + gap
        return output_text
    
    #
    def _write_drug_element_coverage_eta(self, tree_dict, key0, key1):
        output_text = ''
        #maxlen = self._get_max_length_complexity(tree_dict, key1)
        x = str(tree_dict[key0][key1])
        #gap = ' '*(maxlen - len(x) + 2)
        output_text += ' ' + x
        return output_text
    
    #
    def _write_drug_element_fraction(self, tree_dict, key0, key1):
        output_text = ''
        maxlen = self._get_max_length_fraction(tree_dict)
        x = str(tree_dict[key0])
        gap = ' '*(maxlen - len(x) + 1)
        output_text += x + gap
        return output_text
    
    #
    def _write_drug_tree_column_labels(self, augmented_complexity_dict,
                                       drug_complexity_dict):
        output_text = ''
        keys = augmented_complexity_dict.keys()
        maxlen = len(max(keys, key=len))
        gap = ' '*(maxlen + 3)
        output_text += gap + 'P'
        maxlen = self._get_max_length_complexity(augmented_complexity_dict, 'P')
        gap = ' '*(maxlen - len('P') + 1)
        output_text += gap + '|L|'
        maxlen = self._get_max_length_complexity(augmented_complexity_dict, '|L|')
        gap = ' '*(maxlen - len('|L|') + 1)
        output_text += gap + '|Pi|'
        maxlen = self._get_max_length_complexity(augmented_complexity_dict, '|Pi|')
        gap = ' '*(maxlen - len('|Pi|') + 1)
        output_text += gap + 'eta'
        maxlen = self._get_max_length_complexity(augmented_complexity_dict, 'eta')
        gap = ' '*(maxlen - len('eta') + 1)
        output_text += gap + '| |L|'
        maxlen = self._get_max_length_complexity(drug_complexity_dict, '|L|')
        gap = ' '*(maxlen - len('|L|') + 1)
        output_text += gap + '|Pi|'
        maxlen = self._get_max_length_complexity(drug_complexity_dict, '|Pi|')
        gap = ' '*(maxlen - len('|Pi|') + 1)
        output_text += gap + 'eta'
        output_text += '\n'
        return output_text
    
    #
    def _write_drug_tree_csv(self, augmented_complexity_dict, 
                             drug_complexity_dict):
        output_text = ''
        output_text_tmp = \
            self._write_drug_tree_column_labels(augmented_complexity_dict,
                                                drug_complexity_dict)
        output_text_tmp = re.sub(' +', ',', output_text_tmp)
        output_text_tmp = re.sub(',\|,', ',', output_text_tmp)
        output_text += 'DRUG_TREE,TREE_NAME' + output_text_tmp
        keys = augmented_complexity_dict.keys()
        keys = sorted(keys)
        for i in range(len(keys)):
            output_text += 'DRUG_TREE,' + keys[i] + ','
            output_text += str(augmented_complexity_dict[keys[i]]['P']) + ','
            output_text += str(augmented_complexity_dict[keys[i]]['|L|']) + ','
            output_text += str(augmented_complexity_dict[keys[i]]['|Pi|']) + ','
            output_text += str(augmented_complexity_dict[keys[i]]['eta']) + ','
            output_text += str(drug_complexity_dict[keys[i]]['|L|']) + ','
            output_text += str(drug_complexity_dict[keys[i]]['|Pi|']) + ','
            output_text += str(drug_complexity_dict[keys[i]]['eta']) + ','
            output_text += '\n'
            output_text = re.sub(', +', ',', output_text)
        return output_text
    
    #
    def _write_drug_tree_table(self, augmented_complexity_dict, 
                               drug_complexity_dict):
        output_text = ''
        output_text += \
            self._write_drug_tree_column_labels(augmented_complexity_dict,
                                                drug_complexity_dict)
        keys = augmented_complexity_dict.keys()
        keys = sorted(keys)
        maxlen = len(max(keys, key=len))
        for i in range(len(keys)):
            gap = '.' * (maxlen - len(keys[i]) + 2)
            output_text += ' ' + keys[i] + gap
            output_text += self._write_drug_element_complexity(
                augmented_complexity_dict, keys[i], 'P')
            output_text += self._write_drug_element_complexity(
                augmented_complexity_dict, keys[i], '|L|')
            output_text += self._write_drug_element_complexity(
                augmented_complexity_dict, keys[i], '|Pi|')
            output_text += self._write_drug_element_complexity(
                augmented_complexity_dict, keys[i], 'eta')
            output_text += '| '
            output_text += self._write_drug_element_complexity(
                drug_complexity_dict, keys[i], '|L|')
            output_text += self._write_drug_element_complexity(
                drug_complexity_dict, keys[i], '|Pi|')
            output_text += self._write_drug_element_complexity(
                drug_complexity_dict, keys[i], 'eta')
            output_text += '\n'
        return output_text
    
    #
    def _write_target_tree_csv(self, augmented_complexity_dict, 
                               target_complexity_dict):
        output_text = ''
        output_text_tmp = \
            self._write_drug_tree_column_labels(augmented_complexity_dict,
                                                target_complexity_dict)
        output_text_tmp = re.sub(' +', ',', output_text_tmp)
        output_text_tmp = re.sub(',\|,', ',', output_text_tmp)
        output_text += 'TARGET_TREE,TREE_NAME' + output_text_tmp
        keys = augmented_complexity_dict.keys()
        keys = sorted(keys)
        for i in range(len(keys)):
            output_text += 'TARGET_TREE,' + keys[i] + ','
            output_text += str(augmented_complexity_dict[keys[i]]['P']) + ','
            output_text += str(augmented_complexity_dict[keys[i]]['|L|']) + ','
            output_text += str(augmented_complexity_dict[keys[i]]['|Pi|']) + ','
            output_text += str(augmented_complexity_dict[keys[i]]['eta']) + ','
            output_text += str(target_complexity_dict[keys[i]]['|L|']) + ','
            output_text += str(target_complexity_dict[keys[i]]['|Pi|']) + ','
            output_text += str(target_complexity_dict[keys[i]]['eta']) + ','
            output_text += '\n'
            output_text = re.sub(', +', ',', output_text)
        return output_text
    
    #
    def _write_target_tree_table(self, augmented_complexity_dict, 
                                 target_complexity_dict):
        output_text = ''
        output_text += \
            self._write_drug_tree_column_labels(augmented_complexity_dict,
                                                target_complexity_dict)
        keys = augmented_complexity_dict.keys()
        keys = sorted(keys)
        maxlen = len(max(keys, key=len))
        for i in range(len(keys)):
            gap = '.' * (maxlen - len(keys[i]) + 2)
            output_text += ' ' + keys[i] + gap
            output_text += self._write_drug_element_complexity(
                augmented_complexity_dict, keys[i], 'P')
            output_text += self._write_drug_element_complexity(
                augmented_complexity_dict, keys[i], '|L|')
            output_text += self._write_drug_element_complexity(
                augmented_complexity_dict, keys[i], '|Pi|')
            output_text += self._write_drug_element_complexity(
                augmented_complexity_dict, keys[i], 'eta')
            output_text += '| '
            output_text += self._write_drug_element_complexity(
                target_complexity_dict, keys[i], '|L|')
            output_text += self._write_drug_element_complexity(
                target_complexity_dict, keys[i], '|Pi|')
            output_text += self._write_drug_element_complexity(
                target_complexity_dict, keys[i], 'eta')
            output_text += '\n'
        return output_text
    
    #
    def _write_text_list(self, text_list):
        output_text = ''
        for i in range(len(text_list)-1):
            text_list[i] += ','
        for i in range(len(text_list)):
            output_text += ' ' + text_list[i]
            if i % 5 == 4:
                output_text += '\n'
            elif i == len(text_list) - 1:
                output_text = output_text
                output_text += '\n'
        return output_text
    
    #
    def generate_drug_coverage_report(self, filename):
        output_text = self._generate_drug_coverage_report()
        with open(filename + '.txt', 'w') as f:
            f.write(output_text)
        output_text = self._generate_drug_coverage_csv()
        with open(filename + '.csv', 'w') as f:
            f.write(output_text)
    
    #
    def generate_drug_target_coverage_report(self, filename):
        output_text = self._generate_drug_target_coverage_report()
        with open(filename + '.txt', 'w') as f:
            f.write(output_text)
        output_text = self._generate_drug_target_coverage_csv()
        with open(filename + '.csv', 'w') as f:
            f.write(output_text)
            
    #
    def generate_drug_tree_report(self, filename):
        output_text = self._generate_drug_tree_report()
        with open(filename + '.txt', 'w') as f:
            f.write(output_text)
        output_text = self._generate_drug_tree_csv()
        with open(filename + '.csv', 'w') as f:
            f.write(output_text)
            
    #
    def generate_target_tree_report(self, filename):
        output_text = self._generate_target_tree_report()
        with open(filename + '.txt', 'w') as f:
            f.write(output_text)
        output_text = self._generate_target_tree_csv()
        with open(filename + '.csv', 'w') as f:
            f.write(output_text)
            
    #
    def push_data_dict(self, data_dict):
        self.data_dict = data_dict