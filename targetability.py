# -*- coding: utf-8 -*-
"""
<one line to give the program's name and a brief idea of what it does.>
Copyright (C) <year>  <name of author>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

#
import os
import sys

#
path = os.path.dirname(__file__)
sys.path.insert(0, path + '/lib')

#
import argparse
from datetime import datetime
from output import Output_object
from targetability import Targetability_object

#
class Flat_files_object(object):
    
    #
    def __init__(self):
        self.pkl_path = os.path.join(path, 'pkl')
        self.reactome_flat_files_path = os.path.join(path, 'reactome_flat_files')
        self.targetome_flat_files_path = os.path.join(path, 'targetome_flat_files')
        self.all_levels_pathways_filename = 'UniProt2Reactome_All_Levels.txt'
        self.pathways_relation_filename = 'ReactomePathwaysRelation.txt'
        self.pkl_filename = 'PathwayLattice.pkl'
        self.reactome_pathways_filename = 'ReactomePathways.txt'
        self.targetome_evidence_filename = 'Targetome_Evidence_TIPS_101017.txt'
        
    #
    def pull_all_levels_pathways_filename(self):
        return self.all_levels_pathways_filename
    
    #
    def pull_pathways_relation_filename(self):
        return self.pathways_relation_filename
    
    #
    def pull_pkl_filename(self):
        return self.pkl_filename
    
    #
    def pull_pkl_path(self):
        return self.pkl_path
    
    #
    def pull_reactome_pathways_filename(self):
        return self.reactome_pathways_filename
    
    #
    def pull_reactome_flat_files_path(self):
        return self.reactome_flat_files_path
    
    #
    def pull_targetome_evidence_filename(self):
        return self.targetome_evidence_filename
    
    #
    def pull_targetome_flat_files_path(self):
        return self.targetome_flat_files_path

#
def _get_filename(args, filename=None):
    if args.filename:
        filename = args.filename
    else:
        print('--filename not specified')
        now = datetime.now()
        now_str = now.strftime("%Y%m%d_%H%M%S")
        if filename is not None:
            filename = filename + '_' + now_str
        else:
            filename = 'targetability_' + now_str
    return filename

#
def _get_max_assay_value(args):
    if args.max_assay_value:
        max_assay_value = args.max_assay_value
    else:
        print('--max_assay_value not specified')
        max_assay_value = 100
    return max_assay_value

#
def main():
    
    # license information
    text = '''<program>  Copyright (C) <year>  <name of author>
        This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
        This is free software, and you are welcome to redistribute it
        under certain conditions; type `show c' for details.'''
    print(text)
    
    # initialize parser
    parser = argparse.ArgumentParser(description='Optional app description')
    parser.add_argument('--drug_coverage', type=str, nargs='+',
                        help='A string giving a list of drugs in Targetome separated by spaces')
    parser.add_argument('--drug_tree', type=str, nargs='+',
                        help='A string giving a list of drugs in Targetome separated by spaces')
    parser.add_argument('--filename', type = int, nargs=1,
                        help='A string giving a valid filename')
    parser.add_argument('--max_assay_value', type = int, nargs=1,
                        help='An integer indiating the maximum assay value')
    parser.add_argument('--target_coverage', type=str, nargs='+',
                        help='A string giving a list of targets (genes) in Uniprot separated by spaces')
    parser.add_argument('--target_tree', type=str, nargs='+',
                        help='A string giving a list of targets (genes) in Uniprot separated by spaces')
    parser.add_argument('--paper_0', action='store_true',
                        help='Generate paper_0 data set')
    parser.add_argument('--update', action='store_true',
                        help='Update data structures')
    args = parser.parse_args()
    
    # initialize flat_files_obbject
    flat_files_object = Flat_files_object()
        
    # initialize targetability_object
    targetability_object = Targetability_object()
    targetability_object.push_all_levels_pathways_filename(flat_files_object.pull_all_levels_pathways_filename())
    targetability_object.push_pathways_relation_filename(flat_files_object.pull_pathways_relation_filename())
    targetability_object.push_pkl_path(flat_files_object.pull_pkl_path())
    targetability_object.push_reactome_flat_files_path(flat_files_object.pull_reactome_flat_files_path())
    targetability_object.push_reactome_pathways_filename(flat_files_object.pull_reactome_pathways_filename())
    targetability_object.push_targetome_evidence_filename(flat_files_object.pull_targetome_evidence_filename())
    targetability_object.push_targetome_flat_files_path(flat_files_object.pull_targetome_flat_files_path())
    targetability_object.initialize()
       
    # apply line arguments
    if args.drug_coverage and not args.target_coverage:
        drug_list = args.drug_coverage
        print('--target_coverage not specified')
        max_assay_value = _get_max_assay_value(args)
        filename = _get_filename(args, filename='drug_coverage')
        targetability_object.load_pkl_file(flat_files_object.pull_pkl_filename())
        output_object = Output_object()
        data_dict = \
            targetability_object.get_drug_data_dict(drug_list, max_assay_value)
        output_object.push_data_dict(data_dict)
        output_object.generate_drug_coverage_report(filename)
    elif args.drug_coverage and args.target_coverage:
        drug_list = args.drug_coverage
        target_list = args.target_coverage
        max_assay_value = _get_max_assay_value(args)
        filename = _get_filename(args, filename='drug_target_coverage')
        targetability_object.load_pkl_file(flat_files_object.pull_pkl_filename())
        output_object = Output_object()
        data_dict = \
            targetability_object.get_drug_target_data_dict(drug_list,
                                                           target_list,
                                                           max_assay_value)
        output_object.push_data_dict(data_dict)
        output_object.generate_drug_target_coverage_report(filename)
    elif args.drug_tree:
        drug_list = args.drug_tree
        targetability_object.load_pkl_file(flat_files_object.pull_pkl_filename())
        max_assay_value = _get_max_assay_value(args)
        filename = _get_filename(args, filename='drug_tree')
        output_object = Output_object()
        data_dict = \
            targetability_object.get_drug_data_dict(drug_list, max_assay_value)
        output_object.push_data_dict(data_dict)
        output_object.generate_drug_tree_report(filename)
    elif args.paper_0:
        from papers.paper_0.paper_0_output import Paper_0_output_object
        paper_0_output_object = Paper_0_output_object()
        paper_0_output_path = os.path.join(path, 'papers/paper_0')
        targetability_object.load_pkl_file(flat_files_object.pull_pkl_filename())
        data_dict = targetability_object.paper_0()
        paper_0_output_object.push_data_dict(data_dict)
        filename = 'python_data.dat'
        paper_0_output_object.create_python_data_file(
            os.path.join(paper_0_output_path, filename))
    elif args.target_tree:
        target_list = args.target_tree
        max_assay_value = _get_max_assay_value(args)
        filename = _get_filename(args, 'target_tree')
        targetability_object.load_pkl_file(flat_files_object.pull_pkl_filename())
        output_object = Output_object()
        data_dict = \
            targetability_object.get_target_data_dict(target_list,
                                                      max_assay_value)
        output_object.push_data_dict(data_dict)
        output_object.generate_target_tree_report(filename)
    elif args.update:
        targetability_object.update(flat_files_object.pull_pkl_filename())
    else:
        print('invalid arguments')
        
#
if __name__ == "__main__":
    main()