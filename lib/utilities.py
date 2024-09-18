# -*- coding: utf-8 -*-
"""
Created on Thu May  2 11:47:57 2024

@author: haglers
"""

#
def _unique(data_list_list):
    data_list_list = \
        [list(x) for x in set(tuple(x) for x in data_list_list)]
    return data_list_list