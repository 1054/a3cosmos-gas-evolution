#!/usr/bin/env python
# 

from __future__ import print_function

import os, sys, time, re

if sys.version_info.major >= 3:
    long = int
else:
    pass


# 
# prepare function for cross-matching
# 
def flatten(item, keepcls=(), keepobj=()):
    """ Flatten a list, 
        see -- https://stackoverflow.com/questions/2158395/flatten-an-irregular-list-of-lists
    """
    if not hasattr(item, '__iter__') or isinstance(item, keepcls) or item in keepobj:
        yield item
    else:
        for i in item:
            for j in flatten(i, keepcls, keepobj + (item,)):
                yield j


def search_for_matches_in_a_sorted_array(input_array, match_value, start_position = 0, search_direction = +1, output_allmatches = False):
    """ Example: xmatch2 = search_for_matches_in_a_sorted_array(array2, array1[i1], i2, -1)
    """
    i = start_position
    xmatches = []
    while (i >= 0 and i <= len(input_array)-1):
        val = input_array[i]
        if val == match_value:
            xmatches.append(i)
            if not output_allmatches:
                break
        else:
            if search_direction > 0:
                if val > match_value:
                    break
            elif search_direction < 0:
                if val < match_value:
                    break
        i = i + search_direction
    return xmatches


def cross_match_sorted_arrays(input_array_list, output_allmatches = False, output_nonmatches = False):
    """ Return two index array with common items in the two arrays
        We will search for only one match per item, unless the option 'output_allmatches' is set to True. 
        We have the option 'output_nonmatches' to also output all non-matches.
        TODO: what if input_array has duplicates?
    """ 
    xmatches = [] # cross-matched index array for each input_array item
    nonmatches = [] # non-matched index array for each input_array item
    i = []
    if not (type(input_array_list) is list):
        print('Error! The input_array_list should be a list!')
        sys.exit()
    if len(input_array_list) == 0:
        print('Error! The input_array_list should be a non-empty list!')
        sys.exit()
    for input_array in input_array_list:
        if len(input_array) == 0:
            print('Error! The input_array should be non-empty!')
            sys.exit()
        xmatches.append([])
        nonmatches.append([])
        i.append(0) # i[0] = 0
    #i = [0]*len(input_array_list)
    while i[0] < len(input_array_list[0]):
        # take the first input array as the reference array
        input_array_1 = input_array_list[0]
        i1 = i[0]
        val1 = long(input_array_1[i1])
        countxmatches = 1 # account for input_array_1 itself
        if output_allmatches:
            tmpxmatches = [[i1]] # a list of lists, each item is a list of matched item indexes for each input array.
        else:
            tmpxmatches = [i1] # when not outputing all matches, this is a list of index
        # try to find matches in each other input array
        for j in range(1,len(input_array_list)):
            input_array_2 = input_array_list[j]
            i2 = i[j]
            tmpxmatches2 = []
            tmpxmatches2a = search_for_matches_in_a_sorted_array(input_array_2, val1, i2, -1, output_allmatches = output_allmatches)
            if len(tmpxmatches2a) > 0:
                tmpxmatches2.extend(tmpxmatches2a)
            if (len(tmpxmatches2a) == 0 or output_allmatches):
                # if we already got matches and output_allmatches is True, then we stop finding new matches toward the other direction
                tmpxmatches2b = search_for_matches_in_a_sorted_array(input_array_2, val1, i2+1, +1, output_allmatches = output_allmatches)
                if len(tmpxmatches2b) > 0:
                    tmpxmatches2.extend(tmpxmatches2b)
            if len(tmpxmatches2) > 0:
                countxmatches = countxmatches + 1
                if output_allmatches:
                    tmpxmatches.append(tmpxmatches2)
                else:
                    tmpxmatches.append(tmpxmatches2[0]) # take only one match for each item, ignoring other matches
        # if found matches in all input arrays
        if countxmatches == len(input_array_list):
            xmatches[0].append(tmpxmatches[0]) # which is just i1 or [i1]
            for j in range(1,len(input_array_list)):
                xmatches[j].append(tmpxmatches[j])
                if output_allmatches:
                    i[j] = max(tmpxmatches[j]) + 1 # here we increase one because all input arrays are sorted
                else:
                    i[j] = tmpxmatches[j] + 1 # here we increase one because all input arrays are sorted
        # 
        i[0] = i[0] + 1
        #print('i[0] =', i[0])
    # 
    if output_nonmatches:
        input_array_1 = input_array_list[0]
        xmatches1 = list(flatten(xmatches[0]))
        allindex1 = range(len(input_array_1))
        nonmatch1 = list(set(allindex1)-set(xmatches1))
        nonmatches[0] = nonmatch1
        for j in range(1,len(input_array_list)):
            input_array_2 = input_array_list[j]
            xmatches2 = list(flatten(xmatches[j]))
            allindex2 = range(len(input_array_2))
            nonmatch2 = list(set(allindex2)-set(xmatches2))
            nonmatches[j] = nonmatch2
        return xmatches, nonmatches
    # 
    return xmatches


