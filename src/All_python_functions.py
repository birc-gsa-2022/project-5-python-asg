################################################################
# Libraries:

import numpy as np

################################################################
# Functions:

def SuffixArray(string):
    if string == '' or string == None:
        return None
    string += '$'
    index = {v: i for i, v in enumerate(sorted(set(string)))}
    string = [index[v] for v in string]
    rank_list = list(range(len(string)))
    SA = [None]*len(string)
    tuple_list = SA[:]
    M,j = 0,1
    while M < len(string)-1:
        for i in range(len(string)):
            i_j = i+j
            if i_j < len(string):
                key = (string[i], string[i_j])
            else: 
                key = (string[i], string[-1])
            tuple_list[i] = key
        j*=2
        keys = sorted(set(tuple_list))
        ranks = dict(zip(keys, rank_list))
        string = [ranks[tuple_list[i]] for i in range(len(string))]
        M = max(string)
    for i in rank_list:
        SA[string[i]] = i
    return SA


def binary_search(SA, string, pattern):
    if pattern == '' or pattern == None:
        return None
    string+='$'
    SA_positions = []
    # Binary search:
    hi, lo = len(SA), 0, 
    match = None
    B = 0
    while lo < hi and B == 0:
        mid = (lo + hi) // 2
        count = 0
        for i in range(count, len(pattern)):
            if pattern[i] == string[SA[mid]+i]:
                count+=1
            if count == len(pattern):
                match = mid    
                SA_positions.append(SA[match])
                B=1
            elif pattern[i] > string[SA[mid]+i]:
                lo = mid + 1
                break
            elif pattern[i] < string[SA[mid]+i]:
                hi = mid
                break
    # Scan up/down from match pos:
    if match != None:
        k=1
        while match-k >= 0 and string[SA[match-k]:SA[match-k]+len(pattern)] == pattern:
            SA_positions.append(SA[match-k])
            k+=1
        k=1
        while match+k < len(string) and string[SA[match+k]:SA[match+k]+len(pattern)] == pattern:
            SA_positions.append(SA[match+k])
            k+=1  
    return SA_positions


def local_alignment(seq1, seq2, max_edit_dist):

    def create_matrix(seq1, seq2):
        rows = len(seq1) + 1
        cols = len(seq2) + 1
        matrix = np.zeros((rows, cols), dtype=int)
        for i in range(1, rows):
            for j in range(1, cols):
                matrix[i, 0] = i
                matrix[0, j] = j
        for col in range(1, cols):
            for row in range(1, rows):
                cost = 0 if seq1[row - 1] == seq2[col - 1] else 1
                matrix[row, col] = min(matrix[row - 1, col] + 1,
                                       matrix[row, col - 1] + 1,
                                       matrix[row - 1, col - 1] + cost)
        edit_distance = matrix[row, col]
        return edit_distance, matrix

    def backtrace(seq1, seq2, matrix):
        seq1, seq2 = list(seq1), list(seq2)
        aligned1, aligned2 = [], []
        row, col = len(seq1), len(seq2)
        while True:
            cost = 0 if seq1[row - 1] == seq2[col - 1] else 1
            cur = matrix[row, col]
            vertical = matrix[row-1, col]
            diagonal = matrix[row - 1, col - 1]
            horizontal = matrix[row, col - 1]
            if cur == diagonal + cost:
                aligned1 += [seq1[row - 1]]
                aligned2 += [seq2[col - 1]]
                row, col = row - 1, col - 1
            else:
                if cur == vertical + 1:
                    aligned1 += [seq1[row - 1]]
                    aligned2 += ["-"]
                    row, col = row - 1, col
                elif cur == horizontal + 1:
                    aligned1 += ["-"]
                    aligned2 += [seq2[col - 1]]
                    row, col = row, col - 1
            if row == 0 or col == 0:
                return ''.join(aligned1[::-1]), ''.join(aligned2[::-1])

    edit_distance, matrix = create_matrix(seq1, seq2)
    if edit_distance > max_edit_dist:
        return None
    seq1_aligned, seq2_aligned = backtrace(seq1, seq2, matrix)
    return (seq1_aligned, seq2_aligned)


def approx_search(string, pattern, SA, max_edit_dist):
    n_segments = max_edit_dist+1
    segment_size = int(round(len(pattern)/n_segments))
    approx_matches = []
    exp_start = set()
    for k in range(n_segments):
        segment_start = k*segment_size
        segment_end = segment_start + segment_size
        if segment_end > len(pattern): segment_end = len(pattern)
        matches = binary_search(SA, string, pattern[segment_start:segment_end])
        if matches != None:
            for m in matches:
                # print('m_pos:', m, pattern[segment_start:segment_end], 'segment_start:', segment_start)
                start = m-segment_start
                if start >= 0-max_edit_dist and start <= len(string)+max_edit_dist:
                    exp_start.add((start,m))
    # print('exp_start: ', exp_start)
    if exp_start != None:
        for s in exp_start:
            try: 
                if s[0] >= 0: 
                    sub_string = string[s[0]:s[0]+len(pattern)] 
                if s[0] < 0: 
                    sub_string = string[0:s[0]+len(pattern)]
                if s[0]+len(pattern) > len(string):
                    sub_string = string[s[0]:len(string)]
            except: continue
            alignment = local_alignment(pattern, sub_string, max_edit_dist)
            if alignment != None:
                approx_matches.append((s[0], alignment[0], alignment[1]))
    return approx_matches

################################################################
# Usage:
# string =  'ATCAGCATGCTAGC'
# pattern = 'ATCGAT'
# SA = SuffixArray(string)
# print(approx_search(string, pattern, SA, 3))
################################################################



