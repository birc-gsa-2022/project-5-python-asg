import numpy as np
from collections import deque

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


def local_alignment(seq1, seq2, d):

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
        return matrix

    def backtrace(seq1, seq2, matrix, d):
        row, col = len(seq1), len(seq2)
        # Allowing mismatches:
        last_row = matrix[len(seq1),:]
        last_row_min_d_pos = len(last_row)-(d)+np.argmin(last_row[-(d):])
        last_col = matrix[:,len(seq2)]
        last_col_min_d_pos = len(last_col)-(d)+np.argmin(last_col[-(d):])
        t0 = matrix[last_col_min_d_pos, len(seq2)]
        t1 = matrix[last_col_min_d_pos, len(seq2)]       
        t2 = matrix[len(seq1), last_row_min_d_pos]
        edit_distance = min(t0,t1,t2)
        if edit_distance <= d:
            return edit_distance
        else: 
            return None

    matrix = create_matrix(seq1, seq2)
    if backtrace(seq1, seq2, matrix, d) != None:
        edit_distance = backtrace(seq1, seq2, matrix, d)
        return (edit_distance, matrix)
    else: 
        return None


def approx_positions(string, pattern, SA, d):
    n_segments = d+1
    segment_size = int(round(len(pattern)/n_segments))
    approx_pos = set()
    approx_trimmed = set()
    
    # Generate all possible start positions given d+1 possible exact matching segments:
    exp_starts = set()
    for k in range(n_segments):
        segment_start = k*segment_size
        segment_end = segment_start + segment_size
        if segment_end > len(pattern): segment_end = len(pattern)
        matches = binary_search(SA, string, pattern[segment_start:segment_end])
        if matches != None:
            for m in matches:
                # print('m_pos:', m, pattern[segment_start:segment_end], 'segment_start:', segment_start, segment_end)
                exp_start = m-segment_start
                if exp_start >= 0-d and exp_start <= len(string)+d:
                    exp_starts.add(exp_start)
    # print('exp_starts:', exp_starts)
    
    # Generate all possible pattern intervals:
    possible_intervals = set()
    for val in exp_starts:
        for flex in range(-d,d+1,1):
            flex_start = val + flex 
            flex_end = val + len(pattern) + (d-abs(flex))
            if flex_start >= 0 and flex_end <= len(string):
                possible_intervals.add((flex_start, flex_end))
            else:
                if flex_start >= 0 and flex_end > len(string):
                    possible_intervals.add((flex_start, len(string)))
                if flex_start < 0 and flex_end <= len(string):
                    possible_intervals.add((0, flex_end))
    possible_intervals = list(possible_intervals)
    # print('possible_intervals:', possible_intervals)
    
    # Get all matrices with best-fit <= d+chopped_pattern_ends:
    for tup in possible_intervals:
        seq = string[tup[0]:tup[1]]
        if 0 <= abs(len(seq) - len(pattern)) <= d*2:
            alignment = local_alignment(pattern, seq, d*2)
            if alignment != None:

                # Report all paths of matrix within d edits:
                d_max = (d*2)
                if alignment[0] <= d_max:
                    # print('Alignments:', alignment[0], alignment[1])
                    # print('Matrix_and_pos: \n', alignment[3], tup[0])
                    
                    seq1, seq2 = list(pattern), list(seq)
                    matrix = alignment[1]  
                    row, col = len(seq1), len(seq2)
                    stack = set()
                    for i in range(d_max):
                        stack.add(('', '', row-i, col, 0))  # alignment1, alignment2, row, col, mismatches.
                        stack.add(('', '', row, col-i, 0))  # alignment1, alignment2, row, col, mismatches.

                    while len(stack) > 0:
                        cur = stack.pop()
                        row, col = cur[2], cur[3]

                        if row == 0 and col == 0:
                            align1, align2 = ''.join(cur[0])[::-1], ''.join(cur[1])[::-1]
                            if align1.count('-') + align2.count('-') <= d_max:
                                approx_pos.add((tup[0], align1, align2))
                        
                        else: 
                            vertical = matrix[row-1, col]
                            diagonal = matrix[row-1, col-1]
                            horizontal = matrix[row, col-1]
                            
                            if matrix[row,col] == diagonal and diagonal < d_max and row>=0 and col>=0:
                                path_tup = ( cur[0]+seq1[row-1], cur[1]+seq2[col-1], row-1, col-1, cur[4] )
                                if path_tup[4] < d_max:
                                    stack.add(path_tup)
                            if matrix[row,col] == diagonal+1 and diagonal+1 < d_max and row>=0 and col>=0:
                                path_tup = ( cur[0]+seq1[row-1], cur[1]+seq2[col-1], row-1, col-1, cur[4]+1 )
                                if path_tup[4] < d_max:
                                    stack.add((path_tup)) 
                            
                            if matrix[row,col] == vertical and vertical < d_max and row>=0:
                                path_tup = ( cur[0]+seq1[row - 1], cur[1]+"-", row-1, col, cur[4]+1 )
                                if path_tup[4] < d_max:
                                    stack.add(path_tup)
                            if matrix[row,col] == vertical+1 and vertical+1 < d_max and row>=0:
                                path_tup = ( cur[0]+seq1[row - 1], cur[1]+"-", row-1, col, cur[4]+1 )
                                if path_tup[4] < d_max:
                                    stack.add(path_tup)
                            if matrix[row,col] == vertical-1 and vertical-1 < d_max and row>=0:
                                path_tup = ( cur[0]+seq1[row - 1], cur[1]+"-", row-1, col, cur[4]+1 )
                                if path_tup[4] < d_max:
                                    stack.add(path_tup)
                            
                            if matrix[row,col] == horizontal and horizontal < d_max and col>=0:
                                path_tup = ( cur[0]+"-", cur[1]+seq2[col-1], row, col-1, cur[4]+1 )
                                if path_tup[4] < d_max:
                                    stack.add(path_tup)
                            if matrix[row,col] == (horizontal+1) and (horizontal+1) < d_max and col>=0:
                                path_tup = ( cur[0]+"-", cur[1]+seq2[col-1], row, col-1, cur[4]+1 )
                                if path_tup[4] < d_max:
                                    stack.add(path_tup)
                            if matrix[row,col] == (horizontal-1) and (horizontal-1) < d_max and col>=0:
                                path_tup = ( cur[0]+"-", cur[1]+seq2[col-1], row, col-1, cur[4]+1 )
                                if path_tup[4] < d_max:
                                    stack.add(path_tup)

    for alignment in approx_pos:
        start_gaps = 0
        j = 0
        while alignment[1][j] == '-':
            start_gaps+=1
            j+=1
        end_gaps = 0
        j = 0
        while alignment[1][::-1][j] == '-':
            end_gaps+=1
            j+=1
        
        al1 = alignment[1][0+start_gaps:len(alignment[1])-end_gaps]
        al2 = alignment[2][0+start_gaps:len(alignment[2])-end_gaps]
        pos = alignment[0]+(start_gaps)
        # print(pos, al1,al2)
        if len(al1) >= (len(pattern) + al1.count('-')):
            mm = 0
            for i in range(len(al1)):
                if al1[i] != al2[i]:
                    mm+=1
            if mm <= d:
                approx_trimmed.add((pos,al1,al2))
                    
    return list(approx_trimmed)



    
# if in ['-', pattern[-1]] (test if new match right)
###########################################################
# Usage:
# string =  'ttgatgaaacgtcgctgctacataggagattcccggcaggcgctatgccttggatgagactaaaggtcacctactccattcctacttccttcagtggagaacgctgcggtccggaagatttgactgagacccgcttaaagttttccgtgcatatttgtagtactaagcgcggctcgatgatgttacacgcttaatccacagttggaggtcatccatgggtgcaccaatgcgtttaagtcagagttaccgatcgttcttaagtgagcttctcggcgaattgtacggaggtgtgctatcactcgttccataagtggcgtactgattatcttcactgacccgcctagacttgtaagcttcgaacagactccgccaatgagagcgtgcaatggtgtacggcattacggagacggtagcgtccaacggaagggccagagtattagattcatttgaaaagaacactgacttttgctaacaaaagctcgggcgtggtaagcggttca'
# pattern = '                                                                                                       tgcgggat'
# pattern = 'tgcgggat'

# string = 'cgcttaccgttcttaaggctattgcgaaacaggctatgattttgaccctgaccgttcagtcgtcaaaatctcgattctagtacgtgggatttctagacgtcaac-ccccgtgaccgattctgaaaggaattctaggggaatagtcgacttagcgacactggttgaccggggagcacacgatcggacttcccgagacccagatcgaagagcgtctcgtccacgtgccggagacaatagcggggatgaccgggtcgtgcccttctggccgtaactcgaattgtcgaggtggctagcgtccgttccagaagtacttaagtcgataatgagctgggacctagattttcaagcgaacagatgtacacctttcgtggtcgtactaaggcttcgtgtaagtagttggggactgttcacactaatagcatcggaattcttcatgaatcaatcgtatcgattgactacctgcacattgattcactcgtctgaggacgtgcgttaccagataagccgagagtacccgatctatcagaatggatacccctcacagggatctattagtagttacacgggtggtgccaccccagagcgccaacttcttacttgacacgttgggtacgacagaaagttcaacgcgattttccggattagacagttccgtgttaccgctcagacacatcccgcatcccaagttacaaacacgtccaaattgaacctagcgatctccgcggaacatacagttttttgagaaataagtacgtttctcacggccgttgctaactcccgattgaccacctgcgacgttaaaaatcttacagtgcgacactggatcaccaaatggccggtgacatacgcccatacgaataagatcacgggcacttttgccactaggattcaaattccacgtagttcttacaaggaacaacaacaaagataaacttgtacagttgatgcgagcaactgattctacaaaaatccgcatgccagagatgaatttgtaccctatt'
# pattern = '                                                                                                     a-cccccgtga'
# string = 'cgcttaccgttcttaaggctattgcgaaacaggctatgattttgaccctgaccgttcagtcgtcaaaatctcgattctagtacgtgggatttctagacgtcaacccccgtgaccgattctgaaaggaattctaggggaatagtcgacttagcgacactggttgaccggggagcacacgatcggacttcccgagacccagatcgaagagcgtctcgtccacgtgccggagacaatagcggggatgaccgggtcgtgcccttctggccgtaactcgaattgtcgaggtggctagcgtccgttccagaagtacttaagtcgataatgagctgggacctagattttcaagcgaacagatgtacacctttcgtggtcgtactaaggcttcgtgtaagtagttggggactgttcacactaatagcatcggaattcttcatgaatcaatcgtatcgattgactacctgcacattgattcactcgtctgaggacgtgcgttaccagataagccgagagtacccgatctatcagaatggatacccctcacagggatctattagtagttacacgggtggtgccaccccagagcgccaacttcttacttgacacgttgggtacgacagaaagttcaacgcgattttccggattagacagttccgtgttaccgctcagacacatcccgcatcccaagttacaaacacgtccaaattgaacctagcgatctccgcggaacatacagttttttgagaaataagtacgtttctcacggccgttgctaactcccgattgaccacctgcgacgttaaaaatcttacagtgcgacactggatcaccaaatggccggtgacatacgcccatacgaataagatcacgggcacttttgccactaggattcaaattccacgtagttcttacaaggaacaacaacaaagataaacttgtacagttgatgcgagcaactgattctacaaaaatccgcatgccagagatgaatttgtaccctatt'
# pattern = 'acccccgtga'
# SA = SuffixArray(string)
# print(approx_positions(string, pattern, SA, 2)) 






