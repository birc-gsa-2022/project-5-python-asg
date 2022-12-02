############################################################
# Librarys:
import re
import random
import numpy as np

############################################################
# Functions:

def edits_to_cigar(edits: str):
    cdef str cigar
    cdef str a
    def split_blocks(x: str):
        return [m[0] for m in re.findall(r"((.)\2*)", x)]
    element_list = split_blocks(edits)
    cigar=''
    for a in element_list:
        short = '{}{}'.format(len(a),a[0])
        cigar+= short
    return cigar

def get_edits(p: str, q: str):
    cdef int i
    cdef str edits
    assert len(p) == len(q)
    edits = ''
    for i in range(len(p)):
        if p[i] != '-' and q[i] != '-':
            edits += 'M'
        if p[i] == '-' and q[i] != '-':
            edits += 'I'
        if p[i] != '-' and q[i] == '-':
            edits += 'D'
    if len(p) == 0 and len(q) == 0:
        p_out = ''
        q_out = ''
        edits = ''
    else:
        p_out = p.replace('-','')
        q_out = q.replace('-','')
        edits = edits
    return p_out, q_out , edits

def SuffixArray(string):
    cdef int i
    cdef int j
    cdef int i_j
    cdef int M
    cdef str v
    cdef list SA
    cdef list rank_list
    cdef list tuple_list
    cdef dict index
    cdef tuple key
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
    cdef int i
    cdef int mid
    cdef int count
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


def local_matrix(seq1, seq2, d):

    def create_matrix(seq1, seq2):
        cdef int i
        cdef int j
        cdef int row
        cdef int col
        cdef int cost
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

    def check_dist(seq1, seq2, matrix, d):
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
    edit_distance = check_dist(seq1, seq2, matrix, d)
    if edit_distance != None:
        return (edit_distance, matrix)
    else: 
        return None


def approx_positions(string, pattern, SA, d):
    cdef int s
    cdef int m
    cdef int val
    cdef int flex
    cdef int n_segments
    cdef int segment_size
    cdef int row
    cdef int col
    cdef tuple tup
    n_segments = d+1
    segment_size = int(round(len(pattern)/n_segments))
    approx_pos = set()
    approx_trimmed = set()
    
    # Generate all possible start positions given d+1 possible exact matching segments:
    exp_starts = set()
    for s in range(n_segments):
        segment_start = s*segment_size
        segment_end = segment_start + segment_size
        if segment_end > len(pattern): segment_end = len(pattern)
        matches = binary_search(SA, string, pattern[segment_start:segment_end])
        if matches != None:
            for m in matches:
                exp_start = m-segment_start
                if exp_start >= 0-d and exp_start <= len(string)+d:
                    exp_starts.add(exp_start)
    
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
    
    # Get all matrices with best-fit <= d*2:
    for tup in possible_intervals:
        seq, d_max = string[tup[0]:tup[1]], d*2
        if 0 <= abs(len(seq) - len(pattern)) <= d_max:
            matrix_inf = local_matrix(pattern, seq, d_max)
            if matrix_inf != None:
                if matrix_inf[0] <= d_max:
                    # print('Matrix_and_pos: \n', matrix_inf[1], tup[0])
                    
                    seq1, seq2 = list(pattern), list(seq)
                    matrix = matrix_inf[1]  
                    stack = [('', '', 0, 0, 0)]
                    end_row = list(range(len(seq2)+1))[-d_max:]
                    end_col = list(range(len(seq1)+1))[-d_max:]

                    while len(stack) > 0:
                        cur = stack.pop()
                        row, col = cur[2], cur[3]

                        # print(cur)
                        if row < len(seq1) and col < len(seq2):
            
                            diagonal = matrix[row+1, col+1]
                            vertical = matrix[row+1, col]
                            horizontal = matrix[row, col+1]
                            
                            if matrix[row,col] == diagonal and diagonal <= d_max and row+1<=len(seq1) and col+1<=len(seq2):
                                path_tup = ( cur[0]+seq1[row], cur[1]+seq2[col], row+1, col+1, cur[4] )
                                if path_tup[4] <= d_max:
                                    stack.append(path_tup)
                            if matrix[row,col] == diagonal-1 and diagonal+1 <= d_max and row+1<=len(seq1) and col+1<=len(seq2):
                                path_tup = ( cur[0]+seq1[row], cur[1]+seq2[col], row+1, col+1, cur[4]+1 )
                                if path_tup[4] <= d_max:
                                    stack.append(path_tup) 
                            
                            if matrix[row,col] == vertical and vertical <= d_max and row+1<=len(seq1):
                                path_tup = ( cur[0]+seq1[row], cur[1]+"-", row+1, col, cur[4]+1 )
                                if path_tup[4] <= d_max:
                                    stack.append(path_tup)
                            if matrix[row,col] == vertical+1 and vertical+1 <= d_max and row+1<=len(seq1):
                                path_tup = ( cur[0]+seq1[row], cur[1]+"-", row+1, col, cur[4]+1 )
                                if path_tup[4] <= d_max:
                                    stack.append(path_tup)
                            if matrix[row,col] == vertical-1 and vertical-1 <= d_max and row+1<=len(seq1):
                                path_tup = ( cur[0]+seq1[row], cur[1]+"-", row+1, col, cur[4]+1 )
                                if path_tup[4] <= d_max:
                                    stack.append(path_tup)
                            
                            if matrix[row,col] == horizontal and horizontal <= d_max and col+1<=len(seq2):
                                path_tup = ( cur[0]+"-", cur[1]+seq2[col], row, col+1, cur[4]+1 )
                                if path_tup[4] <= d_max:
                                    stack.append(path_tup)
                            if matrix[row,col] == (horizontal+1) and (horizontal+1) <= d_max and col+1<=len(seq2):
                                path_tup = ( cur[0]+"-", cur[1]+seq2[col], row, col+1, cur[4]+1 )
                                if path_tup[4] <= d_max:
                                    stack.append(path_tup)
                            if matrix[row,col] == (horizontal-1) and (horizontal-1) <= d_max and col+1<=len(seq2):
                                path_tup = ( cur[0]+"-", cur[1]+seq2[col], row, col+1, cur[4]+1 )
                                if path_tup[4] <= d_max:
                                    stack.append(path_tup)
                        
                        else:
                            if row == len(seq1) and col in end_col or col == len(seq2) and row in end_row:
                                # align1, align2 = ''.join(cur[0]), ''.join(cur[1])
                                if cur[0].count('-') <= d and cur[1].count('-') <=d:
                                    if cur[4] <= d_max:
                                        approx_pos.add((tup[0], cur[0], cur[1]))

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
        ends_gaps = (start_gaps+end_gaps)
        al1_indels = alignment[1].count('-') - ends_gaps
        pos = alignment[0]+start_gaps
        if len(alignment[1])-ends_gaps >= len(pattern) + al1_indels:
            mm = 0
            for i in range(len(alignment[1])-ends_gaps):
                if alignment[1][i+start_gaps] != alignment[2][i+start_gaps]:
                    mm+=1
                if mm > d: break
            if mm <= d:
                al1 = ''.join(alignment[1][0+start_gaps:len(alignment[1])-end_gaps])
                al2 = ''.join(alignment[2][0+start_gaps:len(alignment[2])-end_gaps])
                edits = get_edits(al2,al1)
                cigar = edits_to_cigar(edits[2])
                approx_trimmed.add((pos,cigar))
    return list(approx_trimmed)

############################################################