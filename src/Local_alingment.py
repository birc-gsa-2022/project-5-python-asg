#####################################################
# Score matrix (global value):
match_scores = {'A': {'A': 0, 'C': -4, 'G': -2, 'T': -4},
                'C': {'A': -4, 'C': 0, 'G': -4, 'T': -2},
                'G': {'A': -2, 'C': -4, 'G': 0, 'T': -4},
                'T': {'A': -4, 'C': -2, 'G': -4, 'T': -4}}

#####################################################
# Functions:

def empty_matrix(len1,len2):
    matrix = []
    for i in range(len1):
        matrix.append([])
        for j in range(len2):
            matrix[i].append(None)
    return matrix

def prepare_matrix(len1,len2,gap_score):
    matrix = empty_matrix(len1,len2)
    for i in range(len(matrix[0])):
        cell_value = i*gap_score
        matrix[0][i] = cell_value
    for j in range(len(matrix)):
        cell_value = j*gap_score
        matrix[j][0]= cell_value
    return(matrix)

def fill_matrix(seq1,seq2,match_scores,gap_score):
    matrix = prepare_matrix(len(seq1)+1,len(seq2)+1,gap_score)
    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            if matrix[i][j] == None:
                matrix[i][j] = 0
            value_dict = {'left': matrix[i][j-1] + gap_score,'above':matrix[i-1][j] + gap_score,'diagonal':matrix[i-1][j-1] + match_scores[seq1[i-1]][seq2[j-1]]}
            matrix[i][j] = max([j for j in value_dict.values()])
    return matrix

def get_traceback_arrow(matrix, row, col, match_score, gap_score):
    score_diagonal = matrix[row-1][col-1]
    score_left = matrix[row][col-1]
    score_up = matrix[row-1][col]
    score_current = matrix[row][col]
    if score_current == score_diagonal + match_score:
        return 'diagonal'
    elif score_current == score_left + gap_score:
        return 'left'
    elif score_current == score_up + gap_score:
        return 'up'

def trace_back(seq1, seq2, matrix, score_matrix, gap_score):
    aligned1 = ''
    aligned2 = ''
    row = len(seq1)
    col = len(seq2)
    while row > 0 and col > 0:
        base1 = seq1[row-1]
        base2 = seq2[col-1]
        match_score = score_matrix[base1][base2]
        traceback_arrow = get_traceback_arrow(matrix, row, col, match_score, gap_score)
        if traceback_arrow == 'diagonal':
            aligned1 = base1 + aligned1
            aligned2 = base2 + aligned2
            row -= 1
            col -= 1
        elif traceback_arrow == 'up':
            aligned1 = base1 + aligned1
            aligned2 = '-' + aligned2
            row -= 1
        elif traceback_arrow == 'left':
            aligned1 = '-' + aligned1
            aligned2 = base2 + aligned2
            col -= 1
    while row > 0:
        base1 = seq1[row-1] 
        aligned1 = base1 + aligned1
        aligned2 = '-' + aligned2
        row -= 1
    while col > 0:
        base2 = seq2[col-1]
        aligned1 = '-' + aligned1
        aligned2 = base2 + aligned2 
        col -= 1
    return [aligned1, aligned2]

def align(seq1,seq2,match_scores,gap_score):
    matrix = fill_matrix(seq1,seq2,match_scores,gap_score)
    return trace_back(seq1, seq2, matrix, match_scores, gap_score)

def get_edits(p: str, q: str):
    edits = ''  # edits to go from p to q
    for i in range(len(p)):
        if p[i] != '-' and q[i] != '-':
            edits += 'M'
        if p[i] == '-' and q[i] != '-':
            edits += 'I'
        if p[i] != '-' and q[i] == '-':
            edits += 'D'
    return edits
####################################################
# Usage:
string1 = 'GATGAT'
string2 = 'GGAAGG'
alignment = align(string1, string2, match_scores, -6)
print(get_edits(alignment[0], alignment[1]))

