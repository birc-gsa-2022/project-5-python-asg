import numpy as np

def local_alignment(seq1, seq2):

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
            if row == 0 and col == 0:
                return ''.join(aligned1[::-1]), ''.join(aligned2[::-1])

    edit_distance, matrix = create_matrix(seq1, seq2)
    seq1_aligned, seq2_aligned = backtrace(seq1, seq2, matrix)
    return (seq1_aligned, seq2_aligned)

# Usage:
# string = 'agctagctagctgc'
# pattern = 'cgatcgatcg'
# print(local_alignment(string,pattern))



