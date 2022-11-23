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


def approx_search(string, pattern, SA, max_edit_dist):
    n_segments = max_edit_dist+1
    segment_size = int(round(len(pattern)/n_segments))
    approx_matches = set()
    exp_start = set()
    
    for k in range(n_segments):
        segment_start = k*segment_size
        segment_end = segment_start + segment_size
        if segment_end > len(pattern): segment_end = len(pattern)
        matches = binary_search(SA, string, pattern[segment_start:segment_end])
        if matches != None:
            for m in matches:
                start = m-segment_start
                if start >= 0 and start <= len(string):
                    exp_start.add((start,m))
    
    if exp_start != None:
        for s in exp_start:
            try: sub_string = string[s[0]:s[0]+len(pattern)] 
            except: continue
            stack = [(0,0,0)]  # [i, j, edits, edit_score].
            while len(stack) > 0:
                sub = stack.pop()
                j = sub[1]
                mm = sub[2]
                for i in range(sub[0], len(sub_string)):
                    if i == len(pattern)-1:
                        approx_matches.add(s[0])
                    if mm+1 <= max_edit_dist and pattern[i] != sub_string[j]:
                        stack.append( (i, j, sub[2]+1) )  # m=mismatch.
                        stack.append( (i+1, j, sub[2]+1) )  # D=deletion.
                        stack.append( (i, j+1, sub[2]+1) )  # I=insertion.
                        continue
                    j+=1                    
    
    return approx_matches

################################################################
# string =  'mississippi'
# pattern = 'mssmsi'
# SA = SuffixArray(string)
# print(approx_search(string, pattern, SA, 3))

