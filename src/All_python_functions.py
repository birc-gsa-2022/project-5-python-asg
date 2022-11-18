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
    segment_size = len(pattern) // n_segments + (len(pattern) % n_segments > 0)
    approx_matches = []
    for j in range(n_segments):
        segment_start = j*segment_size
        segment_end = segment_start + segment_size
        if segment_end > len(pattern): segment_end = len(pattern)
        matches = binary_search(SA, string, pattern[segment_start:segment_end])
        if matches != None:
            for match in matches:
                if match >= segment_start and len(string) >= (match-segment_start+len(pattern)):
                    mm = 0 
                    for i in range(0,segment_start):
                        if pattern[i] != string[match-segment_start+i]:
                            mm+=1
                            if mm > max_edit_dist: break
                    for i in range(segment_end,len(pattern)):
                        if pattern[i] != string[match-segment_start+i]:
                            mm+=1
                            if mm > max_edit_dist:  break
                    if mm <= max_edit_dist:
                        approx_pos = match-segment_start
                        approx_matches.append(approx_pos) if approx_pos not in approx_matches else approx_matches
    return approx_matches

################################################################