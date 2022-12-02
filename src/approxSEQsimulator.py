import random

def simulate_string(m):
    """Simulate a DNA sequence of length m.
    
    Example:
    simulate_string(20)
    
    """
    DNA = 'ACGT'
    nucleotides = [random.choice(DNA) for _ in range(m)]
    return ''.join(nucleotides)


def get_approx_read(DNA_seq:str, n:int, d:int):
    """Simulates an exact read of length n from a DNA sequence with d edits.
    
    Example:
    get_exact_read('GCTGAATCTATTTATGGTGG',5)

    """
    assert n < len(DNA_seq), 'm must be smaller or equal to length of DNA_seq'
    DNA = 'ACGT'
    idx = random.randint(0,len(DNA_seq)-n)
    read = list(DNA_seq[idx:idx+n])

    for i in range(d):
        r = random.randint(1, len(read)-1)
        t = random.randint(0,2)
        if t == 0: 
            letter = read[r]
            while letter == read[r]:
                letter = DNA[random.randint(0,3)]
            read[i] = letter
        if t == 1:
            read.pop(r)
        if t == 2:
            letter = DNA[random.randint(0,3)]
            read.insert(r, letter)
    return ''.join(read)

