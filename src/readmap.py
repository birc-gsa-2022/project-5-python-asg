#########################################
# Libraries:

import argparse
import sys
import os

#########################################
# Functions

from All_python_functions import SuffixArray
from All_python_functions import approx_search
# from Readmapper_Cython_Convert import SuffixArray
# from Readmapper_Cython_Convert import approx_search
from cigar import edits_to_cigar
from align import get_edits

def read_fasta(inFile):
    lines = inFile.readlines()
    record_list = []
    header = ''
    sequence = []
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if header != "":
                record_list.append([header.strip(), ''.join(sequence).strip()])
                sequence = []
            header = line[1:]
        else:
            sequence.append(line)
    record_list.append([header.strip(), ''.join(sequence).strip()])
    return record_list

def read_fastq(inFile):
    lines = inFile.readlines()
    record_list = []
    header = ''
    sequence = []
    for line in lines:
        line = line.strip()
        if line.startswith('@'):
            if header != "":
                record_list.append([header.strip(), ''.join(sequence).strip()])
                sequence = []
            header = line[1:]
        else:
            sequence.append(line)
    record_list.append([header.strip(), ''.join(sequence).strip()])
    return record_list

def write_SA(genome_name, fa_rec, SA):
    os.makedirs('../{}/{}/'.format(genome_name, fa_rec))
    with open('../{}/{}/SA.txt'.format(genome_name, fa_rec), 'w') as f:
        print(SA, file=f)
        
def open_SA(path_to_preprocessed_dir, fa_rec):
    SA = open('{}/Preprocessed_{}/SA.txt'.format(path_to_preprocessed_dir, fa_rec), 'r').read()
    SA = eval(SA)
    return SA


# Usage:
# string  = 'aaaaaaaaaaaaaaa'
# pattern = 'ba'
# SA = SuffixArray(string)
# print(approx_search(string, pattern, SA, 2))

#########################################
def main():
    argparser = argparse.ArgumentParser(
        description="Readmapper",
        usage="\n\treadmap -p genome\n\treadmap -d dist genome reads"
    )
    argparser.add_argument(
        "-p", action="store_true",
        help="preprocess the genome."
    )
    argparser.add_argument(
        "-u", help="Path to dir containing preprocessed genome",
    )
    argparser.add_argument(
        "-d", type=int, metavar="integer",
        default=1, help="max edit distance."
    )
    argparser.add_argument(
        "genome",
        help="Simple-FASTA file containing the genome.",
        type=argparse.FileType('r')
    )
    argparser.add_argument(
        "reads", nargs="?",
        help="Simple-FASTQ file containing the reads.",
        type=argparse.FileType('r')
    )
    args = argparser.parse_args()

    if args.p:
        #print(f"Preprocess {args.genome}")
        fasta_recs = read_fasta(args.genome)
        for fa_rec in fasta_recs:
            ref = fa_rec[1]
            SA = SuffixArray(ref)
            try: genome_name = args.genome.name.split('/')[-1]
            except: genome_name = args.genome.name
            write_SA('{}'.format(genome_name),'Preprocessed_{}'.format(fa_rec[0]), SA)
    
    elif args.u:
        #print(f"Use preprocessed genome {args.u}")
        fasta_recs = read_fasta(args.genome)
        fastq_recs = read_fastq(args.reads)
        for fa_rec in fasta_recs:
            ref = fa_rec[1]
            SA = open_SA(args.u, fa_rec[0])
            for fq_rec in fastq_recs:
                read = fq_rec[1]
                matches = approx_search(ref, read, SA, args.d)
                for match in matches:
                    read_name = fq_rec[0]
                    read_seq = fq_rec[1]
                    edits = get_edits(read_seq, fa_rec[1][match:match+len(fq_rec[1])])
                    cigar = edits_to_cigar(edits[2])
                    output = [read_name,fa_rec[0],str(match+1),cigar,read_seq]
                    print('\t'.join(output))

    else:
        if args.reads is None:
            argparser.print_help()
            sys.exit(1)
        else:
            print(f"Search {args.genome} for {args.reads} within distance {args.d}")
            fasta_recs = read_fasta(args.genome)
            fastq_recs = read_fastq(args.reads)
            for fa_rec in fasta_recs:
                ref = fa_rec[1]
                SA = SuffixArray(ref)
                for fq_rec in fastq_recs:
                    read = fq_rec[1]
                    matches = approx_search(ref, read, SA, args.d)
                    for match in matches:
                        read_name = fq_rec[0]
                        read_seq = fq_rec[1]
                        edits = get_edits(read_seq, fa_rec[1][match:match+len(fq_rec[1])])
                        cigar = edits_to_cigar(edits[2])
                        output = [read_name,fa_rec[0],str(match+1),cigar,read_seq]
                        print('\t'.join(output))

        

#########################################
# Code:

# (e.g. 'python3 readmap.py path_to_fasta path_to_fastq -d 2')
if __name__ == '__main__':
    main()






