from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

#Get the sequence hits from NCBI
fasta_string = open(r"C:\Users\pette\Documents\Python\unknownseq.txt").read()
result_handle=NCBIWWW.qblast("blastn","nt",fasta_string)

#Pretty print them
blast_record = NCBIXML.read(result_handle)
E_VALUE_THRESH = 0.01

for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print(alignment.title)

            