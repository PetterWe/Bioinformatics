#Python for genomics
import Bio
import numpy
import sys

try:
    f = open(r"C:\Users\pette\Documents\Python\dna2.fasta")
except IOError:
    print("oops wrong location!")

seqs={}
for line in f:
    line=line.rstrip() #rstrip removes newline characters like /n
    if line[0] =='>':#FASTA files always starts a new sequence with ">"" in the name line
        words=line.split()#splitting the test in the first line into words separated by space
        name=words[0][1:]#the name of the sequence is the first word, but we do not want the ">" so the slice the [0]:th word with the [1:] which indicate that we start with the second character
        seqs[name]=''#new entry into the dictionary
    else:#when we are not at the header line we do this, i.e. append (+ is append on strings) the sequence to seqs
        seqs[name]=seqs[name]+line

f.close()

#How many records are there in the file?
print("Headers:")
i=0
for name, seq in seqs.items():
    print(name)
    i=i+1

print("Number of headers:")
print(i)
print("\n")

#How long are the sequences?
print("How long are the sequences?")
from Bio import SeqIO
full_seq=[]
for seq_record in SeqIO.parse(r"C:\Users\pette\Documents\Python\dna2.fasta", "fasta"):
    full_seq.append([seq_record.id, str(seq_record.seq), len(seq_record)])
    print("Length of record",seq_record.id)
    print(len(seq_record))
    print("\n")
print("\n")
#Open reading frames
print("Open reading frames, writing to file")
print("\n")

def find_ORF(sequence, read_frame):
    longest_orf_length=0
    for i in range(read_frame,len(sequence),3):
        start_codon=sequence[i:i+3] #Reading 3 nucleotides from i
        if start_codon == "ATG": #Only ATG is a start codon
            start_position=i #remember which position we start at
            for j in range(start_position, len(sequence),3):
                stop_codon=sequence[j:j+3]
                if stop_codon in ["TAA", "TAG", "TGA"]:
                    stop_position=j
                    orf_length=(stop_position-start_position)+3
                    print(orf_length)
                    if orf_length > longest_orf_length:
                        longest_orf_length=orf_length
                        longest_start=start_position
                        longest_stop=stop_position


original_stdout = sys.stdout


for seq_record in SeqIO.parse(r"C:\Users\pette\Documents\Python\dna2.fasta", "fasta"):
    with open(r"C:\Users\pette\Documents\Python\orfs.txt","a") as f:
        sys.stdout = f
        print(seq_record.id, file=f)
        sequence=str(seq_record.seq)
        read_frame=2
        print(find_ORF(sequence,read_frame),file=f)


sys.stdout = original_stdout

#Finding repeats
print("Finding the repeats and writing to file")
from collections import Counter
def get_repeated_strings(input_string, min_str_length, calculate_largest_repeated_string = True ):

    all_substrings = [input_string[start_index:][:end_index + 1]
                      for start_index in range(len(input_string))
                      for end_index in range(len(input_string[start_index:]))]
    counted_substrings = Counter(all_substrings)
    not_counted_final_candidates = [item[0]
                                    for item in counted_substrings.most_common()
                                    if item[1] > 1 and len(item[0])==min_str_length]
    counted_final_candidates = {item: counted_substrings[item] for item in not_counted_final_candidates}
    print(counted_final_candidates)

for seq_record in SeqIO.parse(r"C:\Users\pette\Documents\Python\dna2.fasta", "fasta"):
    with open(r"C:\Users\pette\Documents\Python\repeats.txt","a") as f:
        sys.stdout = f
        print(seq_record.id, file=f)
        sequence=str(seq_record.seq)
        length=6
        print(get_repeated_strings(sequence,length),file=f)
sys.stdout = original_stdout

print("All done!")
