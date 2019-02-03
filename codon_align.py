import subprocess
from Bio import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna

in_file = "/home/ubuntu/Python2/not_codon_aligned.fasta"

#mafft align
#subprocess.call(["mafft-fftns", "--adjustdirection", 
#                "--thread",  "2","--ep", "2", "--op", "3", "--out",
#                "nuc_aligned.fasta", in_file])
subprocess.call(["mafft", "--out", "nuc_aligned.fasta", in_file])

seqs=[]

#add sequences with initial gaps replaced by n's
for seq_record in SeqIO.parse("nuc_aligned.fasta", "fasta"):

    i=0
    while seq_record[i]=='-':
        i+=1
    ns="n"*i
    sequence=Seq.Seq(str(seq_record.seq).replace("-", ""))

    seqs.append(ns+sequence)

#trim to multiples of three 
for ind in range(len(seqs)):
    if len(seqs[ind])%3==2:
        temp=seqs[ind][0:-2]
        seqs[ind]=temp

    elif len(seqs[ind])%3==1:
        temp=seqs[ind][0:-1]
        seqs[ind]=temp

#translate to amino acids
aa_seqs=[]
for seq in seqs:
    aa_seqs.append(seq.translate(table=2))

#cumbersome type conversion
record_arr=[]
for seq in aa_seqs:
    record_arr.append(SeqIO.SeqRecord(seq))

#write to file and mafft align the amino acids
SeqIO.write(record_arr, "not_aligned_aas.fasta", "fasta")
subprocess.call(["mafft", "--out", "aligned_aas.fasta", "not_aligned_aas.fasta"])

#read in aligned amino acids
aligned_aas=[]
for seq_record in SeqIO.parse("aligned_aas.fasta", "fasta"):
    sequence=Seq.Seq(str(seq_record.seq))
    aligned_aas.append(sequence)

print(aligned_aas)
