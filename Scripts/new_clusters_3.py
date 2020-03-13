import glob
from Bio import SeqIO
from Bio.Seq import Seq


gap_folder = "/home/IAME/gaps/*"
new_clusters_folder = "/home/IAME/new_clusters/"
clusters_folder = "/home/SEYDI/cluster_aln_files_step7/"

for gap_file in glob.glob(gap_folder):
    with open(gap_file, "r") as dr:
        cluster = (gap_file.split('/')[-1]).strip("_aln.fasta.out")
        record = {record.id:record for record in SeqIO.parse(clusters_folder+cluster+"_aln.fasta", "fasta")}
        dr.readline()
        refseq = dict()
        sequence = dict()
        for line in dr:
            values = line.split(";")
            if(tuple(values[1:]) not in refseq.keys()):
                refseq[tuple(values[1:])] = values[0]
                sequence[values[0]] = dict()
            aln = record[values[0]]
            aln.seq = Seq(str(aln.seq).replace("-",""))
            sequence[refseq[tuple(values[1:])]][values[0]] = aln
        for key in sequence.keys():
            SeqIO.write(list(sequence[key].values()),new_clusters_folder+key+".fasta","fasta")
