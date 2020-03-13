from Bio import SeqIO
from Bio.Seq import Seq


clusters_folder='/home/SEYDI/cluster_aln_files_step7/'
output_folder = "/home/IAME/usearch/"
output_fasta = output_folder+"db.fasta"
input_fasta = "/home/IAME/test.csv"


with open(output_fasta, "a") as dw:
    with open(input_fasta, "r") as dr:
        dr.readline()
        line = dr.readline()
        values = line.split(";")
        names = [values[1].strip()]
        cluster = values[0]
        for line in dr:
            values = line.split(";")
            if(cluster!=values[0]):
                record = {record.id:record for record in SeqIO.parse(clusters_folder+cluster+"_aln.fasta", "fasta")}
                sequences = dict()
                for seqname in names:
                    sequence = record[seqname]
                    sequence.seq = Seq(str(sequence.seq).replace("-",""))
                    sequences[seqname] = sequence
                SeqIO.write(list(sequences.values()),output_folder+"temp.fasta","fasta")
                with open(output_folder+"temp.fasta", "r") as dt:
                    for line_temp in dt:
                        dw.write(line_temp)
                names = [values[1].strip()]
                cluster = values[0]
            else:
                names.append(values[1].strip())