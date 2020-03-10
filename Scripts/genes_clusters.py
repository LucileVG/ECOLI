import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

path_1 = "/home/SEYDI/gene_cluster_files_step6/singleton_seq_files"
path_2 = "/Users/Lucile/Documents/IAME/ECOLI/Data/Sequences/HPI sequences/Full sequences"

MAX_GAPS_PROP = 1/180
MAX_GAPS_ABS = 5

with open("/Users/Lucile/Documents/IAME/ECOLI/Data/HPI outputs/Unique_HPI_hits.csv", "r") as df:
    df.readline()
    for line in df:
        row = line.strip().split(",")
        gene_id = row[0]
        cluster_id = row[1]
        gene_length = int(row[2])
        fragment = (row[3]=="True")
        singleton = (row[4]=="True")
        bad_algnt = (row[5]=="True")
        print("Processing gene "+gene_id+"/cluster "+cluster_id)
        if(singleton):
            folder = "/Users/Lucile/Documents/IAME/ECOLI/Data/Sequences/HPI sequences/presumed_genes/full_sequence/"
            if(fragment):
                folder = "/Users/Lucile/Documents/IAME/ECOLI/Data/Sequences/HPI sequences/presumed_genes/fragment/"
            with open(folder+"gene_"+gene_id+".fasta", "a") as write_file:
                with open(path_1+"/cluster"+cluster_id+".fasta", "r") as read_file:
                    for line in read_file:
                        write_file.write(line)
        else:
            #clusters_seq = {record.id:record for record in SeqIO.parse(path_2+"/cluster"+cluster_id+"_aln.fasta", "fasta")}
            clusters_seq = SeqIO.to_dict(SeqIO.parse(path_2+"/cluster"+cluster_id+"_aln.fasta", "fasta"))
            characteristics_seq = dict()
            full_algnt_seq = dict()
            fragments_seq = dict()
            others_seq = dict()
            for key in clusters_seq.keys():
                sequence_dash = [int(i=="-") for i in str(clusters_seq[key].seq)]
                a = sequence_dash+[0]
                b = [0]+sequence_dash
                c = [int((a[i]-b[i])==1) for i in range(len(a))]
                gaps = np.sum(np.array(c))
                dashes = np.sum(np.array(sequence_dash))
                characters = len(sequence_dash)-np.sum(np.array(sequence_dash))
                characteristics_seq[key] = [characters,dashes,gaps]
            if(bad_algnt):
                gaps = int(row[6])
                dashes = int(row[7])
                for key in clusters_seq.keys():
                    [characters_seq,dashes_seq,gaps_seq] = characteristics_seq[key]
                    if(dashes_seq>0.9*dashes and dashes_seq<1.1*dashes and gaps_seq>0.9*gaps and gaps_seq<1.1*gaps):
                        sequence = clusters_seq[key]
                        sequence.seq = Seq((str(clusters_seq[key].seq)).replace('-',''))
                        if(fragment):
                            fragments_seq[key] = sequence
                        else:
                            full_algnt_seq[key] = sequence
                    else:
                        sequence = clusters_seq[key]
                        sequence.seq = Seq((str(clusters_seq[key].seq)).replace('-',''))
                        others_seq[key] = sequence
            else:
                for key in clusters_seq.keys():
                    [characters_seq,dashes_seq,gaps_seq] = characteristics_seq[key]
                    sequence = clusters_seq[key]
                    len_seq = len(sequence.seq)
                    sequence.seq = Seq((str(clusters_seq[key].seq)).replace('-',''))
                    nb_characters = len(sequence.seq)
                    if(gaps_seq>max(MAX_GAPS_PROP*len_seq,MAX_GAPS_ABS) and dashes_seq>0.1*len_seq):
                        others_seq[key] = sequence
                    elif(nb_characters<0.9*gene_length):
                        fragments_seq[key] = sequence
                    elif(nb_characters<1.1*gene_length):
                        full_algnt_seq[key] = sequence
                    else:
                        others_seq[key] = sequence
            if(len(others_seq.keys())>0):
                SeqIO.write(list(others_seq.values()),"/Users/Lucile/Documents/IAME/ECOLI/Data/Sequences/HPI sequences/presumed_genes/temp.fasta","fasta")
                folder = "/Users/Lucile/Documents/IAME/ECOLI/Data/Sequences/HPI sequences/presumed_genes/other/"
                with open(folder+"gene_"+gene_id+".fasta", "a") as write_file:
                    with open("/Users/Lucile/Documents/IAME/ECOLI/Data/Sequences/HPI sequences/presumed_genes/temp.fasta", "r") as read_file:
                        for line in read_file:
                            write_file.write(line)
            if(len(full_algnt_seq.keys())>0):
                SeqIO.write(list(full_algnt_seq.values()),"/Users/Lucile/Documents/IAME/ECOLI/Data/Sequences/HPI sequences/presumed_genes/temp.fasta","fasta")
                folder = "/Users/Lucile/Documents/IAME/ECOLI/Data/Sequences/HPI sequences/presumed_genes/full_sequence/"
                with open(folder+"gene_"+gene_id+".fasta", "a") as write_file:
                    with open("/Users/Lucile/Documents/IAME/ECOLI/Data/Sequences/HPI sequences/presumed_genes/temp.fasta", "r") as read_file:
                        for line in read_file:
                            write_file.write(line)
            if(len(fragments_seq.keys())>0):
                SeqIO.write(list(fragments_seq.values()),"/Users/Lucile/Documents/IAME/ECOLI/Data/Sequences/HPI sequences/presumed_genes/temp.fasta","fasta")
                folder = "/Users/Lucile/Documents/IAME/ECOLI/Data/Sequences/HPI sequences/presumed_genes/fragment/"
                with open(folder+"gene_"+gene_id+".fasta", "a") as write_file:
                    with open("/Users/Lucile/Documents/IAME/ECOLI/Data/Sequences/HPI sequences/presumed_genes/temp.fasta", "r") as read_file:
                        for line in read_file:
                            write_file.write(line)