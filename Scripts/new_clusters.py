import glob

output_fasta = "/home/lucile.vigue/IAME/test.csv"
gap_folder = "/home/lucile.vigue/IAME/gaps/*"

with open(output_fasta, "w") as dw:
    dw.write("Cluster;Sequence\n")
    for gap_file in glob.glob(gap_folder):
        with open(gap_file, "r") as dr:
            dr.readline()
            gaps_set = set()
            for line in dr:
                values = line.split(";")
                if(tuple(values[1:]) not in gaps_set):
                    gaps_set.add(tuple(values[1:]))
                    dw.write((gap_file.split('/')[-1]).strip("_aln.fasta.out")+";"+values[0]+"\n")
