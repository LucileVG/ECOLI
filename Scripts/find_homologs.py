import glob
import os

K12_folder = "/home/IAME/usearch/K12/*"
K12_files = glob.glob(K12_folder)
output_path="/home/IAME/usearch/output_files_b6/"

for query in K12_files:
    os.system("vsearch --usearch_global "+query+" --db /home/IAME/usearch/db.udb --id 0.90 --maxaccepts 70000 --blast6out "+output_path+query.split("/")[-1].strip(".fasta")+".b6")
