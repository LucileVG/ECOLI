import glob
import os
cluster_path='/home/IAME/annotated_genes/*.fasta'
output_path='/home/IAME/annotated_genes_aln/'
cluster_file=glob.glob(cluster_path)
for f in cluster_file:
    outf_name=f.split('/')[-1].rstrip('.fasta')+'_aln.fasta'
    os.system('mafft --retree 1 --thread 16 '+f+' > '+output_path+outf_name)