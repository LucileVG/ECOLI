new_clusters_folder = "/home/IAME/new_clusters/"
annotated_genes_folder = "/home/IAME/annotated_genes/"
annotations_csv = "/home/IAME/full_sequences_no_duplicates_test.csv"

with open(annotations_csv, "r") as dr:
    dr.readline()
    line_dr = dr.readline()
    values = line_dr.split(";")
    gene = values[0]
    reference_sequence = [values[4]]
    for line_dr in dr:
        values = line_dr.split(";")
        if(values[0]!=gene):
            with open(annotated_genes_folder+gene+".fasta", "w") as dw:
                for file_name in reference_sequence:
                    with open(new_clusters_folder+file_name+".fasta", "w") as df:
                        for line_df in df:
                            dw.write(line_df)
            gene = values[0]
            reference_sequence = [values[4]]
        else:
            reference_sequence.append(values[4])



