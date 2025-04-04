def pymol_rvas(df_rvas, uniprot_id, reference_directory):
    # make a pymol session with case and control mutations
    # output a gif and a .pse file
    # df_rvas columns: uniprot_id, aa_pos, aa_ref, aa_alt, pdb_file, file_index, ac_case, and ac_control

def pymol_annotation(annot_file, uniprot_id, reference_directory):
    # visualize the annotation

def pymol_scan_test(df_results, uniprot_id, reference_directory):
    # color by case/control ratio of the neighborhood

def pymol_neighborhood(df_results, uniprot_id, reference_directory):
    # for each significant neighborhood, zoom in and show the case and control mutations
    # just in that neighborhood.