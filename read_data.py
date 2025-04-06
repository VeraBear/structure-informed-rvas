def map_to_protein(rvas_data, which_proteins, genome_build, reference_directory):
    '''
    rvas_data is a path to a .tsv.gz file with columns chr, pos, ref, alt, ac_case, ac_control.
    which_proteins is either the name of a protein or file with a list of proteins. we could 
    also make a mapper from gene name to uniprot id and allow this to be a gene name or file with
    multiple gene names.

    the output of this function is a dataframe where the columns are uniprot_id, aa_pos, aa_ref,
    aa_alt, pdb_file, file_index, ac_case, and ac_control

    this function should drop any proteins with insufficient data, and only output data for the
    requested proteins.
    '''

