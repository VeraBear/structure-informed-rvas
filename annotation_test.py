def annotation_test(
        df_rvas,
        annotation_file,
        reference_directory,
        neighborhood_radius,
        filter_file,
    ):
    '''
    perform annotation test. annotation file and filter file have columns uniprot_id,
    aa_pos, aa_ref, aa_alt, which specify the members of the annotation/filter. 
    reference_directory has pdb_files. 

    this function loops over proteins. for each protein, it takes the annotation, uses the 
    pdb files to extend by the neighborhood radius, then filters using the filter file. then 
    performs fisher's exact to compare the resulting set of variants to the background of the 
    whole protein.
    '''