from clinvar_af3 import process_chromosome, get_interactors, create_af3_jobs, get_interfaces, create_annotation_file_per_chr, clean_result, perform_fdr_correction
import argparse
import pandas as pd

if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument(
        '--protein-sequences',
        type=str,
        default=None,
        help='Provide the path to the file with protein sequences in tsv format with columns uniprot_id and sequence.'
    )
    parser.add_argument(
        '--create-raw-clinvar-rvas-files',
        type=str,
        default=None,
        help='create rvas file for each chromosome for clinvar data for the proteins we care. Provide the directory where to store the rvas files.'
    )
    parser.add_argument(
        '--raw-clinvar-dir',
        type=str,
        default='raw_ClinVar',
        help='directory that contains the raw clinvar files. Provide the directory where the raw clinvar files are stored.'
    )
    parser.add_argument(
        '--get-bioplex-interacting-partners',
        type=str,
        default=None,
        help='get interacting partners for each protein in the raw clinvar files using bioplex. Provide the directory where the mapped rvas files is stored.'
    )
    parser.add_argument(
        '--create-af3-jobs',
        type=str,
        default=None,
        help='create af3 jobs for each protein and its interacting partners. Provide the directory where the af3 jobs will be stored.'
    )
    parser.add_argument(
        '--get-interfaces-from-af3-results',
        type=str,
        nargs=2,
        metavar=('af3_result_dir', 'input_dir'),
        default=None,
        help='get interfaces from af3 results for each protein and its interacting partners and an annotation file that contains all chr. Provide the directory where the af3 results are stored and the directory where the combined annotation file and interface file will be stored.'
    )
    parser.add_argument(
        '--create-annotation-file-by-chromosome',
        type=str,
        default=None,
        help='create annotation file by chromosome. Provide the directory where the annotation files by chromosome will be stored.'
    )
    
    parser.add_argument(
        '--reference-dir',
        type=str,
        default=None,
        help='Provide the directory where all input files are stored.'
    )
    parser.add_argument(
        '--clean-result',
        type=str,
        nargs=2,
        metavar=('result_dir', 'input_dir'),
        default=None,
        help='aggregate and clean the results and correct for the fdr. Provide the directory where the preliminary results are stored.'
    )

    args = parser.parse_args()
    chromosomes = [str(i) for i in range(1, 23)] + ['X']
    uniprot_seq_df = pd.read_csv(f'{args.reference_dir}/{args.protein_sequences}', sep='\t') 
    uniprot_seq_dict = uniprot_seq_df.set_index('uniprot_id')['sequence'].to_dict()

    if args.create_raw_clinvar_rvas_files:
        df_map = pd.read_csv(f'{args.reference_dir}/all_missense_variants_gr38.tsv.gz', sep='\t', compression='gzip')
        for chrom in chromosomes:
            process_chromosome(chrom, df_map, uniprot_seq_dict, args.raw_clinvar_dir, args.create_raw_clinvar_rvas_files, args.reference_dir)
    
    if args.get_bioplex_interacting_partners and args.create_af3_jobs:
        BioPlex_network = pd.read_csv(f'{args.reference_dir}/BioPlex_293T_Network_10K_Dec_2019.tsv', sep='\t')
        all_df = get_interactors(BioPlex_network, args.get_bioplex_interacting_partners, args.reference_dir, uniprot_seq_dict, min_count=0)
        print(all_df.shape)
        create_af3_jobs(args.create_af3_jobs, all_df)

    if args.get_interfaces_from_af3_results:
        af3_out_main, out_dir = args.get_interfaces_from_af3_results
        get_interfaces(af3_out_main, out_dir, threshold=0.5)

    if args.create_annotation_file_by_chromosome:
        create_annotation_file_per_chr(args.create_annotation_file_by_chromosome)
    
    if args.clean_result:
        result_dir, input_dir = args.clean_result
        clean_result(result_dir, input_dir, args.reference_dir)
    
    

