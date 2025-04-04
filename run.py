import argparse
from scan_test import scan_test
from read_data import map_to_protein
from annotation_test import annotation_test


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument(
        '--rvas-data',
        type=str,
        help='.tsv.gz file with columns chr, pos, ref, alt, ac_case, ac_control',
    )
    parser.add_argument(
        '--scan-test',
        action='store_true',
        default=False,
        help = 'perform scan test',
    )
    parser.add_argument(
        '--clinvar-test',
        action='store_true',
        default=False,
        help='perform clinvar test',
    )
    parser.add_argument(
        '--annotation_file',
        type=str,
        default=None,
        help='perform a burden test using this annotation',
    )
    parser.add_argument(
        '--neighborhood-radius',
        type=float,
        default=15.0,
        help='neighborhood radius for clinvar or annotation tests',
    )
    parser.add_argument(
        '--filter-file',
        type=str,
        default=None,
        help='file to filter variants after expanding neighborhood'
    )
    parser.add_argument(
        '--genome-build',
        type=str,
        default='hg38',
        help='genome build. must be hg38 or hg37',
    )
    parser.add_argument(
        '--which-proteins',
        type=str,
        default='all',
        help='name of a protein or file with list of proteins'
    )
    parser.add_argument(
        '--reference-dir',
        type=str,
        help='directory with reference files'
    )
    args = parser.parse_args

    # map rvas results onto protein coordinates, linked to pdb files
    df_rvas = map_to_protein(args.rvas_data, args.which_proteins, args.genome_build)

    if args.scan_test:
        scan_test(df_rvas, args.reference_dir, args.neighborhood_radius)
    
    elif args.clinvar_test:
        annotation_file = f'{args.reference_dir}/ClinVar_PLP_uniprot_canonical.tsv.gz'
        filter_file = f'{args.reference_dir}/AlphaMissense_gt_0.9.tsv.gz'
        annotation_test(
            df_rvas,
            annotation_file,
            args.reference_dir,
            args.neighborhood_radius,
            filter_file,
        )

    elif args.annotation_file is not None:
        annotation_test(
            df_rvas,
            args.annotation_file,
            args.reference_dir,
            args.neighborhood_radius,
            args.filter_file,
        )
    
    else:
        print('no analysis specified')