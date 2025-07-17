import os
import json
import pandas as pd
import numpy as np
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
import statsmodels.stats.multitest as multitest

def fetch_sequence(uniprot_id: str, session: requests.Session):
    try:
        resp = session.get(f'https://rest.uniprot.org/uniprotkb/{uniprot_id}', timeout=10)
        resp.raise_for_status()
        seq = resp.json().get('sequence', {}).get('value', None)
        return uniprot_id, seq
    except Exception:
        return uniprot_id, None

def process_chromosome(chrom: str, df_map, uniprot_seq_dict, raw_clinvar_dir, out_dir, reference_dir):

    if chrom == 'X':
        clinvar_path = f'{raw_clinvar_dir}/clinvar_grch38_annotated_2_chrX.tsv.gz'
    else:
        clinvar_path = f'{raw_clinvar_dir}/clinvar_grch38_annotated_2_chr{chrom}.tsv.gz'

    if not os.path.exists(clinvar_path):
        print(f"ClinVar file for chromosome {chrom} does not exist. Skipping.")
        return

    df = pd.read_csv(clinvar_path, sep='\t', compression='gzip')

    df = df[~df['clinical_significance'].str.contains('Uncertain|Conflicting', na=False)]
    df['clinical_significance_group'] = (
        df['clinical_significance']
          .apply(lambda x: 'benign' if ('Likely benign' in x or 'Benign' in x)
                           else ('pathogenic' if ('Pathogenic' in x or 'Likely pathogenic' in x)
                                 else 'other'))
    )

    tmp_df_map = df_map[df_map['chrom'] == f'chr{chrom}']
    df_mapped = pd.merge(df, tmp_df_map, on=['pos', 'ref', 'alt'])
    

    df_mapped = df_mapped[df_mapped['uniprot_id'].notna()].copy()
    df_mapped['sequence'] = df_mapped['uniprot_id'].map(uniprot_seq_dict)

    missing_ids = set(df_mapped.loc[df_mapped['sequence'].isna(), 'uniprot_id'].unique()) 

    if missing_ids:
        with requests.Session() as session:
            with ThreadPoolExecutor(max_workers=10) as executor:
                futures = {
                    executor.submit(fetch_sequence, up_id, session): up_id
                    for up_id in missing_ids
                }
                for fut in as_completed(futures):
                    up_id = futures[fut]
                    uid, seq = fut.result()
                    if seq:
                        df_mapped.loc[df_mapped['uniprot_id'] == uid, 'sequence'] = seq

    num_variants_initial = len(df_mapped)
    num_variants_with_seq = df_mapped['sequence'].notna().sum()
    print(f"chr{chrom:>2}: {num_variants_initial=}  {num_variants_with_seq=}")
    df_mapped = df_mapped[df_mapped['clinical_significance_group'] != 'other']

    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(f'{out_dir}/uniprot_clinvar', exist_ok=True)
    df_mapped.to_csv(f'{out_dir}/uniprot_clinvar/clinvar_grch38_annotated_2_chr{chrom}_with_uniprot.tsv.gz', sep='\t', index=False, compression='gzip')

    df_mapped_rvas = df_mapped[['locus', 'alleles', 'clinvar_variation_id', 'pos', 'clinical_significance', 'clinical_significance_group', 'uniprot_id', 'variant_id']]
    df_mapped_rvas['ac_case'] = df_mapped_rvas['clinical_significance'].apply(
        lambda x: 1 if ('Pathogenic' in x or 'Likely pathogenic' in x) else 0)
    df_mapped_rvas['ac_control'] = df_mapped_rvas['clinical_significance'].apply(
        lambda x: 1 if ('Benign' in x or 'Likely benign' in x) else 0)

    os.makedirs(f'{out_dir}/uniprot_clinvar_rvas', exist_ok=True)
    df_mapped_rvas.to_csv(f'{out_dir}/uniprot_clinvar_rvas/clinvar_grch38_annotated_2_chr{chrom}_rvas.tsv.gz', sep='\t', index=False, compression='gzip')

def get_interactors(network_df, out_dir, reference_dir, uniprot_seq_dict, min_count=0):

    if not os.path.exists(f'{out_dir}/uniprot_clinvar'):
        raise FileNotFoundError(f"Run --create-raw-clinvar-rvas-files to map uniprot id to clinvar data.")

    tsv_fs = [f for f in os.listdir(f'{out_dir}/uniprot_clinvar') if f.endswith('.tsv.gz')]
    # print(tsv_fs)
    all_df = pd.DataFrame(columns=["chromosome", "UniprotA", "UniprotB", "pInt"])

    for f in tsv_fs:
        chrom = f.split('_')[-3]
        df_clean = pd.read_csv(os.path.join(f'{out_dir}/uniprot_clinvar', f), sep='\t', compression='gzip')
        df_clean = df_clean[df_clean['clinical_significance_group'] != 'other']
        # print(df_clean)

        df_counts = (
            df_clean.groupby(['uniprot_id', 'clinical_significance_group'])
            .size()
            .unstack(fill_value=0)
        )

        df_counts['min_count'] = df_counts.min(axis=1)
        df_counts = df_counts[df_counts['min_count'] > min_count].reset_index()

        for uniprot_id in df_counts['uniprot_id']:
            print(f'Processing {uniprot_id} for chromosome {chrom}')
            df_A = network_df[network_df['UniprotA'] == uniprot_id]
            df_B = network_df[network_df['UniprotB'] == uniprot_id]
            print(df_A, df_B)

            pairsA = [(uniprot_id, row['UniprotB'], row['pInt']) for _, row in df_A.iterrows()]
            pairsB = [(uniprot_id, row['UniprotA'], row['pInt']) for _, row in df_B.iterrows()]

            tmp_pairs = pairsA + pairsB
            df_tmp_pairs = pd.DataFrame(tmp_pairs, columns=["UniprotA", "UniprotB", "pInt"])
            df_tmp_pairs['chromosome'] = chrom
            df_tmp_pairs['pInt'] = df_tmp_pairs['pInt'].astype(float)
            print(df_tmp_pairs)
            df_tmp_pairs = df_tmp_pairs[df_tmp_pairs['pInt'] >= 0.9]
            if len(df_tmp_pairs) == 0:
                continue
            if len(df_tmp_pairs) > 10:
                df_tmp_pairs_10 = df_tmp_pairs.sort_values(by='pInt', ascending=False).iloc[:10]
            else:
                df_tmp_pairs_10 = df_tmp_pairs
            all_df = pd.concat([all_df, df_tmp_pairs_10], ignore_index=True)
        
        all_df_clean = all_df[all_df['UniprotB'] != 'UNKNOWN']
        all_df_clean.reset_index(drop=True, inplace=True)
        all_df_clean['SeqA'] = all_df_clean['UniprotA'].apply(lambda x: uniprot_seq_dict.get(x, np.nan))
        all_df_clean['SeqB'] = all_df_clean['UniprotB'].apply(lambda x: uniprot_seq_dict.get(x, np.nan))

        if not os.path.exists(f'{reference_dir}/uniprot_summary.csv'):
            raise FileNotFoundError(f"{reference_dir}/uniprot_summary.csv file not found. This is the form that has isoform sequences for uniprot ids. Please ensure it exists in the current directory.")
        isoform_seq = pd.read_csv(f'{reference_dir}/uniprot_summary.csv')
        isoform_seq_dict = isoform_seq.set_index('uniprot_isoform')['sequence'].to_dict()
        all_df_clean['SeqB'] = all_df_clean.apply(
            lambda x: isoform_seq_dict.get(x['UniprotB'], np.nan) if pd.isna(x['SeqB']) else x['SeqB'],
            axis=1
        )
        all_df_clean['SeqB'] = all_df_clean.apply(lambda x: fetch_sequence(x['UniprotB'], requests.Session()) if pd.isna(x['SeqB']) else x['SeqB'], axis=1)
        all_df_clean['SeqB'] = all_df_clean['SeqB'].apply(lambda x: x[-1] if str(x).startswith('(') else x)
        all_df_clean = all_df_clean[~all_df_clean['SeqB'].astype(str).str.contains('None', na=False)]
        print(all_df_clean)
    return all_df_clean

def create_af3_jobs(af3_job_dir, all_df_clean):
    os.makedirs(af3_job_dir, exist_ok=True)
    print(f'Creating AF3 jobs in {af3_job_dir} for {len(all_df_clean)} pairs.')
    for (chromosome, uniprotA), group_df in all_df_clean.groupby(['chromosome', 'UniprotA']):
        print(f'Group: chromosome={chromosome}, UniprotA={uniprotA}')
        print(group_df)
        jobs = []
        for _, row in group_df.iterrows():
            j = {"name": f"{row['chromosome']}_{row['UniprotA']}_{row['UniprotB']}", 
                "modelSeeds": [],
                "sequences":[
                    {'proteinChain':{
                        'sequence': row['SeqA'],
                        'count' : 1
                        }
                    }, 
                    {'proteinChain':{
                            'sequence': row['SeqB'],
                            'count' : 1
                        }
                    }
                ]}
            jobs.append(j)
        with open(f'{af3_job_dir}/{chromosome}_{uniprotA}.json', "w") as f:
            json.dump(jobs, f, indent=4)

def get_interfaces(af3_out_main, out_dir, threshold=0.5):
    '''
    af3_out_main: str, path to the main directory containing AF3 output files
    gene_to_uniprot: dict, mapping of gene names to UniProt IDs
    threshold: float, contact probability threshold to consider a residue as part of the interface (default is 0.5)
    '''

    af3_out_dir = [f'{af3_out_main}/{d}' for d in os.listdir(af3_out_main) if os.path.isdir(f'{af3_out_main}/{d}')]
    data_file_dirs = {}
    for item in af3_out_dir:
        data_file_dirs[item] = []
        for f in os.listdir(item):
            if 'full_data' in f:
                data_file_dirs[item].append(f'{item}/{f}')

    df_annot = pd.DataFrame(columns=['chrom', 'uniprot_id', 'aa_pos', 'interactor'])
    df_inter = pd.DataFrame(columns=['chrom', 'uniprot_id', 'interactor', 'aa_pos', 'max_contact_prob'])
    
    for k, v in data_file_dirs.items():
        print(k)
        for item in v:
            if not item.endswith('0.json'):
                continue
            gene = item.split('/')[-1].split('_')
 
            
            if gene[0] != 'fold':
                tmp_chrom = gene[0]
                tmp_uniprot = gene[1]
                tmp_interactor = gene[2]
            else:
                tmp_chrom = gene[1]
                tmp_uniprot = gene[2]
                tmp_interactor = gene[3]

            try:
                full_data = json.load(open(item))
            except Exception as error:
                print(error)
                continue
                
            token_chain_ids = np.array(full_data['token_chain_ids'])
            contact_prob = np.array(full_data['contact_probs'])
            chainA_inds_dv = np.where(token_chain_ids == 'A')[0].flatten()
            chainB_inds = np.where(token_chain_ids == 'B')[0].flatten()
            contact_prob = contact_prob[np.ix_(chainA_inds_dv, chainB_inds)]
            print(contact_prob.shape)
            tmp_max_prob = np.max(contact_prob)
            contact_prob_high = np.where(contact_prob>threshold)
            tmp_aa_pos = contact_prob_high[0] + 1
            tmp_aa_pos = list(set(tmp_aa_pos))
            

            if len(tmp_aa_pos) != 0:
                tmp_chroms = [tmp_chrom] * len(tmp_aa_pos)
                tmp_uniprot_id = [tmp_uniprot] * len(tmp_aa_pos)
                tmp_interactors = [tmp_interactor.upper()] * len(tmp_aa_pos)
            
            else:
                flat = contact_prob.flatten()
                sorted_indices = np.argsort(flat)[::-1]
                rows, cols = np.unravel_index(sorted_indices, contact_prob.shape)
                top_coords = list(zip(rows, cols))
                tmp_aa_pos = list(set([i[0] + 1 for i in top_coords[:5]]))
                tmp_chroms = [tmp_chrom] * len(tmp_aa_pos)
                tmp_uniprot_id = [tmp_uniprot] * len(tmp_aa_pos)
                tmp_interactors = [tmp_interactor.upper()] * len(tmp_aa_pos)
            
            for i in range(len(tmp_aa_pos)):
                tmp_row = {'chrom': tmp_chroms[i], 'uniprot_id': tmp_uniprot_id[i], 'aa_pos': tmp_aa_pos[i], 'interactor': tmp_interactors[i]}
                df_annot.loc[len(df_annot)] = tmp_row
        
        tmp_inter_row = {'chrom': tmp_chroms, 'uniprot_id': tmp_uniprot, 'interactor': tmp_interactor, 'aa_pos': tmp_aa_pos, 'max_contact_prob': tmp_max_prob}
        df_inter.loc[len(df_inter)] = tmp_inter_row
        df_annot['uniprot_id'] = df_annot['uniprot_id'].apply(lambda x: x.upper())
        df_annot['annotation_id'] = df_annot['uniprot_id'] + '_' + df_annot['interactor']
        os.makedirs(out_dir, exist_ok=True)
        df_annot.to_csv(f'{out_dir}/af3_interfaces_annot_combined.tsv', sep='\t', index=False)
        df_inter.to_csv(f'{out_dir}/af3_interfaces_combined.tsv', sep='\t', index=False)

    return df_annot, df_inter

def create_annotation_file_per_chr(out_dir):
    if not os.path.exists(out_dir):
        raise FileNotFoundError(f"Reference directory {out_dir} does not exist.")
    if not os.path.exists(f'{out_dir}/af3_interfaces_annot_combined.tsv'):
        raise FileNotFoundError(f"Run --get-interfaces-from-af3-results to create the combined annotation file.")
    
    df_annot = pd.read_csv(f'{out_dir}/af3_interfaces_annot_combined.tsv', sep='\t')
    print(df_annot)
    all_chr = set(df_annot['chrom'].tolist())
    all_chr = [int(item[3:]) for item in list(all_chr)]
    for i in range(1, max(all_chr)+1):
        df = df_annot[df_annot['chrom'] == f'chr{i}'].reset_index(drop=True)
        print(df)
        os.makedirs(f'{out_dir}/annotation_files_by_chr', exist_ok=True)
        df.to_csv(f'{out_dir}/annotation_files_by_chr/af3_interfaces_annot_chr{i}.tsv', index=False, sep='\t')
    
    df_annot_uni = set(df_annot['uniprot_id'].tolist())
    with open(f'{out_dir}/protein_list.txt', 'w') as f:
        for item in df_annot_uni:
            f.write(f'{item}\n')

def perform_fdr_correction(p):
    p = np.array(p)
    mask = np.isfinite(p)
    p_reject1, p_fdr1 = multitest.fdrcorrection(p[mask], alpha=0.05)
    p_fdr = np.full(p.shape, np.nan)
    p_fdr[mask] = p_fdr1
    p_reject = np.full(p.shape, False)
    p_reject[mask] = p_reject1
    print(len(p_fdr), len(p_reject))
    return p_fdr, p_reject

def clean_result(annot_test_dir, out_dir, reference_dir):
    annot_test_files = [f'{annot_test_dir}/{d}' for d in os.listdir(annot_test_dir) if os.path.isdir(f'{annot_test_dir}/{d}')]
    annot_file_dirs = {}

    for item in annot_test_files:
        annot_file_dirs[item] = []
        tsv_paths = [
            os.path.join(item, f)
            for f in os.listdir(item)
            if f.endswith('.tsv')
        ]
        if tsv_paths:
            newest = max(tsv_paths, key=os.path.getmtime)
            annot_file_dirs[item].append(newest)

    print(annot_file_dirs)
    results_pair_df = pd.DataFrame(columns=['chromosome', 'uniprot_id', 'annotation_id', 'in_case', 'out_case', 'in_control', 'out_control', 'or', 'p', 'p_fdr', 'fdr_reject'])
    
    for k, v in annot_file_dirs.items():
        print(k,v)
        chr = k.split('/')[-1].split('_')[-1]
        print(chr)
        for item in v:
            print(item)
            tmp_df = pd.read_csv(item, sep='\t')
            tmp_df = tmp_df.dropna()
            for i in range(len(tmp_df)):
                tmp_row = {
                    'chromosome': chr,
                    'uniprot_id': tmp_df['uniprot_id'][i],
                    'annotation_id': tmp_df['annotation_id'][i],
                    'in_case': tmp_df['in_case'][i],
                    'out_case': tmp_df['out_case'][i],
                    'in_control': tmp_df['in_control'][i],
                    'out_control': tmp_df['out_control'][i],
                    'or': tmp_df['or'][i],
                    'p': tmp_df['p'][i],
                    'p_fdr': tmp_df['p_fdr'][i],
                    'fdr_reject': tmp_df['fdr_reject'][i]
                }
                results_pair_df.loc[len(results_pair_df)] = tmp_row
    
    p_fdr, p_reject = perform_fdr_correction(results_pair_df['p'].values)
    results_pair_df['p_fdr_all_chr'] = p_fdr
    results_pair_df['fdr_reject_all_chr'] = p_reject

    gene_uni = pd.read_csv(f'{reference_dir}/gene_to_uniprot_id.tsv', sep='\t')
    af3_interface = pd.read_csv(f'{out_dir}/af3_interfaces_combined.tsv', sep='\t')

    # results_pair_df['protein'] = results_pair_df['uniprot_id'].apply(lambda x: x.split('-')[0])
    results_pair_df['interactor'] = results_pair_df['annotation_id'].apply(lambda x: x.split('_')[1])
    results_pair_df['gene_symbol'] = results_pair_df['uniprot_id'].apply(
        lambda x: gene_uni[gene_uni['uniprot_id'] == x]['gene_name'].values[0] 
        if (gene_uni['uniprot_id'] == x).any() else None
    )
    results_pair_df['interactor_gene_symbol'] = results_pair_df['interactor'].apply(
        lambda x: gene_uni[gene_uni['uniprot_id'] == x]['gene_name'].values[0] 
        if (gene_uni['uniprot_id'] == x).any() else None
    )

    af3_interface['uniprot_id'] = af3_interface['uniprot_id'].apply(lambda x: x.upper())
    af3_interface['interactor'] = af3_interface['interactor'].apply(lambda x: x.upper())

    results_pair_df['aa_pos'] = results_pair_df.apply(
    lambda x: af3_interface[
        (af3_interface['uniprot_id'] == x['uniprot_id']) &(af3_interface['interactor'] == x['interactor'])
    ]['aa_pos'].values[0]
    if not af3_interface[
        (af3_interface['uniprot_id'] == x['uniprot_id']) & (af3_interface['interactor'] == x['interactor'])
    ].empty else None, axis=1
    )

    results_pair_df['max_contact_prob'] = results_pair_df.apply(
    lambda x: af3_interface[
        (af3_interface['uniprot_id'] == x['uniprot_id']) &(af3_interface['interactor'] == x['interactor'])
    ]['max_contact_prob'].values[0]
    if not af3_interface[
        (af3_interface['uniprot_id'] == x['uniprot_id']) & (af3_interface['interactor'] == x['interactor'])
    ].empty else None, axis=1
    )
    results_pair_df.to_csv(f'{annot_test_dir}/clinvar_af3_results.tsv', sep='\t', index=False)

    return results_pair_df