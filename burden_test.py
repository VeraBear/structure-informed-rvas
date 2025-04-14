import numpy as np
from pygam import LogisticGAM, LinearGAM, s, f, te
import statsmodels.formula.api as smf
import statsmodels.api as sm

def get_ypred(df, vep_list):
    #df has locus, is_case, and vep columns for all vep in vep_list
    df = df[~np.any(df[vep_list].isna(), axis=1)].copy()
    locus_list = np.unique(df.locus)
    n_loci = len(locus_list)
    df['ypred'] = -1.0
    for i, l in enumerate(locus_list):
        if i%100 == 0:
            print(f'locus {i} out of {n_loci}')
        f = df.locus == l
        df2 = df[~f]

        ytrain = df2['is_case'].values - df2['null_prob']
        
        aic = {}
        models = {} 
        for vep in vep_list:
            xtrain = df2[[vep]].values
            models[vep] = LinearGAM(s(0), verbose=False, lam=1000).fit(xtrain, ytrain)
            aic[vep] = models[vep].statistics_['AIC']
        vep = min(aic, key=aic.get)
        model = models[vep]
        xpred = df.loc[f, vep].values.reshape(-1, 1)
        ypred = model.predict(xpred)
        df.loc[f, 'ypred'] = ypred + df.loc[f, 'null_prob']
    return df

def burden_regression(df, column, filter=None):
    # df must have is_case, null_logit, and <column> as columns
    if filter is not None:
        data = df[filter]
    else:
        data = df
        
    model = smf.glm(
        f'is_case ~ {column} - 1',
        data=data,
        offset=data.null_logit,
        family=sm.families.Binomial(),
    ).fit(disp=False)
    if model.params[column] < 0:
        pval = 1
    else:
        pval = model.pvalues[column] 
    return pval

def burden_test_ypred(df_rvas, vep_list):
    uniprot_id_list = np.unique(df_rvas.uniprot_id)
    n_proteins = len(uniprot_id_list)
    for i, uniprot_id in enumerate(uniprot_id_list):
        print('\n', uniprot_id, f'number {i} out of {n_proteins}')
        try:
            df = df_rvas[df_rvas.uniprot_id == uniprot_id]
            df = get_ypred(df, vep_list)
            df['logit_ypred_minus_offset'] = np.log(df['ypred']/(1-df['ypred'])) - df['null_logit']
            burden_regression(df, 'logit_ypred_minus_offset')
        except Exception as e:
            print(f'Error for {uniprot_id}: {e}')
            continue