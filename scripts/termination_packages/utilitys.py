import pyranges as pr
import pandas as pd
import pysam
import numpy as np
import pickle
from scipy.stats import mannwhitneyu

def get_pas_info(pas_bed_in):
    pas_bed = pas_bed_in
    pas_bed = pr.read_bed(pas_bed, as_df=True)

    pas_bed.columns = ['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand', 'ratio']
    pas_bed.loc[:, 'Chromosome'] = pas_bed.loc[:, 'Chromosome'].astype('str')

    mask = pas_bed.loc[:, 'Name'].str.contains('_1')
    len(mask), len(pas_bed[mask])  ### mask is a list of bool

    single_pa_site_gene = pas_bed[mask].loc[:, 'Name'].map(lambda x: x.split('_')[0])

    pas_bed['Name'] = pas_bed['Name'].map(lambda x: x.split('_')[0])

    pas_info = pas_bed.apply(
        lambda x: (x.Chromosome, x.Start, x.Strand, x.Name)
        if x.Strand == '+' else
        (x.Chromosome, x.End, x.Strand, x.Name),
        axis=1)
    pas_info = np.array(list(pas_info))
    # print(pas_info[:3])
    print("Loading successfully!")
    return pas_info, pas_bed, single_pa_site_gene

def lp(pkl):
    with open(pkl, "rb") as fh:
        data = pickle.load(fh)
    return data

def get_header_of_bam(bamin:str):
    with pysam.AlignmentFile(bamin, 'rb') as inbam:
        header_str = inbam.header
    return str(header_str)

def termination_window(df, all_gene_in):
    df = df[df['gene_id'].isin(all_gene_in)]
    return df

def cleaved_ratio(df, all_gene_in):
    df = df.eval('ratio = (`5_cleavage`)/(readthrough+`5_cleavage`)')
    df = df[df['gene_id'].isin(all_gene_in)]
    return df

def cleavage_and_tw(ratios:list,
                    tws:list,
                    all_gene_in:set,
                    key:str = 'data'):
    wt_ratio, mutant_ratio = ratios
    _cl = pd.merge(wt_ratio, 
                   mutant_ratio, 
                   on='gene_id')
    _cl.drop(['readthrough_x','5_cleavage_x','3_cleavage_x','readthrough_y','5_cleavage_y','3_cleavage_y'], axis=1, inplace=True)
    _cl.dropna(inplace=True)
    _cl['log2FoldChange'] = np.log2(_cl['ratio_y']/_cl['ratio_x'])
    _cl = _cl[_cl.replace([np.inf, -np.inf], np.nan).notnull().all(axis=1)] 
    _cl.set_index(["gene_id"], inplace=True)
    
    wt_tw, mutant_tw = tws
    _tw = pd.merge(wt_tw, mutant_tw, 
                   on='gene_id')
    _tw.dropna(inplace=True)
    _tw['p_values'] = _tw.apply(lambda x: mannwhitneyu(
        x[f'{key}_x'], x[f'{key}_y'])[1], axis=1)
    _tw['move'] = _tw['median_y'] - _tw['median_x']
    _tw = _tw[_tw.replace([np.inf, -np.inf], np.nan).notnull().all(axis=1)] 
    
    _result = pd.concat([_cl, _tw],
                        axis=1)
    _result.dropna(inplace=True)
    _result = _result[_result.index.isin(all_gene_in)]
    
    return _cl, _result

def gnf(data, fo):
    with open(fo, "w") as fh:
        fh.write(data)
    print(f'>>{fo} had generated.<<')
