#!/usr/bin/env python
# coding=utf-8

import numpy as np
import pyBigWig
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
import pysam
import scipy
import pandas as pd
import pyranges as pr

from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches

def last_pa_func(last_pa_bed):
    last_pa = pr.read_bed(last_pa_bed, as_df=True)
    last_pa.columns = ['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand', 'ratio']
    last_pa.loc[:, 'Chromosome'] = last_pa.loc[:, 'Chromosome'].astype('str')

    mask = last_pa.loc[:, 'Name'].str.contains('_1')
    single_pa_site_gene = last_pa[mask].loc[:, 'Name'].map(lambda x: x.split('_')[0])

    last_pa['Name'] = last_pa['Name'].map(lambda x: x.split('_')[0])
    last_pa = last_pa.set_index(['Name'])
    return last_pa


def gene_model_func(gene_model_bed):
    gene_model = pr.read_bed(gene_model_bed, as_df=True)
    gene_model = gene_model.set_index(['Name'])
    return gene_model


def downstream_df(infile):
    infile = '/public/home/mowp/db/Arabidopsis_thaliana/intergenic_region/araport11.distance_to_downstream_gene.bed'
    distance_to_downstream = pr.read_bed(infile, as_df=True)

    distance_to_downstream['upstream'] = distance_to_downstream['Name'].map(lambda x: x.split('_')[0])
    distance_to_downstream['downstream'] = distance_to_downstream['Name'].map(lambda x: x.split('_')[1])

    same_strand_gene = set(distance_to_downstream.query('Score == 1')['upstream'])
    different_strand_gene = set(distance_to_downstream.query('Score == 0')['upstream'])

    return distance_to_downstream


def get_target_site(site_type, gene_id, last_pa, gene_model, distance_to_downstream):
    try:
        if site_type == 'PAS':
            return last_pa.at[gene_id, 'End']
        elif site_type == 'TSS':
            try:
                values = tss_bed.loc[gene_id, :].values
                if values[4] == '+':
                    return values[1]
                else:
                    return values[2]
            except KeyError:
                values = gene_model.loc[gene_id, :].values
                if values[4] == '+':
                    return values[1]
                else:
                    return values[2]

        elif site_type == 'TES':
            values = gene_model.loc[gene_id, :].values
            if values[4] == '+':
                return values[2]
            else:
                return values[1]

        elif site_type == 'downstream':
            try:
                downstream_gene = distance_to_downstream.query('upstream == @gene_id')['downstream'].values[0]
                v = distance_to_downstream.query('upstream == @gene_id').values[0]
                start = (v[1]+v[2])//2
                # start = gene_model.loc[downstream_gene, :].values[1]
                values = gene_model.loc[downstream_gene, :].values
                if abs(values[1]-start) < abs(values[2]-start):
                    return values[1]
                else:
                    return values[2]
            except:
                return None

        else:
            raise KeyError
    except KeyError as e:
        print(f"KeyError: {str(e)}")


def get_three_end_pos(infile, gene_id, last_pa, gene_model, distance_to_downstream, before, after):
    
    STRAND_TO_BOOL = {'-': True, '+': False}
    
    chrom, *_, strand = gene_model.loc[gene_id]
    strand_boo = STRAND_TO_BOOL[strand]
    if chrom in {'Pt', 'Mt'}:
        return None
    
    n = 0
    cov_list = []
    read_set = set()

    target_site = get_target_site('PAS', gene_id, last_pa, gene_model, distance_to_downstream)
    downstream = get_target_site('downstream', gene_id, last_pa, gene_model, distance_to_downstream)
    if downstream is None:
        return
    
    three_end_pos = np.zeros(before+after)

    if strand == '+':
        start = target_site-before
        end = target_site+after
    else:
        start = target_site-after
        end = target_site+before
    if start <= 0:
        return

    with pysam.AlignmentFile(infile, 'rb') as inbam:
        for read in inbam.fetch(chrom, start, end):
            read_strand = read.is_reverse
            if strand_boo is not read_strand:
                continue
            if strand == '+':
                read_three_end = read.reference_end-target_site
                read_five_end = read.reference_start
                if read_five_end > downstream - 100:
                    continue
            else:
                read_three_end = target_site-read.reference_start
                read_five_end = read.reference_end
                if read_five_end < downstream + 100:
                    continue
            if -before <= read_three_end < after:
                three_end_pos[read_three_end+before] += 1

    if sum(three_end_pos) >= 1:  ### origin is 15
        return three_end_pos

def get_three_end_pos_turbo(infile, gene_list, 
                            last_pa, gene_model, distance_to_downstream, 
                            before=1000, after=1000, threads=32):
    results = []
    with ProcessPoolExecutor(max_workers=threads) as e:
        chunksize = int(len(gene_list)/threads)
        results = e.map(get_three_end_pos, repeat(infile), gene_list,
                        repeat(last_pa), repeat(gene_model),
                        repeat(distance_to_downstream),
                        repeat(before), repeat(after), chunksize=chunksize)
    n = 0
    three_end_pos = np.zeros(before+after)
    # _i = 0
    for res in results:
        if res is not None and sum(res) > 0:
            three_end_pos += res/sum(res)
            n += 1
    return three_end_pos, n


############
#  rt+3cl  #
############
def draw_comp_PAS_rtcl(a1, a2, b1, b2, label1, label2, window_size = 3):
    plt.figure(figsize=(4.5, 3))
    
    df1 = pd.DataFrame(a1/a2).rolling(window_size, center=True).mean()
    df1.iloc[0:256] = np.nan
    plt.plot(df1, alpha=.65,
             label=label1, 
             color='#377EB8', 
             linewidth=3, 
             zorder=10) ### origin 2.6
    
    df2 = pd.DataFrame(b1/b2).rolling(window_size, center=True).mean()
    df2.iloc[0:256] = np.nan
    plt.plot(df2, alpha=.65, 
             label=label2, linewidth=3, color='#E41A1C')
    
    # plt.ylim(-0.00025, 0.005)
    plt.ylim(0, 0.006)
    # plt.xlim(0)
    
    #ax[1].axvline(200, ls='--', color='#555555')
    plt.ylabel('Normalized density')

    xticks = np.array([0, 200, 400, 600, 800, 1000, 1200])

    plt.xticks(xticks, xticks-200)
    
    yticks = np.array([0, 0.00200, 0.004, 0.006])
    
    plt.yticks(yticks, yticks)
    
    plt.xlabel('Distance relative to PAS site (nt)')

    # plt.xlim(0, 1200)
    plt.xlim(200, 1200)
    plt.legend(frameon=False, prop={'size': 15})
    sns.despine(top=True, right=True)
    
    ax = plt.gca()  # 'get current axes'

    # 将左侧和底部的spines（轴线）移动到数据空间的0位置
    # ax.spines['bottom'].set_position(('data', -0.0002))
    # ax.spines['left'].set_position(('data', 180))
    
    colors_used = ['#EBDBD1']+ ['#E5E1EB'] + ['#CBE0F4']

    # add background of shadow
    
    # rect = mpatches.Rectangle((200, 0),
    #                             width=50, height=0.006, 
    #                             color=colors_used[1], 
    #                             alpha = 0.5,
    #                             #transform=ax.transAxes, 
    #                             clip_on=False)
    
    # plt.gca().add_patch(rect)
    
    return plt.gcf()

############
#    3cl   #
############
def draw_comp_PAS_rt(a1, a2, b1, b2, label1, label2, window_size = 3):
    plt.figure(figsize=(4.5, 3))
    # plt.figure(figsize=(4, 2))
    
    df1 = pd.DataFrame(a1/a2).rolling(window_size, center=True).mean()
    df1.iloc[0:256] = np.nan
    plt.plot(df1, alpha=.65,
             label=label1, color='#377EB8', linewidth=2.4, zorder=10)
    
    df2 = pd.DataFrame(b1/b2).rolling(window_size, center=True).mean()
    df2.iloc[0:256] = np.nan
    plt.plot(df2, alpha=.65, 
             label=label2, linewidth=2.4, color='#E41A1C')
    
    # plt.ylim(-0.00025, 0.005)
    plt.ylim(0, 0.006)
    plt.xlim(0)
    
    #ax[1].axvline(200, ls='--', color='#555555')
    plt.ylabel('Normalized density')

    xticks = np.array([0, 200, 400, 600, 800, 1000, 1200])

    plt.xticks(xticks, xticks-200)

    yticks = np.array([0, 0.00200, 0.004, 0.006])
    
    plt.yticks(yticks, yticks)

    plt.xlabel('Distance relative to PAS site (nt)')
    # plt.xlim(180, 1200)
    plt.xlim(200, 1200)
    plt.legend(frameon=False, prop={'size': 15})
    sns.despine(top=True, right=True)
    
    ax = plt.gca()  # 'get current axes'

    # 将左侧和底部的spines（轴线）移动到数据空间的0位置
    # ax.spines['bottom'].set_position(('data', -0.0002))
    # ax.spines['left'].set_position(('data', 180))
    
    colors_used = ['#EBDBD1']+ ['#E5E1EB'] + ['#CBE0F4']
    
    # add background of shadow
    # rect = mpatches.Rectangle((201, 0.00005),
    #                             width=51, height=0.006, 
    #                             color=colors_used[1], 
    #                             alpha = 0.5,
    #                             #transform=ax.transAxes, 
    #                             clip_on=False)
    # rect.set_zorder(11)
    # plt.gca().add_patch(rect)
    
    return plt.gcf()

############
#    3cl   #
############
def draw_comp_PAS_3cl(a1, a2, b1, b2, label1, label2, window_size = 3):
    plt.figure(figsize=(4.5, 3))
    # plt.figure(figsize=(4, 2))
    
    df1 = pd.DataFrame(a1/a2).rolling(window_size, center=True).mean()
    df1.iloc[0:270] = np.nan
    plt.plot(df1, alpha=.65,
             label=label1, color='#377EB8', linewidth=2.4, zorder=10)
    
    df2 = pd.DataFrame(b1/b2).rolling(window_size, center=True).mean()
    df2.iloc[0:270] = np.nan
    plt.plot(df2, alpha=.65, 
             label=label2, linewidth=2.4, color='#E41A1C')
    
    # plt.ylim(-0.00025, 0.005)
    plt.ylim(0, 0.006)
    plt.xlim(0)
    
     
    #ax[1].axvline(200, ls='--', color='#555555')
    plt.ylabel('Normalized density')

    xticks = np.array([0, 200, 400, 600, 800, 1000, 1200])

    plt.xticks(xticks, xticks-200)

    yticks = np.array([0, 0.00200, 0.004, 0.006])
    
    plt.yticks(yticks, yticks)
    
    plt.xlabel('Distance relative to PAS site (nt)')
    # plt.xlim(180, 1200)
    plt.xlim(200, 1200)
    plt.legend(frameon=False, prop={'size': 15})
    sns.despine(top=True, right=True)
    
    ax = plt.gca()  # 'get current axes'

    # 将左侧和底部的spines（轴线）移动到数据空间的0位置
    # ax.spines['bottom'].set_position(('data', -0.0002))
    # ax.spines['left'].set_position(('data', 180))
    
    colors_used = ['#EBDBD1']+ ['#E5E1EB'] + ['#CBE0F4']
    
    # add background of shadow
    # rect = mpatches.Rectangle((201, 0.00005),
    #                             width=200, height=0.006, 
    #                             color=colors_used[1], 
    #                             alpha = 0.5,
    #                             #transform=ax.transAxes, 
    #                             clip_on=False)
    # rect.set_zorder(11)
    # plt.gca().add_patch(rect)
    
    return plt.gcf()

def get_bin_cov(data: list, bins: int):
    data = np.array(data)
    data = np.nan_to_num(data)
    
    if bins == 1:
        return data
    if bins < 1:
        raise ValueError('bins must be greater than 1')

    results = []
    for i in range(0, len(data), bins):
        bin_data = data[i:i + bins]
        mean_bin_data = np.nanmean(bin_data)
        results.append(mean_bin_data)

    return results


################
# For bw file  #
################

# site-point
def bw_site_cov(
    infile: str, 
    chrom: str, site: int, strand: str,
    before: int = 1000, after : int = 1000,
    bins: int = 100,
    chrom_prefix: str = '',
    normalized: str = 'density',
    exclude_chr = None
    ):
    '''
    Args:
        infile: path to bigWig file
        chrom: chromosome name
        site: target site
        strand: '+' or '-'
        before: distance upstream of the site1 selected
        after: distance downstream of the site2 selected
        regionbody: distance in bases to which all regions will be fit
        bins: length in bases, of the non-overlapping bins for averaging the score over the regions length
        chrom_prefix: prefix of the chromosome name, eg. "chr"
        exclude_chr: chromosomes to be excluded
    
    Return:
        cov: the coverage value of given regions
    '''
    site = int(site)
    chrom = str(chrom)

    if exclude_chr is not None and chrom in exclude_chr:
        return
    chrom = chrom_prefix + chrom
        
    bwfile = pyBigWig.open(infile)
    if strand == '+':
        start = site - before
        end = site + after
    else:
        start = site - after
        end = site + before
    if start < 0 or end > bwfile.chroms()[chrom]:
        # remove Invalid interval
        return

    if strand == '+':
        values = bwfile.values(
            chrom, 
            start,
            end)
        cov = get_bin_cov(values, bins)

    elif strand == '-':
        values = bwfile.values(
            chrom, 
            start,
            end)[::-1]
        cov = get_bin_cov(values, bins)
    else:
        return ValueError('strand must be "+" or "-"')
    
    cov = np.nan_to_num(cov)
    if sum(cov) > 0:
        if normalized == 'density':
            cov = cov / sum(cov)  # density
        return cov

    
def bw_reference_point(
    infile: str, 
    site_info: list,
    before: int = 1000, after : int = 1000,
    bins: int = 100,
    chrom_prefix: str = '',
    normalized: str = 'density',
    exclude_chr = None,
    threads=64):
    '''
    Reference-point refers to a position within a BED region (e.g., the starting point). In this mode, only those genomicpositions before (upstream) and/or after (downstream) of the reference point will be used.

    Args:
        infile: path to bigWig file
        site_info: [(chrom, site1, site2, strand), ...]
        before: distance upstream of the site1 selected
        after: distance downstream of the site2 selected
        bins: length in bases, of the non-overlapping bins for averaging the score over the regions length
        chrom_prefix: prefix of the chromosome name, eg. "chr"
        exclude_chr: chromosomes to be excluded
    
    Return:
        cov: the coverage value of givin regions
    '''
    chrom = site_info[:, 0]
    site = site_info[:, 1]
    strand = site_info[:, 2]
    with ProcessPoolExecutor(max_workers=threads) as e:
        chunksize = int(len(site_info) / threads)
        results = e.map(
            bw_site_cov,
            repeat(infile),
            chrom,
            site,
            strand,
            repeat(before),
            repeat(after),
            repeat(bins),
            repeat(chrom_prefix),
            repeat(normalized),
            repeat(exclude_chr),
            chunksize=chunksize)

    cov = []
    for cov_ in results:
        if cov_ is not None:
            cov.append(cov_)

    cov = np.nanmean(cov, axis=0)
    return cov


# scale region
def bw_scale_cov(
    infile: str, 
    chrom: str, site1: int, site2: int, strand: str,
    before: int = 1000, after : int = 1000, regionbody : int = 1000, 
    bins: int = 100,
    split: bool = False,
    chrom_prefix: str = '',
    normalized: str = 'density',
    exclude_chr = None
    ):
    '''
    Args:
        infile: path to bigWig file
        chrom: chromosome name
        site1: 5' site
        site2: 3' site
        strand: '+' or '-'
        before: distance upstream of the site1 selected
        after: distance downstream of the site2 selected
        regionbody: distance in bases to which all regions will be fit
        bins: length in bases, of the non-overlapping bins for averaging the score over the regions length
        chrom_prefix: prefix of the chromosome name, eg. "chr"
        normalized: normalization method, 'density' or 'count'
        exclude_chr: chromosomes to be excluded
    
    Return:
        cov: the coverage value of givin regions
    '''
    site1 = int(site1)
    site2 = int(site2)
    chrom = str(chrom)

    if exclude_chr is not None and chrom in exclude_chr:
        return

    chrom = chrom_prefix + chrom

    bwfile = pyBigWig.open(infile)
    if strand == '+':
        start = site1 - before
        end = site2 + after
    else:
        start = site1 - after
        end = site2 + before
    if start < 0 or end > bwfile.chroms()[chrom]:
        # remove Invalid interval
        return

    if split:
        # in this mode, regionbody is ignored
        if strand == '+':
            cov_5 = bwfile.values(chrom, site1 - before, site1 + after)
            cov_5 = get_bin_cov(cov_5, bins)
            cov_3 = bwfile.values(chrom, site2 - before, site2 + after)
            cov_3 = get_bin_cov(cov_3, bins)

        elif strand == '-':
            cov_5 = bwfile.values(chrom, site2 - after, site2 + before)[::-1]
            cov_5 = get_bin_cov(cov_5, bins)
            cov_3 = bwfile.values(chrom, site1 - after, site1 + before)[::-1]
            cov_3 = get_bin_cov(cov_3, bins)
        
        cov_5 = np.nan_to_num(cov_5)
        cov_3 = np.nan_to_num(cov_3)

        sum_cov = sum(cov_5)+sum(cov_3)
        if sum_cov > 0:
            # density
            cov_5 = cov_5 / sum_cov
            cov_3 = cov_3 / sum_cov
            return cov_5, cov_3

    else:
        if strand == '+':
            start = site1 - before
            end = site2 + after
        else:
            start = site1 - after
            end = site2 + before
        
        if start < 0 or end > bwfile.chroms()[chrom]:
            return

        if strand == '+':
            cov_5 = bwfile.values(chrom, start, site1)
            cov_5 = get_bin_cov(cov_5, bins)
            cov_3 = bwfile.values(chrom, site2, end)
            cov_3 = get_bin_cov(cov_3, bins)
            # gene_body_region
            cov_gb = bwfile.values(chrom, site1, site2)
            cov_gb = scipy.ndimage.zoom(
                cov_gb,
                regionbody / len(cov_gb),
                order=0,
                mode='nearest')
            cov_gb = get_bin_cov(cov_gb, bins)

        elif strand == '-':
            cov_5 = bwfile.values(chrom, site2, end)[::-1]
            cov_5 = get_bin_cov(cov_5, bins)
            cov_3 = bwfile.values(chrom, start, site1)[::-1]
            cov_3 = get_bin_cov(cov_3, bins)
            # gene_body_region
            cov_gb = bwfile.values(chrom, site1, site2)[::-1]
            cov_gb = scipy.ndimage.zoom(
                cov_gb,
                regionbody / len(cov_gb),
                order=0,
                mode='nearest')
            cov_gb = get_bin_cov(cov_gb, bins)
        else:
            raise ValueError('strand must be "-" or "+"')

        cov = np.concatenate([cov_5, cov_gb, cov_3])
        cov = np.nan_to_num(cov)

        if sum(cov) > 0:
            if normalized == 'density':
                cov = cov / sum(cov)  # density
            return cov


def bw_scale_regions(
    infile: str, 
    site_info: list,
    before: int = 1000, after : int = 1000, regionbody : int = 1000, 
    bins: int = 100,
    split: bool = False,
    chrom_prefix: str = '',
    normalized: str = 'density',
    exclude_chr = None,
    threads=64):
    '''
    In the scale-regions mode, all regions in the BED file are stretched or shrunken to the length (in bases) indicated by the user.

    Args:
        infile: path to bigWig file
        site_info: [(chrom, site1, site2, strand), ...]
        before: distance upstream of the site1 selected
        after: distance downstream of the site2 selected
        regionbody: distance in bases to which all regions will be fit
        bins: length in bases, of the non-overlapping bins for averaging the score over the regions length
        chrom_prefix: prefix of the chromosome name, eg. "chr"
        normalized: normalization method, 'density' or 'count'
        exclude_chr: chromosomes to be excluded
    
    Return:
        cov: the coverage value of givin regions
    '''
    site_info = np.array(site_info)
    chrom = site_info[:, 0]
    site1 = site_info[:, 1]
    site2 = site_info[:, 2]
    strand = site_info[:, 3]

    with ProcessPoolExecutor(max_workers=threads) as e:
        chunksize = int(len(site_info) / threads)
        results = e.map(
            bw_scale_cov,
            repeat(infile),
            chrom,
            site1,
            site2,
            strand,
            repeat(before),
            repeat(after),
            repeat(regionbody),
            repeat(bins),
            repeat(split),
            repeat(chrom_prefix),
            repeat(normalized),
            repeat(exclude_chr),
            chunksize=chunksize)

    if split:
        cov_5, cov_3 = [], []
        for res in results:
            if res is not None:
                cov_5_, cov_3_ = res
                cov_5.append(cov_5_)
                cov_3.append(cov_3_)

        cov_5 = np.nanmean(cov_5, axis=0)
        cov_3 = np.nanmean(cov_3, axis=0)
        return cov_5, cov_3

    else:
        cov = []
        for cov_ in results:
            if cov_ is not None:
                cov.append(cov_)

        cov = np.nanmean(cov, axis=0)
        return cov


################
# For bam file #
################

# site-point
def bam_site_cov(
    infile: str, 
    chrom: str, site: int, strand: str, gene_id: str,
    before: int = 1000, after : int = 1000,
    bins: int = 100,
    min_counts: int = 1,
    chrom_prefix: str = '',
    exclude_chr: set = None,
    return_raw: bool = False,
):
    """
    BAM file for tagged FLEP-seq data
    Ignore splicing junction

    Args:
        infile: path to bigWig file
        chrom: chromosome name
        site: target site
        strand: '+' or '-'
        gene_id: gene id
        before: distance upstream of the site1 selected
        after: distance downstream of the site2 selected
        bins: length in bases, of the non-overlapping bins for averaging the score over the regions length
        min_counts: minimum number of the reads within given region
        chrom_prefix: prefix of the chromosome name, eg. "chr"
        exclude_chr: chromosomes to be excluded
    
    Return:
        cov: the coverage value of given region
    """
    site = int(site)
    chrom = str(chrom)
    strand_is_reverse = False if strand == '+' else True

    if exclude_chr is not None and chrom in exclude_chr:
        return
    chrom = chrom_prefix + chrom
    
    n = 0

    cov = np.zeros(before+after)

    if strand == '+':
        start = site-before
        end = site+after
    else:
        start = site-after
        end = site+before
    

    with pysam.AlignmentFile(infile, 'rb') as inbam:
        if start < 0 or end > inbam.get_reference_length(chrom):
            return
            
        for read in inbam.fetch(chrom, start, end):
            # 判断是否跟基因是同个方向，针对于链特异文库
            if read.is_reverse != strand_is_reverse:
                continue
            
            read_gene_id = read.get_tag('gi')
            if read_gene_id not in {gene_id, 'None'}:
                continue
                
            if strand == '+':
                read_five_end = read.reference_start
                read_three_end = read.reference_end
                cov_start = read_five_end-start if read_five_end-start >= 0 else 0
                cov_end = read_three_end-start if read_three_end-start <= before+after else end-start
            else:
                read_five_end = read.reference_end
                read_three_end = read.reference_start
                cov_start = end-read_five_end if end-read_five_end >= 0 else 0
                cov_end = end-read_three_end if end-read_three_end <= before+after else end-start

            cov[cov_start: cov_end] += 1            
            n += 1
    
    if return_raw:
        return cov, n

    if n > min_counts:
        if bins > 1:
            cov = get_bin_cov(cov, bins)
        cov = cov / sum(cov)
        return cov


def bam_reference_point(
    infile: str, 
    site_info: list,
    before: int = 1000, after : int = 1000,
    bins: int = 100,
    min_counts: int = 1,
    chrom_prefix: str = '',
    exclude_chr = None,
    return_raw: bool = False,
    threads=64):
    '''
    Reference-point refers to a position within a BED region (e.g., the starting point). In this mode, only those genomicpositions before (upstream) and/or after (downstream) of the reference point will be used.

    Args:
        infile: path to bigWig file
        site_info: [(chrom, site, strand, gene_id), ...]
        before: distance upstream of the site1 selected
        after: distance downstream of the site2 selected
        bins: length in bases, of the non-overlapping bins for averaging the score over the regions length
        min_counts: minimum number of the reads
        chrom_prefix: prefix of the chromosome name, eg. "chr"
        exclude_chr: chromosomes to be excluded
    
    Return:
        cov: the coverage value of givin regions
    '''
    chrom = site_info[:, 0]
    site = site_info[:, 1]
    strand = site_info[:, 2]
    gene_id = site_info[:, 3]

    with ProcessPoolExecutor(max_workers=threads) as e:
        chunksize = int(len(site_info)/threads)
        results = e.map(
            bam_site_cov, 
            repeat(infile),
            chrom,
            site,
            strand,
            gene_id,
            repeat(before),
            repeat(after),
            repeat(bins),
            repeat(min_counts),
            repeat(chrom_prefix),
            repeat(exclude_chr),
            repeat(return_raw),
            chunksize=chunksize)

    if return_raw:
        cov = []
        for gene, res in zip(gene_id, results):
            if res is not None:
                res, n = res
                cov.append((gene, res, n))
        res = pd.DataFrame([x[1] for x in cov], index=[x[0] for x in cov], columns=range(-before, after))
        res['N'] = [x[2] for x in cov]
        return res
    else:
        cov = []
        n = 0
        for res in results:
            if res is not None:
                cov.append(res)
        
        cov = np.nanmean(cov, axis=0)
    return cov


# scale region
def bam_scale_cov(
    infile: str, 
    chrom: str, site1: int, site2: int, strand: str, gene_id: str,
    before: int = 1000, after : int = 1000, regionbody : int = 1000, 
    bins: int = 100,
    split: bool = False,
    min_counts: int = 1,
    chrom_prefix: str = '',
    exclude_chr = None
    ):
    '''
    Args:
        infile: path to BAM file
        chrom: chromosome name
        site1: 5' site
        site2: 3' site
        strand: '+' or '-'
        before: distance upstream of the site1 selected
        after: distance downstream of the site2 selected
        regionbody: distance in bases to which all regions will be fit
        bins: length in bases, of the non-overlapping bins for averaging the score over the regions length
        split: split mode
        min_count: minimum number of reads in the region
        chrom_prefix: prefix of the chromosome name, eg. "chr"
        exclude_chr: chromosomes to be excluded
    
    Return:
        cov: the coverage value of givin regions
    '''
    site1 = int(site1)
    site2 = int(site2)
    region_len = site2-site1
    chrom = str(chrom)
    strand_is_reverse = False if strand == '+' else True

    if exclude_chr is not None and chrom in exclude_chr:
        return

    chrom = chrom_prefix + chrom

    if split:
        cov1 = bam_site_cov(infile, chrom, site1, strand, gene_id, before=before, after=after, bins=bins, return_raw=True)
        cov2 = bam_site_cov(infile, chrom, site2, strand, gene_id, before=before, after=after, bins=bins, return_raw=True)

        if cov1 is None or cov2 is None:
            return
            
        sum_cov = sum(cov1) + sum(cov2)
        if sum_cov > 0:
            cov1 = cov1 / sum_cov
            cov2 = cov2 / sum_cov
            if strand == '+':
                return cov1, cov2
            elif strand == '-':
                return cov2, cov1
            
    else:
        if strand == '+':
            start = site1-before
            end = site2+after
        else:
            start = site1-after
            end = site2+before

        cov = np.zeros(before+region_len+after)

        with pysam.AlignmentFile(infile, 'rb') as inbam:
            if start < 0 or end > inbam.get_reference_length(chrom):
                return
            n = 0
            for read in inbam.fetch(chrom, start, end):
                # 判断是否跟基因是同个方向，针对于链特异文库
                if read.is_reverse != strand_is_reverse:
                    continue
                    
                read_gene_id = read.get_tag('gi')
                if read_gene_id not in {gene_id, 'None'}:
                    continue
                    
                if strand == '+':
                    read_five_end = read.reference_start
                    read_three_end = read.reference_end
                    cov_start = read_five_end-start if read_five_end-start >= 0 else 0
                    cov_end = read_three_end-start if read_three_end-start <= before+after+region_len else end-start
                else:
                    read_five_end = read.reference_end
                    read_three_end = read.reference_start
                    cov_start = end-read_five_end if end-read_five_end >= 0 else 0
                    cov_end = end-read_three_end if end-read_three_end <= before+after+region_len else end-start

                cov[cov_start: cov_end] += 1            
                n += 1
            
            if n > min_counts:
                cov_gb = cov[before: before+region_len]
                cov_gb = scipy.ndimage.zoom(
                    cov_gb,
                    regionbody / len(cov_gb),
                    order=0,
                    mode='nearest')
                cov_5 = cov[ : before]
                cov_3 = cov[before+region_len : ]
                cov = np.concatenate([cov_5, cov_gb, cov_3])

                if bins > 1:
                    cov = get_bin_cov(cov, bins)
                cov = cov / sum(cov)
                
                return cov


def bam_scale_region(
    infile: str, 
    site_info: list,
    before: int = 1000, after: int = 1000, regionbody: int = 1000,
    bins: int = 100,
    split: bool = False,
    min_counts: int = 1,
    chrom_prefix: str = '',
    exclude_chr = None,
    threads=64):
    '''
    Reference-point refers to a position within a BED region (e.g., the starting point). In this mode, only those genomicpositions before (upstream) and/or after (downstream) of the reference point will be used.

    Args:
        infile: path to BAM file
        site_info: [(chrom, site1, site2, strand, gene_id), ...]
        before: distance upstream of the site1 selected
        after: distance downstream of the site2 selected
        regionbody: distance in bases to which all regions will be fit
        bins: length in bases, of the non-overlapping bins for averaging the score over the regions length
        split: split mode
        min_counts: minimum number of the reads
        chrom_prefix: prefix of the chromosome name, eg. "chr"
        exclude_chr: chromosomes to be excluded
    
    Return:
        cov: the coverage value of givin regions
    '''
    chrom = site_info[:, 0]
    site1 = site_info[:, 1]
    site2 = site_info[:, 2]
    strand = site_info[:, 3]
    gene_id = site_info[:, 4]

    with ProcessPoolExecutor(max_workers=threads) as e:
        chunksize = int(len(site_info)/threads)
        results = e.map(
            bam_scale_cov, 
            repeat(infile),
            chrom,
            site1,
            site2,
            strand,
            gene_id,
            repeat(before),
            repeat(after),
            repeat(regionbody),
            repeat(bins),
            repeat(split),
            repeat(min_counts),
            repeat(chrom_prefix),
            repeat(exclude_chr),
            chunksize=chunksize)
    

    if split:
        cov_5, cov_3 = [], []
        for res in results:
            if res is not None:
                cov_5_, cov_3_ = res
                cov_5.append(cov_5_)
                cov_3.append(cov_3_)

        cov_5 = np.nanmean(cov_5, axis=0)
        cov_3 = np.nanmean(cov_3, axis=0)
        
        return cov_5, cov_3

    else:
        cov = []
        for res in results:
            if res is not None:
                cov.append(res)
        
        cov = np.nanmean(cov, axis=0)
        return cov


def get_bam_total_readcounts(infile: str):
    """
    This function takes a bam file and returns the total number of reads in the file.
    
    Args:
      infile (str): the input bam file
    """

    return eval('+'.join([line.split('\t')[2] for line in pysam.idxstats(infile).rstrip().split('\n')]))



def plot(ax, cov, n, before=2000, after=2000, target_site=0, label=None, ylabel=None):
    """
    画metaplot
    """
    ax.plot(cov/n, label=label)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    ax.set_xticks([0, before, before+after])
    ax.set_xticklabels([f'-{before//1000} Kb', target_site, f'{after//1000} Kb'], rotation=90)
    
    ax.axvline(before, ls='--', color='#555555')
    if label is not None:
        ax.legend(frameon=False)


def set_ax(
    ax, 
    bins,
    b1: int = None, a1: int = None, 
    b2: int = None, a2: int = None, 
    site1: str = 0, site2: str = 0,
    ylabel=None
    ):
    """
    This function takes in a matplotlib axis object, a list of bins, and two sets of bin numbers and
    site names, and returns a histogram of the bins with the two sets of bins highlighted
    
    Args:
      ax: the axis to plot on
      bins: the bins for the histogram
      b1 (int): int = None, a1: int = None,
      a1 (int): int = None, b1: int = None, a2: int = None, b2: int = None,
      b2 (int): int = None, a2: int = None,
      a2 (int): int = None,
      site1 (str): str = 0, site2: str = 0,. Defaults to 0
      site2 (str): str = 0,. Defaults to 0
      ylabel: str
    """
    if type(ax) is not np.ndarray:
        ax = [ax]
    
    if b1 is not None and a1 is not None:
        ax[0].spines['right'].set_visible(False)
        ax[0].spines['top'].set_visible(False)
        ax[0].set_ylabel(ylabel)
        ax[0].set_xticks([0, b1//bins, (a1+b1)//bins])
        ax[0].set_xticklabels([f'-{b1//1000} kb', site1, f'{a1//1000} kb'], rotation=90)
        ax[0].axvline(b1//bins, ls='--', color='#555555')
    
    if b2 is not None and a2 is not None:
        ax[1].spines['right'].set_visible(False)
        ax[1].spines['left'].set_visible(False)
        ax[1].spines['top'].set_visible(False)
        ax[1].yaxis.set_ticks_position('none')
        ax[1].set_xticks([0, b2//bins, (a2+b2)//bins])
        ax[1].set_xticklabels([f'-{b2//1000} kb', site2, f'{a2//1000} kb'], rotation=90)
        ax[1].axvline(b2//bins, ls='--', color='#555555')