
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author       : windz
Date         : 2022-03-20 18:03:56
LastEditTime : 2022-04-14 19:01:13
LastEditors  : windz
FilePath: flepseq.py
'''

'''
Modified_Author   : Juzuo
Date              : 2024-03-12
LastEditTime      : 2024-03
FilePath          : flepseq.py

Modified part: Function about classifying the readthrough/ 3′cleavage products/ 5′cleavage products.
The method is identical with the GB's description bellow:
Non-poly(A) reads with 5′ ends located at gene body ( > 50 nt upstream of poly(A) site) and 
3′ end located more than 50 nt downstream of the poly(A) site are considered as "readthrough transcripts.
Non-poly(A) reads with 5′ends located at gene body as well as 3′ ends within 50 nt upstream and 
downstream of poly(A) sites are considered as "5′ cleavage products.
Non-poly(A) reads with 5′ end within 200 nt downstream of poly(A) sites are regarded as "3′ cleavage products.
'''

import pysam
import numpy as np
import pandas as pd
import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from collections import Counter


def get_read_polya_len(infile: str, min_len: int = 15):
    '''
    Get read polyA length from BAM file
    Args:
        infile: path to BAM file
        min_len: minimum of the polyA tail length
    
    Return:
        polya_len_list: polya length of all reads in np.array
    '''
    polya_len_list = []
    with pysam.AlignmentFile(infile, 'rb') as bam:
        for read in bam:
            # require unique mapped read
            if read.is_unmapped or read.is_supplementary or read.is_secondary:
                continue
            polya_len = read.get_tag('pa')
            if polya_len >= min_len:
                polya_len_list.append(polya_len)
        
        polya_len_list = np.array(polya_len_list)

    return polya_len_list


def get_read_polya_len_with_splice_info(
    infile: str, 
    min_len: int = 15, 
    exclude_chr: tuple = ('Pt', 'Mt'),
):
    '''
    Get polyA length of spliced and unspliced reads 
    Args:
        infile: path to BAM file
        min_len: minimum tail length of the polyA
        exclude_chr: the chromosomes to be exclude
    
    Return:
        splice_res, incompletely_splice_res: gene polya length info with splice info
    '''
    from collections import defaultdict


    exclude_chr = set(exclude_chr)
    polya_len_incompletely_spliced = []
    polya_len_spliced = []
    
    with pysam.AlignmentFile(infile, 'rb') as bam:
        for read in bam:
            if read.reference_name in exclude_chr:
                continue
                
            if read.is_unmapped or read.is_supplementary or read.is_secondary:
                continue
                
            polya_len = read.get_tag('pa')
            gene_id = read.get_tag('gi')
            if gene_id == 'None':
                continue
            if polya_len >= min_len:
                if read.get_tag('rn') == 0:
                    polya_len_spliced.append(polya_len)
                else:
                    polya_len_incompletely_spliced.append(polya_len)
    
    polya_len_spliced = np.array(polya_len_spliced)
    polya_len_incompletely_spliced = np.array(polya_len_incompletely_spliced)
            
    return polya_len_spliced, polya_len_incompletely_spliced


def get_gene_polya_len(infile, min_len=15, min_count=10, exclude_chr=('Pt', 'Mt')):
    '''
    Get gene median polyA length
    Args:
        infile: path to BAM file
        min_len: minimum of the polyA tail length
        min_count: minimum read counts support the gene polyA length
        exclude_chr: the chromosomes to be exclude
    
    Return:
        gene_polya_len: gene polya length in np.array
                        (('gene_id', length), ...)
        gene_polya_len_raw: raw polyA length of the gene, dict
    '''
    from collections import defaultdict

    exclude_chr = set(exclude_chr)
    gene_polya_len_raw = defaultdict(lambda : []) 
    with pysam.AlignmentFile(infile, 'rb') as bam:
        for read in bam:
            if read.reference_name in exclude_chr:
                continue
                
            if read.is_unmapped or read.is_supplementary or read.is_secondary:
                continue
                
            polya_len = read.get_tag('pa')
            gene_id = read.get_tag('gi')
            if gene_id == 'None':
                continue
                
            if polya_len >= min_len:
                gene_polya_len_raw[gene_id].append(polya_len)
    
    gene_polya_len = []
    for gene_id in gene_polya_len_raw:
        if len(gene_polya_len_raw[gene_id]) >= min_count:
            gene_polya_len.append(np.median(gene_polya_len_raw[gene_id]))
    
    gene_polya_len = np.array(gene_polya_len)
            
    return gene_polya_len, gene_polya_len_raw


def get_gene_polya_len_with_splice_info(
    infile: str, 
    min_len: 15,
    min_count: int = 15, 
    exclude_chr: tuple = ('Pt', 'Mt'),
):
    '''
    Get polyA length between spliced and unspliced read of the gene 
    Args:
        infile: path to BAM file
        min_len: minimum read length of the polyA length
        min_count: minimum read count to support the median polyA length
        exclude_chr: the chromosomes to be exclude
    
    Return:
        splice_res, incompletely_splice_res: gene polya length info with splice info
    '''
    from collections import defaultdict


    exclude_chr = set(exclude_chr)
    gene_polya_len_unspliced = defaultdict(lambda : [])
    gene_polya_len_spliced = defaultdict(lambda : []) 
    
    with pysam.AlignmentFile(infile, 'rb') as bam:
        for read in bam:
            if read.reference_name in exclude_chr:
                continue
                
            if read.is_unmapped or read.is_supplementary or read.is_secondary:
                continue
                
            polya_len = read.get_tag('pa')
            gene_id = read.get_tag('gi')
            if gene_id == 'None':
                continue
            if polya_len >= min_len:
                if read.get_tag('rn') == 0:
                    gene_polya_len_spliced[gene_id].append(polya_len)
                else:
                    gene_polya_len_unspliced[gene_id].append(polya_len)
                    
    spliced_res, incompletely_spliced_res = {}, {}
    for k, v in gene_polya_len_spliced.items():
        if len(v) >= min_count:
            spliced_res[k] = np.median(v)
            
    for k, v in gene_polya_len_unspliced.items():
        if len(v) >= min_count:
            incompletely_spliced_res[k] = np.median(v)
            
    return spliced_res, incompletely_spliced_res

    
def get_idxstats(inbam: str):
    '''
    Get total read number from BAM file
    Args:
        inbam: path to BAM file
    
    Return:
        total_count: total count in BAM file
        mapped_count: mapped count in BAM file
        unmapped_count: unmapped count in BAM file
        count_per_chr: mapped count per chromosomes
    '''
    count_per_chr = {}
    total_count, mapped_count, unmapped_count = 0, 0, 0
    for l in pysam.idxstats(inbam).split('\n')[:-1]:
        l = l.rstrip().split('\t')
        total_count += eval('+'.join(l[2:]))
        if l[0] != '*':
            count_per_chr[l[0]] = int(l[2])
            mapped_count += int(l[2])
        else:
            unmapped_count += int(l[3])
    
    return total_count, mapped_count, unmapped_count, count_per_chr


def get_gene_counts(
    bed_intersect: str = None, 
    inbam : str = None, repr_gene : str = None, 
    exclude_chr: set = {'Pt', 'Mt'}
    ) -> dict:
    '''
    计算有reads覆盖的基因个数
    Args:
        bed_intersect: bamfile intersect with gene model
            # bedtools intersect -abam {input} -b {params.repr_gene} -s -wo -split -bed > {output}

        **如果没有提供bed_intersect结果, 可以调用pybedtools进行计算, 二选一**
        bam_file: path to BAM file
        repr_gene: gene model in BED6 format
        exclude_chr: chromosomes to be excluded
    
    Return:
        gene_counts: read count per gene, Counter
    '''
    from collections import Counter
    from utils.get_overlap_genes import exclude_genes

    
    NAMES=[
    'Chromosome', 'Start', 'End', 'read_id', 'Score', 'Strand', 
    'ThickStart', 'ThickEnd', 'ItemRGB', 'BlockCount', 'BlockSizes', 'BlockStarts', 
    'geneChromosome', 'geneStart', 'geneEnd', 'gene_id', 'geneScore', 'geneStrand', 
    'cov'
    ]

    USECOLS = [
        'Chromosome', 'Start', 'End', 'read_id', 'Strand',
        'geneStart', 'geneEnd', 'gene_id', 'geneStrand', 'cov',
    ]
    
    if bed_intersect is not None:
        df = pd.read_csv(
            bed_intersect, 
            sep='\t', 
            names=NAMES,
            usecols=USECOLS,
            header=None,
            dtype={"Chromosome": str}
            )
    elif inbam is not None and repr_gene is not None:
        from pybedtools import BedTool

        bam = BedTool(inbam)
        if not bam._isbam:
            raise ValueError('inbam is not BAM file')
        res = bam.intersect(repr_gene, s=True, wo=True, split=True, bed=True)
        df = res.to_dataframe(
            disable_auto_names=True, 
            names=NAMES, 
            usecols=USECOLS, 
            dtype={"Chromosome": str}
            )
    else:
        raise ValueError('Please no input found')

    gene_counts = Counter()
    for item in df.itertuples():
        if item.Chromosome in exclude_chr:
            continue
        if item.Strand != item.geneStrand:
            continue
        if item.gene_id in exclude_genes:
            continue

        if item.Strand == '+' and item.Start >= item.geneStart-100:
            gene_counts[item.gene_id] += 1
        elif item.Strand == '-' and item.End <= item.geneEnd+100:
            gene_counts[item.gene_id] += 1
    
    return gene_counts


def termination_window(
    infile: str,
    chrom: str,
    pas: int,
    strand: str,
    gene_id: str,
    method: str = 'median',
    min_counts: int = 15,
    exclude_chr: set = None,
    verbose: bool = False
    ):
    '''> For each gene, find the median distance between the 3' end of the read and the PAS site
    
    Parameters
    ----------
    infile : str
        the bam file
    chrom : str
        chromosome name
    pas : int
        the polyA site
    strand : str
        '+' or '-'
    gene_id : str
        the gene id of the gene you want to find the termination site for
    method : {'median', 'raw'}
        default: 'median'
        method to return the termination window
         - 'median': the median distance between the 3' end of the read and the PAS site
         - 'raw': the raw data of the distance list
    min_counts : int, optional
        minimum number of reads to call a termination site
    exclude_chr : set
        set = None
    verbose : bool, optional
        print out the read names and distances
    
    Returns
    -------
        raw_data or median readthrough length
    '''
    pas = int(pas)
    chrom = str(chrom)

    if exclude_chr is not None and chrom in exclude_chr:
        exclude_chr = set(exclude_chr)
        if chrom in exclude_chr:
            return

    strand_is_reverse = False if strand == '+' else True  # gene strand

    distance_list = []
    with pysam.AlignmentFile(infile, 'rb') as inbam:

        if strand_is_reverse:
            start = pas - 200 if pas > 200 else 1
            end = pas
        else:
            start = pas
            end = pas + 200 if pas + 200 < inbam.get_reference_length(chrom) else inbam.get_reference_length(chrom)

        for read in inbam.fetch(chrom, start, end):
            if read.is_reverse != strand_is_reverse:
                continue
            
            if read.is_unmapped or read.is_supplementary or read.is_secondary:
                continue
            
            read_gene_id = read.get_tag('gi')
            polya_len = read.get_tag('pa')
            if polya_len >= 15 or read_gene_id not in {gene_id, 'None'}:
                continue
            
            if strand_is_reverse:
                five_end = read.reference_end*-1
                three_end = read.reference_start*-1
                pa_site = pas*-1
            else:
                five_end = read.reference_start
                three_end = read.reference_end
                pa_site = pas
            
            distance3 = abs(three_end - pa_site)
            distance5 = abs(five_end - pa_site) 
            distance = three_end - pa_site
            if distance3 >= 50: ###and distance5 >= 50:
                distance_list.append(distance)
                if verbose:
                    print(f'{read.query_name} {distance}', file=sys.stderr)
    
    if len(distance_list) < min_counts:
        return
    
    if method == 'median':
        read_through_len = np.median(distance_list)
    elif method == 'raw':
        return gene_id, distance_list
    else:
        raise NameError(f'{method} not found')
    
    if strand_is_reverse:
        tts = pas - read_through_len
    else:
        tts = pas + read_through_len
    
    return read_through_len, gene_id, chrom, pas, tts, strand


def termination_window_turbo(
    infile: str,
    polya_info: list,
    method: str = 'raw',
    min_counts: int = 15,
    exclude_chr: set = None,
    threads: int = 10,
    ) -> pd.DataFrame:
    '''
    Get readthrough raw data from the BAM file
    
    Parameters
    ----------
    infile : str
        str = path to the BAM file
    polya_info : list
        list of tuples, each tuple is (chromosome, start, end, strand)
    method : str, optional
        'raw' or 'median'
    min_counts : int, optional
        minimum number of counts to consider a gene
    exclude_chr : set
        set = None,
    threads : int, optional
        number of threads to use
    
    '''
    polya_info = np.array(polya_info)  # the polyA site information
    chrom = polya_info[:, 0]
    pas = polya_info[:, 1]
    strand = polya_info[:, 2]
    gene_id = polya_info[:, 3]

    with ProcessPoolExecutor(max_workers=threads) as e:
        chunksize = int(len(polya_info)/threads)
        results = e.map(
            termination_window, 
            repeat(infile),
            chrom,
            pas,
            strand,
            gene_id,
            repeat(method),
            repeat(min_counts),
            repeat(exclude_chr),
            chunksize = chunksize
            )

    readthrough_data = []
    for res in results:
        if res is not None:
            readthrough_data.append(res)

    if method == 'raw':
        readthrough_data = pd.DataFrame(readthrough_data, columns=['gene_id', 'data'])
    else:
         readthrough_data = pd.DataFrame(readthrough_data, columns=['tw_size', 'gene_id', 'chrom', 'pas', 'tts', 'strand'])
    
    return readthrough_data



def count_rna_intermediates(
    infile: str,
    chrom: str,
    pas: int,
    strand: str,
    gene_id: str,
    min_counts: int = 15,
    exclude_chr: set = None,
    ):
    '''Counts the number of RNA intermediates in a given region of a genome
    
    Parameters
    ----------
    infile : str
        the name of the BAM file
    chrom : str
        chromosome id
    pas : int
        the position of the polya site
    strand : str
        "+" or "-"
    gene_id : str
        the gene id of the gene you want to look at
    min_counts : int, optional
        minimum number of reads required
    exclude_chr : set
        set = None,
    
    '''
    pas = int(pas)
    chrom = str(chrom)

    if exclude_chr is not None and chrom in exclude_chr:
        exclude_chr = set(exclude_chr)
        if chrom in exclude_chr:
            return

    strand_is_reverse = False if strand == '+' else True  # gene strand

    rna_intermediates = Counter()
    n = 0
    with pysam.AlignmentFile(infile, 'rb') as inbam:
        if strand_is_reverse:
            start = pas - 200 if pas > 200 else 1
            end = pas + 200
        else:
            start = pas - 200
            end = pas + 200 if pas + 200 < inbam.get_reference_length(chrom) else inbam.get_reference_length(chrom)

        for read in inbam.fetch(chrom, start, end):
            if read.is_reverse != strand_is_reverse:
                continue
            
            if read.is_unmapped or read.is_supplementary or read.is_secondary:
                continue
            
            read_gene_id = read.get_tag('gi')
            polya_len = read.get_tag('pa')
            if polya_len >= 15 or read_gene_id not in {gene_id, 'None'}:
                continue
            
            if strand_is_reverse:
                five_end = read.reference_end*-1
                three_end = read.reference_start*-1
                pa_site = pas*-1
            else:
                five_end = read.reference_start
                three_end = read.reference_end
                pa_site = pas
            
            distance = three_end - pa_site
            distance3 = three_end - pa_site
            distance5 = five_end - pa_site
            if five_end >= pa_site:
                rna_intermediates['3_cleavage'] += 1
                n += 1
            elif (distance3 > 50) and (distance5 < -50):
                rna_intermediates['readthrough'] += 1
                n += 1
            elif -50 <= distance3 <= 50:
                if distance5 <= -50:
                    rna_intermediates['5_cleavage'] += 1
    
    if n >= min_counts:
        return gene_id, rna_intermediates['readthrough'], rna_intermediates['5_cleavage'], rna_intermediates['3_cleavage']


def count_rna_intermediates_turbo(
    infile: str,
    polya_info: list,
    min_counts: int = 15,
    exclude_chr: set = None,
    threads: int = 10,
    ) -> pd.DataFrame:

    '''
    Counts the number of reads that map to each RNA intermediate in a set of RNA-seq data
    
    Parameters
    ----------
    infile : str
        the path to the BAM file
    polya_info : list
        list of tuples, each tuple is the polyA information: (chromosome, start, end, strand)
    min_counts : int, optional
        minimum number of counts for a gene to be considered
    exclude_chr : set
        set = None,
    threads : int, optional
        number of threads to use
    '''
    polya_info = np.array(polya_info)  # the polyA site information
    chrom = polya_info[:, 0]
    pas = polya_info[:, 1]
    strand = polya_info[:, 2]
    gene_id = polya_info[:, 3]

    with ProcessPoolExecutor(max_workers=threads) as e:
        chunksize = int(len(polya_info)/threads)
        results = e.map(
            count_rna_intermediates, 
            repeat(infile),
            chrom,
            pas,
            strand,
            gene_id,
            repeat(min_counts),
            repeat(exclude_chr),
            chunksize=chunksize
            )

    rna_intermediates = []
    for res in results:
        if res is not None:
            rna_intermediates.append(res)
    
    rna_intermediates = pd.DataFrame(rna_intermediates, columns=['gene_id', 'readthrough', '5_cleavage', '3_cleavage'])

    return rna_intermediates


def get_ir_ratio(infile: str, chrom: str, start: int, end: int, intron_id: str, strand: str):
    '''
    Counts the IR ratio for the given intron 
    (only test for Araport11 annotation)
    
    Parameters
    ----------
    infile : str
        the path to the BAM file
    chrom : str
        chromosome id
    start : int
        the start of the intron
    end : int
        the end of the intron
    intron_id : str
        the intron id, eg. AT1G01010.1_intron1
    strand : str
        "+" or "-"
    '''
    # infile = '/public/home/mowp/test/nanopore_test/20210324_col_nuclear/elongating_data/20210324_col_nuclear.elongating.bam'
    if intron_id.startswith('TE-'):
        # eg. TE-04844e9f-6b3f-4f3e-b4d0-435d42787293_intron1
        # compatible for TE transcripts
        gene_id = intron_id.split('_')[0]
    else:
        # eg. AT1G01010.1_intron1
        gene_id = intron_id.split('.')[0]
    intron_num = intron_id.split('intron')[1]

    STRAND_TO_BOOL = {'-': True, '+': False}
    totol_count, retain_count = 0, 0
    with pysam.AlignmentFile(infile, 'rb') as inbam:
        for read in inbam.fetch(chrom, start, end):
            read_gene_id = read.get_tag('gi')
            gene_strand = STRAND_TO_BOOL[strand]
            if read.is_reverse == gene_strand and read_gene_id == gene_id:  # 判断read是否属于这个基因
                if read.reference_start < start-10 and read.reference_end > end+10: # 判断read是否跨过这个intron
                    totol_count += 1
                    retain_intron_id = set(read.get_tag('ri').split(':'))
                    if intron_num in retain_intron_id:
                        retain_count += 1
    if totol_count > 0:
        return intron_id, totol_count, retain_count, retain_count/totol_count


def get_ir_ratio_turbo(inbam: str, repr_intron: str, threads: int = 64):
    '''This function takes a bam file and a bed file of introns and returns the intron retention ratio
    
    Parameters
    ----------
    inbam : str
        the input bam file
    repr_intron : str
        a file containing the representative introns
    threads : int, optional
        number of threads to use
    
    '''

    # repr_intron: mowp_scripts/pipelines/FLEP_seq_preprocessing_pipeline/ath_data/exon_intron_pos.repr.bed

    repr_intron = np.array(repr_intron)
    chrom = repr_intron[:, 0]
    start = repr_intron[:, 1]
    end = repr_intron[:, 2]
    strand = repr_intron[:, 3]
    intron_id = repr_intron[:, 4]

    with ProcessPoolExecutor(max_workers=threads) as e:
        chunksize = int(len(repr_intron) / threads)
        results = e.map(
            get_ir_ratio,
            repeat(inbam),
            chrom,
            start,
            end,
            intron_id,
            strand,
            chunksize=chunksize
        )
    
    ir_ratio = []
    for res in results:
        if res is not None:
            ir_ratio.append(res)

    ir_ratio = pd.DataFrame(ir_ratio, columns=['intron_id', 'totol_count', 'retain_count', 'IR_ratio'])
    
    return ir_ratio