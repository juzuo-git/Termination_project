import pysam
import os
import sys
from pathlib import Path
import argparse
from collections import Counter
import pyranges as pr
import numpy as np

"""Modifided 'or read.is_supplementary or read.is_secondary' in line114"""
parser = argparse.ArgumentParser(description="""\nSeperate 3'cleavage products from bam\n \
- 3'cleavage products\n \
- 5'cleavage products\n \
- read_through\n""", add_help=True)

required = parser.add_argument_group('required arguments')
required.add_argument('-s', '--sample_name', type=str, help='Sample name', required=True)
required.add_argument('-b', '--bam_in', help='Input bam file, form Nanopore basecalling', required=True)
required.add_argument('-c', '--min_cover', type=int, help='Input required depth', required=True)
required.add_argument('-o', '--out_fold', type=str, help='Out put result fold', required=True)

args = parser.parse_args()

def get_header_of_bam(bamin:str):
    with pysam.AlignmentFile(bamin, 'rb') as inbam:
        header_str = inbam.header
    return str(header_str)

def gnf(data, fo):
    with open(fo, "w") as fh:
        fh.write(data)
    print(f'>>{fo} had generated.<<')

def count_rna_intermediates(infile: str, chrom: str, pas: int, strand: str, gene_id: str, \
                            min_counts: int = 15, exclude_chr: set = None, \
                            reads_id_bool: bool = False, methods: str = "raw"):

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
     
    '''Condition 1: 3'cleavage
    `--- means some bases of (ATCG)`
    `PAS means polyadenylation signal`
    `||| means cover this region`
    5`-----------------PAS---------------------------------3`      `Picked Gene region`
                               ||||||||||||||||||||    
                               TGATGATGATGATGATGATG                `read`
    '''

    '''Condition 2: read through
    `--- means intergenic region`
    `PAS means polyadenylation signal`

    5`-----------------PAS---------------------------------3`    `Picked Gene region`
        ||||||||||||||||||||||||||||||||||||||||||||    
        TGATGATGATGATGPASATGTGATGATGATGATGATGATGTGAT             `read`
    '''

    '''Condition 3: 5'cleavage
    `--- means intergenic region`
    `PAS means polyadenylation signal`

    5`-----------------PAS---------------------------------3`    `Picked Gene region`
      |||||||||||||||||||    
      TGTGATGATGATGATATGA                                        `read`
    '''
    
    pas = int(pas)
    chrom = str(chrom)

    three_cleavage_reads, readthrough_reads, five_cleavage_reads = [], [], []
    if exclude_chr is not None and chrom in exclude_chr:
        exclude_chr = set(exclude_chr)
        if chrom in exclude_chr:
            return chrom

    strand_is_reverse = False if strand == '+' else True  # gene strand

    three_cleavage_reads_str = []
    five_cleavage_reads_str = []
    readthrough_reads_str = []
    
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

            if read.is_unmapped: ### or read.is_supplementary or read.is_secondary:
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

            # distance = three_end - pa_site
            distance3 = three_end - pa_site
            distance5 = five_end - pa_site
            
            if five_end >= pa_site:
                rna_intermediates['3_cleavage'] += 1
                three_cleavage_reads.append(read.qname)
                three_cleavage_reads_str.append(read.to_string())
                n += 1
            elif (distance3 > 50) and (distance5 < -50):
                rna_intermediates['readthrough'] += 1
                readthrough_reads.append(read.qname)
                readthrough_reads_str.append(read.to_string())
                n += 1
            elif -50 < distance3 < 50:
                if distance5 <= -50:
                    rna_intermediates['5_cleavage'] += 1
                    five_cleavage_reads.append(read.qname)
                    five_cleavage_reads_str.append(read.to_string())
                    n += 1
    if n >= min_counts:
        if reads_id_bool == True:
            if methods == "dic":
                return {"gene_id":gene_id, "three_cleavage_reads":three_cleavage_reads, \
                        'readthrough_reads':readthrough_reads, 'five_cleavage_reads':five_cleavage_reads}

            elif methods == "raw":
                return {"gene_id":gene_id, "three_cleavage_reads":three_cleavage_reads_str, \
                        'readthrough_reads':readthrough_reads_str, 'five_cleavage_reads':five_cleavage_reads_str}
                
        else:
            return gene_id, rna_intermediates['readthrough'], rna_intermediates['5_cleavage'], rna_intermediates['3_cleavage']

def get_pas_info(pas_bed_in):
    pas_bed = pas_bed_in
    # pas_bed = '/public/home/mowp/workspace/termination/cbRNA_pool/polya_sites/cbRNA.last_polya_cluster_summit.bed'
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
    # pas_info.head(3)
    pas_info = np.array(list(pas_info))
    print(pas_info[:3])

    return pas_info, pas_bed, single_pa_site_gene

def split_bam(out_fold:str, bam_in:str, sample_name:str, min_cover:int):
    
    
    
    ### get polyA site
    pas_info, _, _ = get_pas_info('/public/home/mowp/workspace/termination/cbRNA_pool/polya_sites/cbRNA.last_polya_cluster_summit.bed')

    pol3_black_list = ['AT4G04565', 'AT4G04595', 'AT1G05163','AT1G47720',
                    'AT2G04955', 'AT5G15022', 'AT3G06365', 'AT4G30993']

    polya_sites = [x for x in pas_info if x[3] not in pol3_black_list]
    polya_sites = np.array(polya_sites)
    try:
        Path(f'{out_fold}/{sample_name}').mkdir()
    except:
        print(f'>>{sample_name} had already exist!<<')
        
    # out_fold="/data/Zhaijx/lijz/Termination_project/xh_terminate_project/0_data/seperate_bam"

    # _three_cleavage_data, _readthrough_data, _five_cleavage_data = "", "", ""
    # _data_list = [_three_cleavage_data, _readthrough_data, _five_cleavage_data]
    _data_list = ["", "", ""]

    header_info = get_header_of_bam(bam_in)
    
    for x in range(0, 3, 1):
        _data_list[x]+= str(header_info)

    _three_cleavage_data, _readthrough_data, _five_cleavage_data = _data_list
    
    
    samfile = pysam.AlignmentFile(bam_in, "rb")
    three_reads = pysam.AlignmentFile(f"{out_fold}/{sample_name}/three_cl.bam", "wb", template=samfile)
    five_reads = pysam.AlignmentFile(f"{out_fold}/{sample_name}/five_cl.bam", "wb", template=samfile)
    readthrough_reads = pysam.AlignmentFile(f"{out_fold}/{sample_name}/readthrough_cl.bam", "wb", template=samfile)
    
    dic_reads_all = {'three_cleavage_reads':[], 'readthrough_reads':[], 'five_cleavage_reads':[]}
    for pas in polya_sites:
        chrom, pa_site, strand, gene_id = pas
        _dic = count_rna_intermediates(bam_in, chrom, pa_site, strand, gene_id, int(min_cover), None, True, "dic")
        if _dic:
            for x in ['three_cleavage_reads', 'readthrough_reads', 'five_cleavage_reads']:
                dic_reads_all[x].extend(_dic[x])
    
    dic_reads_all_v2 = {}
    for k in dic_reads_all.keys():
        dic_reads_all_v2[k] = {read_id:"" for read_id in dic_reads_all[k]}
    dic_reads_all = dic_reads_all_v2
    
    # del dictionary 
    dic_reads_all_v2 = ""

    for read in samfile.fetch():
        if read.qname in dic_reads_all['three_cleavage_reads']:
                three_reads.write(read)
        elif read.qname in dic_reads_all['readthrough_reads']:
                readthrough_reads.write(read)
        elif read.qname in dic_reads_all['five_cleavage_reads']:
                five_reads.write(read)
        else:
            pass
    
    three_reads.close()
    readthrough_reads.close()
    five_reads.close()
    samfile.close()
    
### main function###
if __name__ == "__main__":
    split_bam(args.out_fold, args.bam_in, args.sample_name, args.min_cover)
    print(f"Job Done\n>>{args.sample_name} had split.<<")
    