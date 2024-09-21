#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import matplotlib.patches as mp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import pysam
import seaborn as sns
from typing import Optional, Sequence, Tuple, Mapping
import pyBigWig

######
# plot gene model
######


def plot_gene_model(
    ax,
    gene_models: pd.DataFrame,
    fig_start: int,
    fig_end: int,
    gene_color: str = "k",
    y_space: int = 1,
    plot_gene_id: bool = True
):
    """plot gene model in the axis

    Args:
    -----
        ax (matplotlib.axes): An axis object to plot.

        gene_models (pd.DataFrame): A bed12 like DataFrame:

            chrom     start       end      gene_id score strand thickStart  thickEnd
                2  17922018  17924542  AT2G43110.1     0      +   17922066  17924291

            rgb blockCount              blockSizes                blockStarts y_pos
              0          6  361,196,73,50,114,372,  0,527,823,1802,1952,2152,     0

        fig_start (int): figure xlim start.

        fig_end (int): figure xlim end.

        gene_color (str, optional): gene color. Defaults to 'k'.

        y_space (int, optional): the spaces between gene in y direction. Defaults to 1.
    """
    ylim = 0  # ax ylim的下限
    height = 3  # gene model 高度
    y_space = y_space + height * 2
    for gene_model in gene_models.values:
        (
            chrom,
            chromStart,
            chromEnd,
            gene_id,
            _,
            strand,
            thickStart,
            thickEnd,
            _,
            blockCount,
            blockSizes,
            blockStarts,
            y_pos,
        ) = gene_model
        y_pos = -y_space * y_pos
        ylim = min(y_pos, ylim)

        # 数据类型转化
        chromStart = int(chromStart)
        chromEnd = int(chromEnd)
        thickStart = int(thickStart)
        thickEnd = int(thickEnd)
        blockSizes = np.fromstring(blockSizes, sep=",", dtype="int")
        blockStarts = np.fromstring(blockStarts, sep=",", dtype="int") + chromStart

        # 画转录起始位点及方向箭头
        small_relative = 0.06 * (fig_end - fig_start)  # 箭头突出部分相对长度
        arrowprops = dict(arrowstyle="-|>", connectionstyle="angle", color=gene_color)
        if strand == "+":
            ax.annotate(
                "",
                xy=(chromStart + small_relative, height * 2 + y_pos),
                xytext=(chromStart, y_pos),
                arrowprops=arrowprops,
            )
        else:
            ax.annotate(
                "",
                xy=(chromEnd - small_relative, height * 2 + y_pos),
                xytext=(chromEnd, y_pos),
                arrowprops=arrowprops,
            )

        line = mp.Rectangle(
            (chromStart, y_pos - height / 8),
            chromEnd - chromStart,
            height / 4,
            color=gene_color,
            linewidth=0,
        )  # 基因有多长这条线就有多长
        ax.add_patch(line)

        # draw exons
        cds_split = [thickStart]
        for exonstart, size in zip(blockStarts, blockSizes):
            exon = mp.Rectangle(
                    (exonstart, y_pos - height / 2),
                    size,
                    height,
                    color=gene_color,
                    linewidth=0,
                )
            ax.add_patch(exon)
            for exon_edge in [exonstart, exonstart + size]:
                if thickStart < exon_edge and exon_edge < thickEnd:
                    cds_split.append(exon_edge)
        cds_split.append(thickEnd)

        thick_list = [(cds_split[i], cds_split[i+1]) for i in range(0, len(cds_split), 2)]

        for thick_pair in thick_list:
            thick = mp.Rectangle((thick_pair[0], y_pos-height), thick_pair[1]-thick_pair[0], height*2, color=gene_color, linewidth=0)
            ax.add_patch(thick)

        # for exonstart, size in zip(blockStarts, blockSizes):
        #     if exonstart == chromStart and exonstart + size == chromEnd:
        #         utr_size = thickStart - chromStart
        #         utr = mp.Rectangle(
        #             (exonstart, y_pos - height / 2),
        #             utr_size,
        #             height,
        #             color=gene_color,
        #             linewidth=0,
        #         )
        #         ax.add_patch(utr)
        #         utr_size = chromEnd - thickEnd
        #         utr = mp.Rectangle(
        #             (thickEnd, y_pos - height / 2),
        #             utr_size,
        #             height,
        #             color=gene_color,
        #             linewidth=0,
        #         )
        #         ax.add_patch(utr)
        #         exon = mp.Rectangle(
        #             (thickStart, y_pos - height),
        #             thickEnd - thickStart,
        #             height * 2,
        #             color=gene_color,
        #             linewidth=0,
        #         )
        #         ax.add_patch(exon)

        #     elif exonstart + size <= thickStart:
        #         # 只有5'/ 3'UTR
        #         utr = mp.Rectangle(
        #             (exonstart, y_pos - height / 2),
        #             size,
        #             height,
        #             color=gene_color,
        #             linewidth=0,
        #         )
        #         ax.add_patch(utr)

        #     elif exonstart < thickStart and exonstart + size > thickEnd:
        #         # cds surrounded by utrs
        #         utr = mp.Rectangle((exonstart, y_pos-height/2), size, height, color=gene_color, linewidth=0)
        #         exon = mp.Rectangle((thickStart, y_pos-height), thickEnd-thickStart, height*2, color=gene_color, linewidth=0)
        #         ax.add_patch(utr)
        #         ax.add_patch(exon)

        #     elif exonstart < thickStart and exonstart + size > thickStart:
        #         # 带有5' / 3' UTR的exon
        #         utr_size = thickStart - exonstart
        #         utr = mp.Rectangle(
        #             (exonstart, y_pos - height / 2),
        #             utr_size,
        #             height,
        #             color=gene_color,
        #             linewidth=0,
        #         )
        #         exon = mp.Rectangle(
        #             (exonstart + utr_size, y_pos - height),
        #             size - utr_size,
        #             height * 2,
        #             color=gene_color,
        #             linewidth=0,
        #         )
        #         ax.add_patch(utr)
        #         ax.add_patch(exon)

        #     elif exonstart >= thickStart and exonstart + size <= thickEnd:
        #         # 普通exon
        #         exon = mp.Rectangle(
        #             (exonstart, y_pos - height),
        #             size,
        #             height * 2,
        #             color=gene_color,
        #             linewidth=0,
        #         )
        #         ax.add_patch(exon)

        #     elif exonstart < thickEnd and exonstart + size > thickEnd:
        #         # 带有3' / 5' UTR的exon
        #         utr_size = exonstart + size - thickEnd
        #         utr = mp.Rectangle(
        #             (thickEnd, y_pos - height / 2),
        #             utr_size,
        #             height,
        #             color=gene_color,
        #             linewidth=0,
        #         )
        #         exon = mp.Rectangle(
        #             (exonstart, y_pos - height),
        #             size - utr_size,
        #             height * 2,
        #             color=gene_color,
        #             linewidth=0,
        #         )
        #         ax.add_patch(utr)
        #         ax.add_patch(exon)

        #     elif exonstart >= thickEnd:
        #         # 只有3'/ 5'UTR
        #         utr = mp.Rectangle(
        #             (exonstart, y_pos - height / 2),
        #             size,
        #             height,
        #             color=gene_color,
        #             linewidth=0,
        #         )
        #         ax.add_patch(utr)

        if plot_gene_id:
            ax.annotate(
                gene_id, xy=((chromStart + chromEnd) / 2, y_pos + height * 1.5), ha="center"
            )

    # set ax
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.yaxis.set_major_locator(ticker.NullLocator())
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.xaxis.set_ticks_position("none")

    ax.set_ylim(ylim - height, height + y_space)
    ax.set_xlim(fig_start, fig_end)


def get_gene_model(chrom: str, start: int, end: int, bed_path: str, lambda_for_parse_geneId:Optional[str] = None) -> pd.DataFrame:
    """Get gene model information from the bed file.
        The bed file should be indexed by tabix.
        For details, see http://www.htslib.org/doc/tabix.html

    Args:
        chrom (str): chromosome id.

        start (int): the start of the region.

        end (int): the end of the region.

        bed_path (str): the PATH of the bed file.

        lambda_for_parse_geneId(str) : lambda expression for parse geneId. Defalut is x: x

    Returns:
        pd.DataFrame: A bed12 dataframe.
            Examples
            --------
            chrom     start       end      gene_id score strand thickStart  thickEnd
                2  17922018  17924542  AT2G43110.1     0      +   17922066  17924291

            rgb blockCount              blockSizes                blockStarts
              0          6  361,196,73,50,114,372,  0,527,823,1802,1952,2152,
    """
    tbx = pysam.TabixFile(bed_path)
    gene_models = pd.DataFrame(
        [gene_model.split("\t") for gene_model in tbx.fetch(chrom, start, end)],
        columns=[
            "chrom",
            "start",
            "end",
            "gene_id",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "rgb",
            "blockCount",
            "blockSizes",
            "blockStarts",
        ],
    )
    gene_models.sort_values(["start", "end"], inplace=True)
    if lambda_for_parse_geneId:
        lambda_for_parse_geneId_fc = eval(f"lambda {lambda_for_parse_geneId}")
        gene_models = gene_models.assign(gene_id = lambda df:df['gene_id'].map(lambda_for_parse_geneId_fc))
    return gene_models


######
# other function
######


def is_overlap(
    gene_a: Tuple[str, int, int], gene_b: Tuple[str, int, int], threshold: int = 0
) -> bool:
    """To judge whether two region is overlap.

    Args:
        gene_a (tuple): (chrom, start, end)
        gene_b (tuple): (chrom, start, end)
        threshold (int, optional): the minimum space between two region. Defaults to 0.

    Returns:
        bool: if overlap True, else False
    """
    minn = max(int(gene_a[1]), int(gene_b[1]))
    maxn = min(int(gene_a[2]), int(gene_b[2]))

    if maxn - minn >= -threshold:
        return True
    else:
        return False


def get_y_pos_discontinous(df: pd.DataFrame, gene_list: Optional[Sequence[str]] = None):
    """Get the y position of each region. Save the results in the df['y_pos'].

    Args:
        df (pd.DataFrame): a DataFrame must include first four columns:
            chrom, start, end, gene_id.

            Examples
            --------
            chrom    start      end    gene_id strand
                4  9105672  9106504  AT4G16100      +

        gene_list (set, optional): a set contain which gene to plot.
            When gene_list is None, plot all item in the df. Defaults to None.

    Returns:
        int: item counts in the y position.
    """

    # only plot reads/gene_model in the gene_list
    if gene_list is not None:
        df.query("gene_id in @gene_list", inplace=True)

    df.reset_index(drop=True, inplace=True)
    df["y_pos"] = df.index

    return len(df)


def filter_bam(
    df,
    strand=None,
    start_before=None,
    start_after=None,
    end_before=None,
    end_after=None,
):
    if strand is not None:
        df.query("strand == @strand", inplace=True)
    if start_before is not None:
        df.query("start <= @start_before", inplace=True)
    if start_after is not None:
        df.query("start >= @start_after", inplace=True)
    if end_before is not None:
        df.query("end <= @end_before", inplace=True)
    if end_after is not None:
        df.query("end >= @end_after", inplace=True)


def get_y_pos_continuous(df, gene_list=None, threshold=0):
    """Get the y position of each region. Save the results in the df['y_pos'].

    Args:
        df (pd.DataFrame): a DataFrame must include first four columns:
            chrom, start, end, gene_id.

            Examples
            --------
            chrom    start      end    gene_id strand
                4  9105672  9106504  AT4G16100      +

        gene_list (set, optional): a set contain which gene to plot.
            When gene_list is None, plot all item in the df. Defaults to None.

        threshold (int, optional): the minimum space between two region. Defaults to 0.

    Returns:
        int: item counts in the y position.
    """
    # initialization of y_pos columns
    df["y_pos"] = None

    read_list = []
    for index_id, item in enumerate(df.values):
        chrom, start, end, gene_id, *_ = item

        # only plot reads/gene_model in the gene_list
        if gene_list is not None and gene_id not in gene_list:
            df.drop(index_id, inplace=True)
            continue

        current_read = (chrom, start, end)
        is_add_read = False
        y_pos = -1
        for y_pos in range(len(read_list)):
            y_pos_read = read_list[y_pos]
            if not is_overlap(current_read, y_pos_read, threshold=threshold):
                read_list[y_pos] = current_read
                df.at[index_id, "y_pos"] = y_pos
                is_add_read = True
                break
        if not is_add_read:
            read_list.append(current_read)
            y_pos += 1
            df.at[index_id, "y_pos"] = y_pos

    return len(read_list)


def get_splice_sites(x):
    splice_sites = []
    for i in x:
        splice_sites.append(i[0])
        splice_sites.append(i[0]+i[1])
    
    return tuple(splice_sites[1:])


######
# plot bam
######


def find_exon(read):
    """Find exon position in the pysam.read.

    Args:
        read (pysam.read): <pysam.libcalignedsegment.AlignedSegment>

    Returns:
        zip(exon_start, exon_size): the exon position in the read.
            exon_start: exon start coordinate, 0-based
            exon_size: exon length

        Examples
        --------
        [(9105672, 457), (9106400, 142)]
    """
    BAM_CREF_SKIP = 3  # BAM_CREF_SKIP
    blockStart = []
    blockSize = []
    match_or_deletion = {
        0,
        2,
        7,
        8,
    }  # only M/=/X (0/7/8) and D (2) are related to genome position
    exon_start = read.reference_start
    length = 0
    for op, nt in read.cigartuples:
        if op in match_or_deletion:
            length += nt
        elif op == BAM_CREF_SKIP:
            blockStart.append(exon_start)
            blockSize.append(length)
            exon_start += length + nt
            length = 0
    blockStart.append(exon_start)
    blockSize.append(length)
    return zip(blockStart, blockSize)


def convert_bam(
    chrom,
    start,
    end,
    strand,
    infile,
    subsample=None,
    gene_list=None,
    method="continuous",
    filter_strand=None,
    start_before=None,
    start_after=None,
    end_before=None,
    end_after=None,
    force_tag_check: bool = True,
    tag_required_dict:Mapping[str, str]=dict(gene_id='gi', polya_len='pa', span_intron_count='sn', unsplice_count='rn', unsplice_intron='ri'),
    chrom_prefix='',
    color='#3183bd',
):
    """Conver pysam.reads into a DataFrame.

    Args:
    -----

        chrom (str): chromosome id

        start (int): start position of the region

        end (int): end position of the region

        strand (str): '+' or '-'

        infile (str): the PATH of the bam file

        subsample (int or float, optional): <number>|<frac> of items from axis to return.
            Defaults to None.

        gene_list (set, optional): a set contain which gene to plot.
            When gene_list is None, plot all item in the df. Defaults to None.

        method ('continuous' | '3_end' | '5_end')

    Returns:
        DataFrame: A dataframe contain bam reads information
            Examples
            --------
            Return: pd.DataFrame
            chrom    start      end    gene_id strand
                4  9105672  9106504  AT4G16100      +

                                         read_id  gap  polya_len
            37c91155-0ef5-47be-9e1b-d3a6753982f4  0.0          0

                        exon     y_pos
            [(9105672, 832)]      None
    """

    bam_data = []
    chrom = chrom_prefix+chrom
    with pysam.AlignmentFile(infile, "rb") as inbam:
        for read in inbam.fetch(chrom, start, end):
            gap = 0  # useless value

            try:
                polya_len = read.get_tag('pa')
            except KeyError:
                polya_len = 0

            try:
                gene_id = read.get_tag('gi')
            except KeyError:
                gene_id = 0

            try:
                span_intron_count = read.get_tag('sn')  # span_intron_num
                if type(span_intron_count) == str:
                    span_intron_count = len(span_intron_count.split(':'))
                unsplice_count = read.get_tag('rn')  # retention_intron_num
                unsplice_intron = read.get_tag('ri')  # retention_introns
                # only plot reads in gene_list
            except KeyError:
                span_intron_count, unsplice_count, unsplice_intron = None, None, None
            
            '''新增标记是否为novel isoform的tag，0为False，1为True'''
            try:
                novel_isoform = read.get_tag('ni')
            except KeyError:
                novel_isoform = 0

            '''新增标记是否chimeric read的tag，0为False，1为True'''
            try:
                chimeric_read = read.get_tag('cr')
            except KeyError:
                chimeric_read = 0


            if read.is_supplementary or read.is_unmapped:
                continue

            read_strand = "-" if read.is_reverse else "+"
            exon = find_exon(read)

            bam_data.append(
                (
                    read.reference_name,
                    read.reference_start,
                    read.reference_end,
                    gene_id,
                    read_strand,
                    read.query_name,
                    gap,
                    polya_len,
                    list(exon),
                    span_intron_count,
                    unsplice_count,
                    unsplice_intron,
                    novel_isoform,
                    chimeric_read,
                )
            )
    bam_data = pd.DataFrame(
        bam_data,
        columns=[
            "chrom",
            "start",
            "end",
            "gene_id",
            "strand",
            "read_id",
            "gap",
            "polya_len",
            "exon",
            "span_intron_count",
            "unsplice_count",
            "unsplice_intron",
            "novel_isoform",
            "chimeric_read",
        ],
    )

    bam_data['color'] = color  # read color

    # 进行subsample
    if subsample is not None:
        if type(subsample) is int:
            bam_data = bam_data.sample(n=subsample, random_state=42)
        elif type(subsample) is float and subsample < 1:
            bam_data = bam_data.sample(frac=subsample, random_state=42)

    if method == "continuous":
        # 类似于igv那样连续的排
        # add polya
        filter_bam(
            bam_data,
            strand=filter_strand,
            start_before=start_before,
            start_after=start_after,
            end_before=end_before,
            end_after=end_after,
        )
        
        bam_data["polya_len"] = bam_data["polya_len"].map(
            lambda x: int(x) if x >= 15 else 0
        )
        bam_data = bam_data.assign(
            start=np.select(
                [bam_data["strand"] == "+", bam_data["strand"] == "-"],
                [bam_data["start"], bam_data["start"] - bam_data["polya_len"]],
                -1,
            ),
            end=np.select(
                [bam_data["strand"] == "+", bam_data["strand"] == "-"],
                [bam_data["end"] + bam_data["polya_len"], bam_data["end"]],
                -1,
            ),
        )
        if strand == "+":
            bam_data.sort_values(["start", "end"], inplace=True)
        else:
            bam_data.sort_values(["end", "start"], ascending=False, inplace=True)
        bam_data.reset_index(drop=True, inplace=True)
        get_y_pos_continuous(bam_data, gene_list=gene_list)

    elif method == "5_end":
        if strand == "+":
            bam_data.sort_values(["start", "end"], inplace=True)
        else:
            bam_data.sort_values(["end", "start"], ascending=False, inplace=True)

        # add polya
        bam_data["polya_len"] = bam_data["polya_len"].map(
            lambda x: int(x) if x >= 15 else 0
        )
        bam_data = bam_data.assign(
            start=np.select(
                [bam_data["strand"] == "+", bam_data["strand"] == "-"],
                [bam_data["start"], bam_data["start"] - bam_data["polya_len"]],
                -1,
            ),
            end=np.select(
                [bam_data["strand"] == "+", bam_data["strand"] == "-"],
                [bam_data["end"] + bam_data["polya_len"], bam_data["end"]],
                -1,
            ),
        )

        filter_bam(
            bam_data,
            strand=filter_strand,
            start_before=start_before,
            start_after=start_after,
            end_before=end_before,
            end_after=end_after,
        )
        get_y_pos_discontinous(bam_data, gene_list=gene_list)

    elif method == "3_end":
        if strand == "+":
            bam_data.sort_values(["end", "start"], inplace=True)
        else:
            bam_data.sort_values(["start", "end"], ascending=False, inplace=True)

        # add polya
        bam_data["polya_len"] = bam_data["polya_len"].map(
            lambda x: int(x) if x >= 15 else 0
        )
        bam_data = bam_data.assign(
            start=np.select(
                [bam_data["strand"] == "+", bam_data["strand"] == "-"],
                [bam_data["start"], bam_data["start"] - bam_data["polya_len"]],
                -1,
            ),
            end=np.select(
                [bam_data["strand"] == "+", bam_data["strand"] == "-"],
                [bam_data["end"] + bam_data["polya_len"], bam_data["end"]],
                -1,
            ),
        )

        filter_bam(
            bam_data,
            strand=filter_strand,
            start_before=start_before,
            start_after=start_after,
            end_before=end_before,
            end_after=end_after,
        )
        get_y_pos_discontinous(bam_data, gene_list=gene_list)

    elif method == "gene_id":  # 按照gene_id 3' end排序
        if strand == "+":
            bam_data.sort_values(["gene_id", "end", "start"], inplace=True)
        else:
            bam_data.sort_values(
                ["gene_id", "start", "end"], ascending=False, inplace=True
            )

        # add polya
        bam_data["polya_len"] = bam_data["polya_len"].map(
            lambda x: int(x) if x >= 15 else 0
        )
        bam_data = bam_data.assign(
            start=np.select(
                [bam_data["strand"] == "+", bam_data["strand"] == "-"],
                [bam_data["start"], bam_data["start"] - bam_data["polya_len"]],
                -1,
            ),
            end=np.select(
                [bam_data["strand"] == "+", bam_data["strand"] == "-"],
                [bam_data["end"] + bam_data["polya_len"], bam_data["end"]],
                -1,
            ),
        )
        get_y_pos_discontinous(bam_data, gene_list=gene_list)

    elif method == "spliced":  # 只画完全剪切的reads
        bam_data.query("span_intron_count > 0 and unsplice_count == 0", inplace=True)
        if strand == "+":
            bam_data.sort_values(["end", "start"], inplace=True)
        else:
            bam_data.sort_values(["start", "end"], ascending=False, inplace=True)

        # add polya
        bam_data["polya_len"] = bam_data["polya_len"].map(
            lambda x: int(x) if x >= 15 else 0
        )
        bam_data = bam_data.assign(
            start=np.select(
                [bam_data["strand"] == "+", bam_data["strand"] == "-"],
                [bam_data["start"], bam_data["start"] - bam_data["polya_len"]],
                -1,
            ),
            end=np.select(
                [bam_data["strand"] == "+", bam_data["strand"] == "-"],
                [bam_data["end"] + bam_data["polya_len"], bam_data["end"]],
                -1,
            ),
        )
        filter_bam(
            bam_data,
            strand=filter_strand,
            start_before=start_before,
            start_after=start_after,
            end_before=end_before,
            end_after=end_after,
        )
        get_y_pos_discontinous(bam_data, gene_list=gene_list)

    elif method == "partially_spliced":  # 只画不完全剪切的reads
        bam_data.query(
            "span_intron_count > 0 and unsplice_count > 0 and unsplice_count != span_intron_count ",
            inplace=True,
        )
        if strand == "+":
            bam_data.sort_values(["end", "start"], inplace=True)
        else:
            bam_data.sort_values(["start", "end"], ascending=False, inplace=True)

        # add polya
        bam_data["polya_len"] = bam_data["polya_len"].map(
            lambda x: int(x) if x >= 15 else 0
        )
        bam_data = bam_data.assign(
            start=np.select(
                [bam_data["strand"] == "+", bam_data["strand"] == "-"],
                [bam_data["start"], bam_data["start"] - bam_data["polya_len"]],
                -1,
            ),
            end=np.select(
                [bam_data["strand"] == "+", bam_data["strand"] == "-"],
                [bam_data["end"] + bam_data["polya_len"], bam_data["end"]],
                -1,
            ),
        )
        filter_bam(
            bam_data,
            strand=filter_strand,
            start_before=start_before,
            start_after=start_after,
            end_before=end_before,
            end_after=end_after,
        )
        get_y_pos_discontinous(bam_data, gene_list=gene_list)

    elif method == "unspliced":  # 只画不完全剪切的reads
        bam_data.query(
            "span_intron_count > 0 and unsplice_count == span_intron_count",
            inplace=True,
        )
        if strand == "+":
            bam_data.sort_values(["end", "start"], inplace=True)
        else:
            bam_data.sort_values(["start", "end"], ascending=False, inplace=True)

        # add polya
        bam_data["polya_len"] = bam_data["polya_len"].map(
            lambda x: int(x) if x >= 15 else 0
        )
        bam_data = bam_data.assign(
            start=np.select(
                [bam_data["strand"] == "+", bam_data["strand"] == "-"],
                [bam_data["start"], bam_data["start"] - bam_data["polya_len"]],
                -1,
            ),
            end=np.select(
                [bam_data["strand"] == "+", bam_data["strand"] == "-"],
                [bam_data["end"] + bam_data["polya_len"], bam_data["end"]],
                -1,
            ),
        )
        filter_bam(
            bam_data,
            strand=filter_strand,
            start_before=start_before,
            start_after=start_after,
            end_before=end_before,
            end_after=end_after,
        )
        get_y_pos_discontinous(bam_data, gene_list=gene_list)

    elif method == "ir":
        bam_data['exon1'] = bam_data['exon'].map(get_splice_sites)
        if strand == "+":
            bam_data.sort_values(['novel_isoform', 'exon1', "end", "start"], inplace=True)
        else:
            bam_data.sort_values(['exon1', "start", "end"], ascending=False, inplace=True)
            bam_data.sort_values(['novel_isoform'], inplace=True)
        
        bam_data['color'] = bam_data.apply(lambda x: '#A6CEE3' if x.novel_isoform == 0 else x.color, axis=1)

        # add polya
        bam_data["polya_len"] = bam_data["polya_len"].map(
            lambda x: int(x) if x >= 15 else 0
        )
        bam_data = bam_data.assign(
            start=np.select(
                [bam_data["strand"] == "+", bam_data["strand"] == "-"],
                [bam_data["start"], bam_data["start"] - bam_data["polya_len"]],
                -1,
            ),
            end=np.select(
                [bam_data["strand"] == "+", bam_data["strand"] == "-"],
                [bam_data["end"] + bam_data["polya_len"], bam_data["end"]],
                -1,
            ),
        )

        filter_bam(
            bam_data,
            strand=filter_strand,
            start_before=start_before,
            start_after=start_after,
            end_before=end_before,
            end_after=end_after,
        )
        get_y_pos_discontinous(bam_data, gene_list=gene_list)

    elif method == "chimeric":
        if strand == "+":
            
            bam_data.sort_values(["chimeric_read", "end", "start"], inplace=True)
        else:
            bam_data['chimeric_read'] = -1 * bam_data['chimeric_read']
            bam_data.sort_values(["chimeric_read", 'start', 'end'], ascending=False, inplace=True)
        
        bam_data['color'] = bam_data.apply(lambda x: '#A6CEE3' if x.chimeric_read == 0 else x.color, axis=1)

        # add polya
        bam_data["polya_len"] = bam_data["polya_len"].map(
            lambda x: int(x) if x >= 15 else 0
        )
        bam_data = bam_data.assign(
            start=np.select(
                [bam_data["strand"] == "+", bam_data["strand"] == "-"],
                [bam_data["start"], bam_data["start"] - bam_data["polya_len"]],
                -1,
            ),
            end=np.select(
                [bam_data["strand"] == "+", bam_data["strand"] == "-"],
                [bam_data["end"] + bam_data["polya_len"], bam_data["end"]],
                -1,
            ),
        )

        filter_bam(
            bam_data,
            strand=filter_strand,
            start_before=start_before,
            start_after=start_after,
            end_before=end_before,
            end_after=end_after,
        )
        get_y_pos_discontinous(bam_data, gene_list=gene_list)

    return bam_data


def plot_bam(
    ax,
    bam_data,
    fig_start,
    fig_end,
    read_color="#3183bd",  # 已经无用
    polya_color="lightcoral",
    y_space=1.5,  # the space between reads in yaxis
    gene_list=None,
    pal=None,
    height: int = 0.5, # reads的高度
):
    """plot bam information to the ax
    ### modified 1027, 1030, 1033 these lines
    Args:
    -----
        ax (matplotlib.axes): An axis object to plot

        bam_data (pd.DataFrame): A dataframe contain bam reads information
            Examples
            --------
            Return: pd.DataFrame
            chrom    start      end    gene_id strand
                4  9105672  9106504  AT4G16100      +

                                         read_id  gap  polya_len
            37c91155-0ef5-47be-9e1b-d3a6753982f4  0.0          0

                        exon     y_pos
            [(9105672, 832)]         0

        fig_start (int): figure xlim start.

        fig_end (int): figure xlim end.

        read_color (str, optional): read color. Defaults to '#5D93C4'.

        polya_color (str, optional): polya color. Defaults to 'lightcoral'.

        y_space (int, optional): the spaces between reads in y direction. Defaults to 1.5.

        pal (seaborn.palettes._ColorPalette, optional)
    """
    if pal is None:
        pal = sns.color_palette("Paired")

    gene_color = {}
    gene_color_index = 0

    ylim = 0  # the start position of yaxis
    for item in bam_data.itertuples():
        (
            chrom,
            start,
            end,
            gene_id,
            strand,
            read_id,
            gap,
            polya_len,
            exon,
            *_,
            ypos,
        ) = item[1:]
        ypos = -y_space * ypos
        ylim = min(ypos, ylim)
        line = mp.Rectangle(
            (start, ypos - height / 4),
            end - start,
            height / 2,
            color="#A6A6A6",
            linewidth=0,
        )
        ax.add_patch(line)
        for block_start, block_size in exon:

            # set gene color
            # index为色板中的序列编号
            # 如果gene_list存在 则给gene_list里面都基因上色
            if gene_list is not None:
                if gene_id in gene_list:
                    if gene_id not in gene_color:
                        gene_color[gene_id] = gene_color_index
                        gene_color_index += 1
                    read_color_ = read_color ### pal[gene_color[gene_id]]
                else:
                    # 不在gene_list里面的reads都设置成灰色
                    read_color_ = read_color ### "#5D5D5D"
            # 如果没有则统一颜色
            else:
                read_color_ = read_color ### item.color

                # TODO: 不同链的reads不同颜色
            exon = mp.Rectangle(
                (block_start, ypos - height),
                block_size,
                height * 2,
                color=read_color_,
                linewidth=0,
            )
            ax.add_patch(exon)

        # plot polya
        if polya_len > 15:
            if strand == "+":
                polya_tail = mp.Rectangle(
                    (block_start + block_size, ypos - height),
                    polya_len,
                    height * 2,
                    color=polya_color,
                    linewidth=0,
                )
            else:
                polya_tail = mp.Rectangle(
                    (start, ypos - height),
                    polya_len,
                    height * 2,
                    color=polya_color,
                    linewidth=0,
                )
            ax.add_patch(polya_tail)

    ax.set_ylim(ylim * 1.1, height + y_space)
    # ax.set_xlim(fig_start, fig_end)


def set_ax(ax, plot_xaxis=False):
    """Set axes

    Args:
        ax (matplotlib.axes): An axis object to plot

        plot_xaxis (bool, optional): if True plot x_ticks and x_locator.
            Defaults to False.
    """
    ax.spines["right"].set_visible(False)  # 去线
    ax.spines["left"].set_visible(False)  # 去线
    ax.spines["top"].set_visible(False)  # 去线
    ax.yaxis.set_major_locator(ticker.NullLocator())  # 去y数字

    if not plot_xaxis:
        ax.xaxis.set_major_locator(ticker.NullLocator())  # 去x数字
        ax.xaxis.set_ticks_position("none")  # 去x刻度


def plot_bw_track(ax, track_data, start, end, track_color: str = None, data_range: tuple = (None, None)):
    '''This function takes in a matplotlib axis, a dataframe of track data, a start and end time, and a
    color, and plots the track on the axis
    
    Parameters
    ----------
    ax
        the axis to plot the track on
    track_data
        a list of tuples, each tuple is a point in the track
    start
        the start time of the track
    end
        the end of the track
    track_color
        the color of the track
    
    '''
    ax.fill_between(np.linspace(start, end, end-start), track_data, color=track_color)
    # set ax
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.yaxis.set_major_locator(ticker.NullLocator())
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.xaxis.set_ticks_position("none")
    # ax.set_xlim(start, end)
    if not(None in data_range):
        ax.set_ylim(data_range)


######
# main Class
######


class IGV:
    """A IGV class

    Usage:
    ------
        chrom, start, end, strand, force_tag_check = '5', 5371627, 5375616, '+', True

        igv_plot = igv.IGV(chrom, start, end, strand=strand)

        araport11_isoform_path = '/public/home/mowp/db/Arabidopsis_thaliana/representative_gene_model/araport11.representative.gene_model.bed.gz'
        igv_plot.add_gene_model(
            araport11_isoform_path,
            gene_list = {'AT5G16440.1', 'AT5G16450.1'},
        )

        infile = '/public/home/mowp/workspace/termination/cbRNA_pool/elongating_data/elongating_data.bam'
        igv_plot.add_bam(
            infile,
            gene_list = {'AT5G16440', 'AT5G16450'},
            read_color='r',
            bam_title='one bam'
        )

        igv_plot.plot(height=4, width=8)
    ------
    tag_required_dict:Mapping[str, str]:
        tag required. Key is the variable name, and value is corresponding tag.
        if force_tag_check is True, all tag should be stored in bam; if not, will assign values to variables in order, and when error is triggered, left variables will be ignored.
    """

    def __init__(
        self,
        chrom,
        start,
        end,
        strand="+",
        no_x = False,
        tag_required_dict:Mapping[str, str]=dict(gene_id='gi', polya_len='pa', span_intron_count='sn', unsplice_count='rn', unsplice_intron='ri'),
        gene_id_use:str = None
    ):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.bam_list = []
        self.tagRequiered = tag_required_dict
        self.colorList = []
        self.bamTitleList = []
        self.gene_models = None
        self.gene_model_ylim = 0
        self.bw_track_list = []

        self.gene_model_n = 0
        self.bw_tracks_n = 0
        self.bam_list_n = 0
        self.no_x = no_x
        self.gene_id_use = gene_id_use
        
    def add_bam(
        self,
        bam_path,
        gene_list=None,
        subsample=None,
        method="continuous",
        filter_strand=None,
        start_before=None,
        start_after=None,
        end_before=None,
        end_after=None,
        read_color='#3183bd',
        bam_title='',
        chrom_prefix='',
        force_tag_check: bool = True,
    ):

        '''
        chrom_prefix:
            the chromosome name, such as "chr", "Chr", default is ""
        force_tag_check:
            if True, reads which don't have required tags will be ignored
        '''

        bam_data = convert_bam(
            self.chrom,
            self.start,
            self.end,
            self.strand,
            bam_path,
            subsample=subsample,
            gene_list=gene_list,
            method=method,
            filter_strand=filter_strand,
            start_before=start_before,
            start_after=start_after,
            end_before=end_before,
            end_after=end_after,
            force_tag_check=force_tag_check,
            tag_required_dict=self.tagRequiered,
            chrom_prefix=chrom_prefix,
        )
        self.bam_list.append(bam_data)
        self.gene_list = gene_list
        self.colorList.append(read_color)
        self.bamTitleList.append(bam_title)
        self.bam_list_n += 1

        return bam_data

    def add_bw(
        self,
        infile: str,
        color='#3183bd',
        chrom_prefix='',
        data_range: tuple = (None, None),
    ):
        bw = pyBigWig.open(infile)
        chrom = chrom_prefix + self.chrom
        bw_data = bw.values(chrom, self.start, self.end)
        bw_data = np.nan_to_num(bw_data)

        self.bw_track_list.append((bw_data, color, data_range))
        self.bw_tracks_n += 1


    def add_gene_model(self, anno:str, gene_list=None, auto_create_gz=False, lambda_for_parse_geneId:Optional[str]=None, plot_gene_id: bool = True):
        """
        auto_create_gz: if provided bed is not compressed, the bgzip will be used to created the gz format file, corresponding index will also be built.
        lambda_for_parse_geneId: If not None, it will be used to parse the `gene_id` column in bed. None is equal to `x:x`
        """
        if anno.endswith('.gz'):
            self.anno = anno
        else:
            self.anno = anno + '.gz'
            if not os.path.exists(self.anno):
                assert auto_create_gz, "WRONG suffix, and <auto_create_gz> is False"
                print(f"WARNING: Try to create file {self.anno} and build index")
                os.system(f'sort -k1,1 -k2,2n {anno} > temp____ && bgzip -c  temp____ > {self.anno} && tabix -p bed {self.anno}')
        self.gene_models = get_gene_model(self.chrom, self.start, self.end, self.anno, lambda_for_parse_geneId=lambda_for_parse_geneId)
        self.gene_model_ylim = get_y_pos_continuous(
            self.gene_models, gene_list=gene_list, threshold=8
        )  # 不重叠y轴的基因数目
        self.plot_gene_id = plot_gene_id
        self.gene_model_n = 1

    def paras_out(self):
        return self.colorList
        
    def plot(
        self,
        b,
        height=5,
        width=10,
        gene_track_height=0.5,
        bw_track_height=1,
        bam_track_height=4,
        gene_color="k",
        extend_xlim_start=False,
        extend_xlim_end=False,
        polya_site=None,
        read_height: int = 0.5,
        relative_pos: bool = True,
    ):
        '''It plots the data.
        
        Parameters
        ----------
        height, optional
            the height of the plot
        width, optional
            the width of the plot
        gene_track_height
            height of the gene track
        bw_track_height, optional
            height of the bigwig track
        bam_track_height, optional
            the height of the bam track
        gene_color, optional
            color of the gene track
        extend_xlim_start, optional
            If True, the x-axis will be extended to the left by the length of the left-most read.
        extend_xlim_end, optional
            If True, the x-axis will be extended to the right by the length of the right-most read.
        polya_site
            the position of the polyA site. If not provided, it will be plot the polyA site.
        read_height : int
            the height of the reads in the bam track
        relative_pos : bool, optional
            bool = True
        
        '''

        gridspec_kw = [gene_track_height]*self.gene_model_n + [bw_track_height]*self.bw_tracks_n + [bam_track_height]*self.bam_list_n
        nrows = self.gene_model_n + self.bw_tracks_n + self.bam_list_n

        fig, axes = plt.subplots(
            figsize=(width, height), 
            nrows=nrows,
            gridspec_kw={"height_ratios": gridspec_kw},
            sharex=True,
        )

        if nrows == 1:
            axes = [axes]

        axes_index = 0
        # plot gene_model
        if self.gene_models is not None:
            ax = axes[0]
            plot_gene_model(
                ax,
                self.gene_models,
                self.start,
                self.end,
                gene_color=gene_color,
                y_space=6,
                plot_gene_id=self.plot_gene_id,
            )
            ax.patch.set_alpha(0)
            axes_index += 1
        
        # plot bw_tracks
        for bw_data, color, data_range in self.bw_track_list:
            plot_bw_track(axes[axes_index], bw_data, self.start, self.end, color, data_range)
            axes[axes_index].patch.set_alpha(0)
            axes_index += 1
            

        # plot bam files
        for bam_data, read_color, bamTitle in zip(self.bam_list, self.colorList, self.bamTitleList):
            ax = axes[axes_index]

            plot_bam(ax, bam_data, self.start, self.end, gene_list=self.gene_list, read_color=read_color, height=read_height)
            plot_xaxis = axes_index == nrows - 1
            set_ax(ax, plot_xaxis=plot_xaxis)

            # add polya site
            end_ = None
            end_dict = {}
            
            if polya_site is not None:
                tbx = pysam.TabixFile(polya_site)
                for item in tbx.fetch(self.chrom, self.start, self.end):
                    chrom_, start_, end_, gene_id_, _, strand_, _ = item.split("\t")
                    gene_id_ = gene_id_.split("_")[0]
                    end_dict[gene_id_] = end_
                    if self.gene_list is not None:
                        if gene_id_ in self.gene_list:
                            ax.axvline(int(end_), ls="--", color="#555555")
                    else:
                        ax.axvline(int(end_), ls="--", color="#555555")
            plt.sca(ax)
            plt.text(1.0, 0.5, bamTitle, transform=ax.transAxes,
                fontsize=12, va='center', ha='left')
            ax.patch.set_alpha(0)
            axes_index += 1

        if relative_pos:
            # set last xaxis
            maxn = self.end
            minn = self.start
            for bam_data in self.bam_list:
                if len(bam_data) == 0:
                    continue

                minn_ = min(bam_data["start"])
                maxn_ = max(bam_data["end"])

                if extend_xlim_start:
                    if self.strand == "+":
                        minn = min(minn_, minn)
                    else:
                        maxn = max(maxn_, maxn)

                if extend_xlim_end:
                    if self.strand == "+":
                        maxn = max(maxn_, maxn)
                    else:
                        minn = min(minn_, minn)

            if maxn - minn > 400:
                step = (maxn - minn) // 400 * 100  # 坐标轴步长
            else:
                step = (maxn - minn) // 40 * 10

            ax = axes[nrows - 1]
            
            # b = np.array([-200, 0, 500, 1000, 1500])
            if self.gene_id_use != None:
                pas_s = end_dict[self.gene_id_use]
            else:
                _end_list = list(end_dict.keys())
                if self.strand == "+":
                    pas_s = _end_list[0]
                else:
                    pas_s = _end_list[-1] 

            if self.strand == "+":
                xticks = int(pas_s) + b
                ax.set_xticks(xticks)
                ax.set_xticklabels(b)
                left_site = int(pas_s) + b[0]
                right_site = int(pas_s) + b[-1]
                # left_site = max(int(pas_s) + b[0], minn)
                # right_site = min(int(pas_s) + b[-1], maxn)
                ax.set_xlim(left_site, right_site)
            else:
                xticks = int(pas_s) - b ### np.arange(maxn, minn - step, -step)
                ax.set_xticks(xticks)
                ax.set_xticklabels(b)
                # left_site = max(int(pas_s) - b[0], minn)
                # right_site = min(int(pas_s) - b[-1], maxn)
                left_site = int(pas_s) - b[0]
                right_site = int(pas_s) - b[-1]
                ax.set_xlim(minn, maxn)
                ax.set_xlim(right_site, left_site)
                ax.invert_xaxis()

        # else:
        #     # plot genome coordinate
        #     start = (self.start + 100_000) // 100_000 * 100_000
        #     step = (self.end - start) // 3
        #     step = step // 1000 * 1000
        #     ax = axes[nrows - 1]
        #     xticks = np.arange(start, self.end, step)
        #     ax.set_xticks(xticks)
        #     if start > 10_000:
        #         ax.set_xticklabels(map(lambda x: f'{x:,d} kb', xticks//1000))
        #     ax.set_xlim(self.start, self.end)
        
        ax.tick_params(axis=u'both', which=u'both',length=0)
        if self.no_x == True:
            plt.xticks(color = "white")
        return axes, pas_s, left_site, right_site


from jpy_tools.otherTools import F

class FigConcate(object):
    def __init__(self, fig):
        """
        初始化函数，将输入的Matplotlib图形转换为numpy数组，并保存在self.figAr中

        参数：
        fig：Matplotlib图形对象
        """
        self.fig = fig
        self.figAr = self.figureToArray(fig)

    def __or__(self, other):
        """
        重载 | 运算符，将两个FigConcate对象沿着水平方向拼接

        参数：
        other：另一个FigConcate对象

        返回：
        拼接后的FigConcate对象
        """
        ar_concate = self.padAndConcate([self.figAr, other.figAr], axis=1)
        fig = plt.figure()
        plt.imshow(ar_concate)
        plt.axis('off')
        plt.close()
        fig_concate = FigConcate(fig)
        fig_concate.figAr = ar_concate
        return fig_concate
    
    def __truediv__(self, other):
        """
        重载 / 运算符，将两个FigConcate对象沿着垂直方向拼接

        参数：
        other：另一个FigConcate对象

        返回：
        拼接后的FigConcate对象
        """
        ar_concate = self.padAndConcate([self.figAr, other.figAr], axis=0)
        fig = plt.figure()
        plt.imshow(ar_concate)
        plt.axis('off')
        plt.close()
        fig_concate = FigConcate(fig)
        fig_concate.figAr = ar_concate
        return fig_concate
    
    def show(self, figsize=(10, 10), figure_name:str = "tmp.svg"):
        fig = plt.figure(figsize=figsize)
        plt.imshow(self.figAr)
        plt.axis('off')
        plt.savefig(figure_name)
        plt.close()
        return fig

    @staticmethod
    def figureToArray(fig):
        """
        将Matplotlib图形转换为numpy数组

        参数：
        fig：Matplotlib图形对象

        返回：
        numpy数组
        """
        from io import BytesIO
        import PIL

        buf = BytesIO()
        fig.savefig(buf, format="png", bbox_inches='tight')
        buf.seek(0)
        img = PIL.Image.open(buf)
        return np.array(img)

    @staticmethod
    def padAndConcate(ls_arrs, axis):
        """
        将多个numpy数组沿着指定轴拼接，并进行补齐

        参数：
        ls_arrs：包含多个numpy数组的列表
        axis：指定拼接的轴，0表示垂直方向，1表示水平方向

        返回：
        拼接后的numpy数组
        """
        from functools import reduce

        def concatenate_along_axis(arr1, arr2, axis=axis):
            # 获取 arr1 和 arr2 的形状
            paddingAxis = 1-axis
            shape1, shape2 = arr1.shape, arr2.shape
            n1 = shape1[paddingAxis]
            n2 = shape2[paddingAxis]
            max_n = max(n1, n2)

            # 确定需要补齐的维度和补齐的大小
            pad_shape1 = [(0, 0)] * len(shape1)
            pad_shape2 = [(0, 0)] * len(shape1)

            padding_n1 = max_n - n1
            padding_n2 = max_n - n2

            pad_shape1[paddingAxis] = (padding_n1//2, padding_n1 - padding_n1//2)
            pad_shape2[paddingAxis] = (padding_n2//2, padding_n2 - padding_n2//2)

            # 补齐 arr1 和 arr2
            arr1 = np.pad(arr1, pad_shape1, mode='constant', constant_values=255)
            arr2 = np.pad(arr2, pad_shape2, mode='constant', constant_values=255)

            # 沿着指定轴拼接 arr1 和 arr2
            return np.concatenate((arr1, arr2), axis=axis)

        return reduce(concatenate_along_axis, ls_arrs)

class FigConcateWrap(object):
    def __init__(self):
        self.lsFig = []

    def addFig(self, FigConcate):
        self.lsFig.append(FigConcate)

    def wrapAndGenerate(self, wrap=4) -> FigConcate:
        from more_itertools import chunked
        from functools import reduce
        _ls = []
        if wrap > 1:
            for ls in chunked(self.lsFig, wrap):
                _ls.append(reduce(lambda x, y: x|y, ls))
            return reduce(lambda x, y: x/y, _ls)
        else:
            return reduce(lambda x, y: x/y, self.lsFig)

def cld(df, p, lvl_order=None):
    '''The function `cld` performs a post-hoc analysis using the multcomp package in R to determine significant differences between groups based on a given p-value threshold.

    Parameters
    ----------
    df
        The parameter `df` is a pandas DataFrame that contains the data for which you want to perform the cld (comparisons of least significant difference) analysis.
    p
        The parameter "p" represents the significance level for the p-values. It is used to determine which p-values are considered statistically significant.
    lvl_order
        The `lvl_order` parameter is used to specify the order of levels in the output. It is a list or tuple that contains the levels in the desired order. If `lvl_order` is not provided, the levels will be sorted in ascending order.

    Returns
    -------
        a dictionary where the keys are the names of the comparisons and the values are the corresponding letters.

    '''
    from rpy2 import robjects as ro
    from rpy2.robjects.packages import importr
    from .rTools import py2r

    multcomp = importr('multcomp')
    R = ro.r
    if lvl_order is None:
        lvl_order = (df['row'] >> F(set) | df['col'] >> F(set)) >> F(sorted)

    signif = (df['p-adj'] < p).values
    mycomps = df[['row', 'col']].values
    res = multcomp.insert_absorb(R.c(*signif.tolist()), decreasing=False, comps=py2r(mycomps), lvl_order=lvl_order)
    dt_cld = R("""\(x){
    x$Letters
    }
    """)(res) >> F(lambda _: {x:y for x,y in zip(_.names, _)})
    return dt_cld


def bxt(pas_out, _igv_plot, strand, input_regions, num_used:int = 0, svg:str = "boxplot.tmp.svg"):
    if strand == "-":
        readthrough_len_ = int(pas_out) - np.array(_igv_plot.bam_list[num_used]['start'])
    else:
        readthrough_len_ = np.array(_igv_plot.bam_list[num_used]['end']) - int(pas_out)

    fig, ax = plt.subplots(figsize=(10, 1.6))

    flierprops = dict(marker='o', markersize=3, markerfacecolor='k', markeredgecolor='none')

    ax.boxplot(readthrough_len_, widths=.5, showfliers = True, flierprops=flierprops, vert=False)
    plt.xlabel('Readthrough distance (nt)')

    sns.despine(top=True, right=True, left=True)

    b = input_regions
    # b = np.array([-200, 0, 500, 1000, 1500])
    xticks_str = [str(x) for x in b]
    plt.xticks(b, xticks_str )

    plt.yticks([])

    ax.axvline(0, ls="--", color="#555555")
    # ax.margins(x=0, y=0)
    
    plt.xlim(b[0], b[-1])
    
    plt.tight_layout(pad=1.08)
    plt.savefig(svg)
    return plt.gcf()


def bxt_custom(pas_out, _igv_plot, strand, input_regions, num_used:int = 0, 
         boxcolor:int = "black",
         svg:str = "boxplot.tmp.svg"):
    if strand == "-":
        readthrough_len_ = int(pas_out) - np.array(_igv_plot.bam_list[num_used]['start'])
    else:
        readthrough_len_ = np.array(_igv_plot.bam_list[num_used]['end']) - int(pas_out)

    fig, ax = plt.subplots(figsize=(10, 1.6))

    flierprops = dict(marker='o', markersize=3, markerfacecolor='k', markeredgecolor='none')
    boxprops = dict(color = boxcolor, linewidth = 2.5, linestyle = '-')  ### "#fcb847"
    
    whiskerprops = dict(color = boxcolor, linewidth = 2.5, linestyle = '-')

    medianprops=dict(color="black", linewidth = 1.5)
    
    capprops = dict(color = boxcolor, linewidth = 2.5, linestyle = '-')
    
    ax.boxplot(readthrough_len_, widths=.5, showfliers = True, 
               whiskerprops = whiskerprops,
               medianprops = medianprops,
               capprops = capprops,
               flierprops=flierprops, boxprops = boxprops, vert=False)
    
    plt.xlabel('Readthrough distance (nt)')

    sns.despine(top=True, right=True, left=True)

    b = input_regions
    # b = np.array([-200, 0, 500, 1000, 1500])
    xticks_str = [str(x) for x in b]
    plt.xticks(b, xticks_str )

    plt.yticks([])

    ax.axvline(0, ls="--", color="#555555")
    # ax.margins(x=0, y=0)
    
    plt.xlim(b[0], b[-1])
    
    plt.tight_layout(pad=1.08)
    plt.savefig(svg)
    return plt.gcf()

