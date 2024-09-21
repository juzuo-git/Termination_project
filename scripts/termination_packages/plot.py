#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author       : windz
FilePath     : plot.py
'''

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import Normalize
from scipy.interpolate import interpn
from matplotlib.patches import PathPatch


def color_pal(key):
    '''
    Custom color palette
    
    To preview: 
        for pal in color_palette.values():
            sns.palplot(pal) 
    '''

    color_palette = {
        'set1': ['#6e9ece', '#e6928f', '#4e9595', '#84574d', '#8d6ab8', '#efdbb9', '#76ba80'], 
        'set2': ['#46a1cd', '#ce3c35', '#4258a1', '#57b058', '#7b4b99', '#f2df52', '#a9a9a9'],
        'default': ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    }

    return color_palette[key]


def despine(ax, top=True, right=True, left=False, bottom=False, xaxis=False, yaxis=False):
    if top:
        ax.spines['top'].set_visible(False)
    if right:
        ax.spines['right'].set_visible(False)
    if left:
        ax.spines['left'].set_visible(False)
    if bottom:
        ax.spines['bottom'].set_visible(False)

    if xaxis:
        ax.get_xaxis().set_visible(False)
    if yaxis:
        ax.get_yaxis().set_visible(False)



def plot_scatter_plot(
    data1: pd.DataFrame,
    data2: pd.DataFrame,
    on: str = 'gene_id',
    key: str = 'data',
    data1_label: str = None,
    data2_label: str = None,
    data1_font_style: str = None,
    data2_font_style: str = None,
    s: int = 5,
    p_values: float = .001,
    figsize: tuple = (3, 3),
    xlim: tuple = (0, 1000),
    ylim: tuple = (0, 1000)
):
    '''
    Plot scatterplot with p_values
    Args:
        data1, data2: input DataFrame with raw data
        on: Column or index level names to join on
        key: Column names of the raw data
        s: dot size
        p_values: threshold of the significant

    Return:
        data: merge Dataframe
        ax: pre-existing axes for the plot
    '''
    from scipy.stats import mannwhitneyu

    data = pd.merge(data1, data2, on=on)
    data.dropna(inplace=True)

    data['p_values'] = data.apply(lambda x: mannwhitneyu(
        x[f'{key}_x'], x[f'{key}_y'])[1], axis=1)
    data['data_x'] = data[f'{key}_x'].map(np.median)
    data['data_y'] = data[f'{key}_y'].map(np.median)

    data['color'] = data['p_values'].map(lambda x: 1 if x < p_values else 0)

    plt.figure(figsize=figsize)
    plt.scatter(
        x='data_x',
        y='data_y',
        data=data.query('p_values > @p_values'),
        s=s,
        c='#D2D2D2',
        alpha=.5
    )

    plt.scatter(
        x='data_x',
        y='data_y',
        data=data.query('p_values < @p_values and data_y > data_x'),
        s=s,
        c='#3182BD',
        alpha=.5
    )

    plt.scatter(
        x='data_x',
        y='data_y',
        data=data.query('p_values < @p_values and data_y < data_x'),
        s=s,
        c='#C00000',
        alpha=.5
    )

    plt.plot(xlim, ylim, ls='--', color='#555555')
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xlabel(data1_label, style=data1_font_style)
    plt.ylabel(data2_label, style=data2_font_style)
    plt.xticks(np.linspace(xlim[0], xlim[1], 5, dtype=int))
    plt.yticks(np.linspace(ylim[0], ylim[1], 5, dtype=int))

    ax = plt.gca()
    ax.patch.set_visible(False)

    sns.despine(top=True, right=True)

    return data, ax


def boxplot_with_jitter(
    data: list, 
    figsize: tuple = (3, 3), 
    widths = .5, 
    subsample = 1,
    ax = None,
    labels=None
):
    '''
    带有扰动的boxplot
    Args: 
        data: boxplot data, eg. [[...], [...], ...]

    Return:
        ax: pre-existing axes for the plot
    '''
    from random import sample

    if subsample > 1:
        raise ValueError(f'subsample={subsample}, must <= 1')

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    ax.boxplot(
        data,
        labels=labels,
        widths=widths,
        showfliers=False
    )

    ylim = ax.get_ylim()
    for i in range(len(data)):
        y = sample(list(data[i]), int(subsample*len(data[i])))
        x = np.random.normal(1+i, 0.04, size=len(y))
        ax.scatter(x, y, s=1, color='grey', alpha=.5)

    ylim = ax.set_ylim(ylim)

    despine(ax, top=True, right=True)

    return ax


def step_histplot(
    *data: list,
    bins = 'auto',
    labels = None,
    stat = 'density',
    xlabel = None,
    figsize=(4, 3),
    bbox_to_anchor: tuple = (1, 1),
    ax = None,
):
    '''This function takes a list of data and plots a histogram of the data.
    
    Parameters
    ----------
    data : list
        list of data to be plotted
    bins, optional
        number of bins to use
    labels
        list of strings, optional
    stat, optional
        'density' or 'count'
    xlabel
        str
    figsize
        tuple = (width, height)
    bbox_to_anchor : tuple
        tuple = (1, 1)
    
    '''
    if labels is None:
        labels = [None]*len(data)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    for _data, label in zip(data, labels): 
        sns.histplot(
            _data,
            fill=False,
            stat=stat,
            element='step',
            bins=bins,
            label=label
        )

    ax.legend(frameon=False, bbox_to_anchor=bbox_to_anchor)

    plt.xlabel(xlabel)
    sns.despine(top=True, right=True)

    return ax


def read_count_per_gene(
    gene_counts: dict,
    bs = .2,
    max_counts = 3,
    all_gene_counts = None,
    ylabel = 'Gene counts',
    ):
    '''This function plots the number of reads per gene for a given sample
    
    Parameters
    ----------
    gene_counts : dict
        a dictionary of gene counts
    bs
        the bin size for the histogram
    max_counts, optional
        the maximum number of counts to plot.
    all_gene_counts
        a list of all the counts for all genes.
    ylabel, optional
        the label for the y-axis
    
    '''
    if all_gene_counts is None:
        # 如果没有给出基因总数，则默认注释文件里面的基因总数
        from utils.get_overlap_genes import _, total_gene_counts
        all_gene_counts = total_gene_counts

    df = pd.DataFrame(gene_counts.items(), columns=['gene_id', 'counts'])
    read_count_per_gene = list(np.log10(df['counts']+1))
    gene_with_null_read = all_gene_counts-len(df)  # 在文库中没有被测到的基因总数
    read_count_per_gene.extend([0]*gene_with_null_read)
    read_count_per_gene = np.array(read_count_per_gene)

    counts = read_count_per_gene

    i = 0
    x, y = [], []
    while i <= max_counts:
        mask = (counts <= i) & (counts > i-bs)
        y.append(len(counts[mask]))
        x.append(i)
        i += bs

    plt.figure(figsize=(4, 3))
    plt.bar(x, y, width=bs)
    plt.xlabel('$\log_{10}\mathrm{(read\ counts + 1)}$')
    plt.ylabel(ylabel)
    sns.despine(top=True, right=True)

    ax = plt.gca()

    return ax


def density_scatter(x, y, sort=True, bins=20, figsize=(4, 4), **kwargs):
    '''It takes in two lists of numbers, and returns a scatter plot with the points colored by the density of points in the 2D histogram
    
    Parameters
    ----------
    x
        x-coordinates of the points
    y
        y-coordinates of the points
    sort, optional
        If True, the points are colored by density, so that the densest points are plotted last.
    bins, optional
        The number of bins to use for the histogram.
    figsize
        the size of the figure
    
    '''

    fig, ax = plt.subplots(figsize=figsize)
    data, x_e, y_e = np.histogram2d(x, y, bins=bins, density=True)
    z = interpn((0.5*(x_e[1:] + x_e[:-1]), 0.5*(y_e[1:]+y_e[:-1])),
                data, np.vstack([x, y]).T, method="splinef2d", bounds_error=False)

    # To be sure to plot all data
    z[np.where(np.isnan(z))] = 0.0

    # Sort the points by density, so that the densest points are plotted last
    if sort:
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    ax.scatter(x, y, c=z, **kwargs)

    norm = Normalize(vmin=np.min(z), vmax=np.max(z))

    # The dimensions [left, bottom, width, height] of the new Axes. 
    # All quantities are in fractions of figure width and height.
    cax = fig.add_axes([0.95, 0.1, 0.03, .8])
    
    cbar = fig.colorbar(
        cm.ScalarMappable(norm=norm), 
        cax, orientation='vertical')
    cbar.ax.set_ylabel('Density')

    return ax


def karyotype(
        data: dict, 
        genome_size: dict,
        height: float = .9,
        spacing: float = .9,
        marker_linewidth: float = .5,
        marker_color: str = 'k',
        marker_alpha: float = 1,
        chrom_color: str = '#EEEEEE',
        formatter: str = None,
        annotation_region: dict = None,
        ax = None,
        figsize: tuple =(12, 8)):
    '''This function plots the karyotype of a genome.
    
    Parameters
    ----------
    data : dict
        a dictionary of chromosome names and the position of site, and extend length, eg:
        {'chr1': (1, 1), 'chr2': (100, 200), 'chr3': (50, 1)}
    genome_size : dict
        a dictionary of chromosome names and their sizes
    height : float
        the height of the chromosome
    spacing : float
        the spacing between chromosomes
    marker_linewidth : float
        the width of the line that separates the chromosomes
    marker_color : str, optional
        color of the marker lines
    chrom_color : str, optional
        The color of the chromosome.
    formatter : str
        str = None,
    ax: matplotlib.axes.Axes, optional
        the axis to plot on
    figsize
        tuple, optional
    
    '''

    from utils.annotation import natural_keys

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    yticks, yticklabels, ymin = [], [], 0
    for chrom in sorted(genome_size, key=natural_keys):
        ax.broken_barh(
            [(0, genome_size[chrom])], (ymin, height), 
            facecolors=chrom_color, edgecolor='k', linewidth=1)
        
        if chrom in data:
            if annotation_region is not None:
                xranges = annotation_region[chrom]
                ax.broken_barh(xranges, (ymin, height), 
                               facecolors='k', edgecolor='k', linewidth=marker_linewidth, alpha=.25)

            xranges = data[chrom]
            ax.broken_barh(
                xranges, (ymin, height),
                facecolors=marker_color, edgecolor=marker_color, linewidth=marker_linewidth, alpha=marker_alpha)
            
        yticks_pos = ymin + height/2
        yticks.append(yticks_pos)
        yticklabels.append(chrom)
                        
        ymin += height + spacing

    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)

    if formatter is not None:
        if formatter == 'G':
            fmt = lambda x, pos: '%1.0f' % (x * 1e-9)
        if formatter == 'M':
            fmt = lambda x, pos: '%1.0f' % (x * 1e-6)
        elif formatter == 'K':
            fmt = lambda x, pos: '%1.0f' % (x * 1e-3)
        ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(fmt))
        ax.set_xlabel(f'chromsome size ({formatter}b)')
    else:
        ax.set_xlabel('chromsome size')

    return ax



def adjust_box_widths(g, fac):
    '''Adjust the widths of a seaborn-generated boxplot.
    
    Parameters
    ----------
    g
        the graph object
        eg. g = plt.figure(figsize=(5, 3))
    fac
        the factor by which to adjust the width of the boxes
    
    Description
    -----------
        Adding space between boxes in grouped boxplot
        https://github.com/mwaskom/seaborn/issues/1076
    
    '''
    # iterating through Axes instances
    for ax in g.axes:

        # iterating through axes artists:
        for c in ax.get_children():

            # searching for PathPatches
            if isinstance(c, PathPatch):
                # getting current width of box:
                p = c.get_path()
                verts = p.vertices
                verts_sub = verts[:-1]
                xmin = np.min(verts_sub[:, 0])
                xmax = np.max(verts_sub[:, 0])
                xmid = 0.5*(xmin+xmax)
                xhalf = 0.5*(xmax - xmin)

                # setting new width of box
                xmin_new = xmid-fac*xhalf
                xmax_new = xmid+fac*xhalf
                verts_sub[verts_sub[:, 0] == xmin, 0] = xmin_new
                verts_sub[verts_sub[:, 0] == xmax, 0] = xmax_new

                # setting new width of median line
                for l in ax.lines:
                    if np.all(l.get_xdata() == [xmin, xmax]):
                        l.set_xdata([xmin_new, xmax_new])


def adjust_bar_widths(ax, new_value) :
    '''
    Adjust the widths of a seaborn-generated barplot.
    
    https://stackoverflow.com/questions/34888058/changing-width-of-bars-in-bar-chart-created-using-seaborn-factorplot
    '''
    for patch in ax.patches :
        current_width = patch.get_width()
        diff = current_width - new_value

        # we change the bar width
        patch.set_width(new_value)

        # we recenter the bar
        patch.set_x(patch.get_x() + diff * .5)