import matplotlib as mpl
from matplotlib.backends.backend_pgf import FigureCanvasPgf
import matplotlib.pyplot as plt
import matplotlib as mpl
from cycler import cycler
import seaborn as sns
import numpy as np
import pandas as pd

# -- Colours
def adjust_lightness(color, amount=1):
    # amount > 1 -> lighten
    # amount < 1 -> darken
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])

def create_my_palette(lightness=1):
    my_palette = [(55, 125, 185), (230, 111, 0), (75, 175, 75), (150, 80, 165), (230, 219, 90), (222, 24, 28)]
    # transform to 0-1 instead of 0-255
    for color_idx, color in enumerate(my_palette):
        rgb = [0, 0, 0]
        for component_idx, component in enumerate(color):
            rgb[component_idx] = component / 255
        my_palette[color_idx] = adjust_lightness((rgb[0], rgb[1], rgb[2]), lightness)
    return my_palette  

# -- Figure and ticks
def figure_dimensions(width_ratio=1, aspect_ratio=0.61803):
    # Computes the size of a figure in inches from the relative size on a paper
    # Specific to the dimensions of the publication latex template.
    page_width = 6.14 #inch
    return [width_ratio * page_width, width_ratio*aspect_ratio*page_width]

def scale_ticks(ticks, axis_range):
    # Scales ticks for imagesc plots
    ticks = np.array(ticks)
    base_ticks = (ticks - min(ticks)) / (max(ticks) - min(ticks)) # between 0 and 1
    axis_ticks = base_ticks * (axis_range[1] - axis_range[0]) + axis_range[0]
    return axis_ticks

# -- Matplotlib settings
def configure_matplotlib():

    mpl.style.use('ggplot')
    mpl.rcParams['axes.prop_cycle'] = cycler(color=create_my_palette())
    mpl.rcParams['image.cmap'] = 'magma_r'
    mpl.rcParams['figure.figsize'] = figure_dimensions()
    mpl.rcParams['axes.facecolor'] = adjust_lightness('#E5E5E5', 1.065)
    mpl.rcParams['text.color'] = 'black'
    mpl.rcParams['axes.labelcolor'] = 'black'
    mpl.rcParams['xtick.color'] = 'black'
    mpl.rcParams['ytick.color'] = 'black'

    # -- font settings for latex figures
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['font.serif'] = ["TeX Gyre Pagella"]    # may require manual installation
    mpl.rcParams['font.size'] = 10                       # LaTeX default is 10pt font.
    mpl.rcParams['axes.labelsize'] = 10          
    mpl.rcParams['legend.fontsize'] = 8
    mpl.rcParams['xtick.labelsize'] = 8
    mpl.rcParams['ytick.labelsize'] = 8

    # -- save to pgf settings
    mpl.backend_bases.register_backend('pdf', FigureCanvasPgf)
    mpl.rcParams['savefig.format'] = 'pdf'
    mpl.rcParams['savefig.dpi'] = 300
    mpl.rcParams['axes.unicode_minus'] = True
    mpl.rcParams['text.usetex'] = False
    mpl.rcParams['pgf.texsystem'] = 'lualatex'     # Only applied at saving
    mpl.rcParams['pgf.preamble'] = r'\usepackage[T1]{fontenc}'+ \
        r'\usepackage{fontspec}'+ \
        r'\usepackage[utf8]{inputenc}'+ \
        r'\usepackage{unicode-math}'+ \
        r'\usepackage{siunitx}'+ \
        r'\usepackage{amsmath}'+ \
        r'\usepackage{bm}'
    
# -- Plotting (used by SNR)
def plot_spatial_heatmap(x_col, y_col, value_col, value_levels,
                         aggfunc=np.median,
                         xticks=[-0.5, -0.2, 0, 0.2, 0.5],
                         yticks=[0, 0.1, 0.3, 0.5],
                         cmap='viridis',
                         **kwargs):
    data = kwargs.pop('data')
    if not data.empty:
        d = data.pivot_table(index=y_col, columns=x_col, values=value_col, aggfunc=aggfunc)
        plt.contourf(d, levels=value_levels, vmin=min(value_levels), vmax=max(value_levels), cmap=cmap)

    ax = plt.gca()
    ax.set_aspect('equal','box')
    ax.set_xticks(scale_ticks(xticks, ax.get_xlim()))
    ax.set_xticklabels([f'{l}'.replace('-', '\N{MINUS SIGN}') for l in xticks])
    ax.set_yticks(scale_ticks(yticks, ax.get_ylim()))
    ax.set_yticklabels([f'{l}'.replace('-', '\N{MINUS SIGN}') for l in yticks], rotation=0)
    ax.set_xlabel(r'\(x\) (\si{m})')
    ax.set_ylabel(r'\(y\) (\si{m})')
    
def plot_orienation_heatmap(x_col, y_col, value_col, value_levels,
                         aggfunc=np.median,
                         xticks=[0, 0.25, 0.75, 1, 1.25, 1.75],
                         xticklabels=[r'\(0\pi\)', r'\(\frac{1}{4}\pi\)', r'\(\frac{3}{4}\pi\)', 
                                      r'\(1\pi\)', r'\(1\frac{1}{4}\pi\)', r'\(1\frac{3}{4}\pi\)'],
                         yticks=[0.1, 0.3, 0.5],
                         yticklabels = None,
                         cmap='viridis',
                         labels=True,
                         **kwargs):
    data = kwargs.pop('data')
    if not data.empty:
        d = data.pivot_table(index=y_col, columns=x_col, values=value_col, aggfunc=aggfunc)
        r = np.sort(data[y_col].unique())[np.newaxis]
        theta = np.sort(data[x_col].unique())[np.newaxis]
        r, theta = np.meshgrid(theta, r)    
        plt.contourf(r, theta, d, levels=value_levels, vmin=min(value_levels), vmax=max(value_levels), cmap=cmap)
        
    ax = plt.gca()
    ax.set_aspect('equal','box')
    ax.set_yticks([])
    ax.set_xticks([])
    # Overlay an axes for the ticks
    rect = ax.get_position()
    ax2 = plt.gcf().add_axes(rect, polar=True)
    ax2.set_xlim(ax.get_xlim())
    ax2.set_ylim(ax.get_ylim())
    
    ax2.set_yticks(yticks)
    if yticklabels is None:
        ax2.set_yticklabels([str(t) for t in yticks], 
                            color='white',
                            fontsize=mpl.rcParams['ytick.labelsize'])
    else:
        ax2.set_yticklabels(yticklabels, 
                            color='white',
                            fontsize=mpl.rcParams['ytick.labelsize'])
    ax2.set_xticks(np.array(xticks) * np.pi)
    ax2.set_xticklabels(xticklabels)
    ax2.tick_params(axis='x', which='major', pad=-3)
    ax2.tick_params(axis='both', grid_linewidth = 0.2)
    if labels:
        ax2.set_xlabel(r'\(\varphi\) (\si{rad})')
        ax2.set_ylabel(r'\(y\) (\si{m})', labelpad=15)
    ax2.set_facecolor([0, 0, 0, 0])
    
#     bbox = dict(boxstyle="round", ec="white", fc="white",
#                 alpha=0.5, pad = 0.1)
#     plt.setp(ax2.get_yticklabels(), bbox=bbox)
    
# makes subplot adjust apply to the grid-axes of the polar plots
def fix_polar_grid(fig, is_grid=False):
    if not is_grid:
        mask = [isinstance(a, mpl.projections.polar.PolarAxes) for a in fig.axes]
        indices = np.nonzero(mask)[0].tolist()
        # skip elke tweede
        for idx in indices[1::2]:
            fig.axes[idx].set_position(fig.axes[idx-1].get_position())
        return fig
    
    offset = int(len(fig.axes) / 2);
    for idx in range(offset, len(fig.axes)):
        fig.axes[idx].set_position(fig.axes[idx-offset].get_position())
    return fig
    
    
# -- Bar plot functions
def contourColors(levels, cmap=None):
    """Returns the colors used in a contour plot.
    
    mpl.contourf uses the middle of each level to compute the colour.
    
    Parameters
    ----------
    levels: list
        Values of the contour levels
    cmap: optional,
        Name of the colormap to use
    """
    if cmap is None:
        cmap = mpl.cm.get_cmap(mpl.rcParams['image.cmap'])
    else:
        cmap = mpl.cm.get_cmap(cmap)
    if levels[0] != 0:
        levels = [0] + levels
    levels = np.array(levels)
    normalizer = mpl.colors.Normalize(vmin=min(levels), vmax=max(levels))
    color_values = levels[0:-1] + (levels[1::] - levels[0:-1]) / 2
    color_values = np.r_[color_values, max(levels)]
    return cmap(normalizer(color_values))

def groupedStackedBars(data, group, bar, stack, value,
                       group_order=None, bar_order=None, stack_order=None,
                       group_labels=None, bar_labels=None, stack_labels=None,
                       stack_title=None, colors=None, edge_color='white',
                       group_width_factor=0.8, bar_width_factor=0.8,
                       group_label_offset=-30,
                       bar_label_rotation=0):
    """Plots groups of stacked bars.
    
    Parameters
    ----------
    data: Pandas dataframe
        
    group: string,
        Name of column to use as group
    bar: string,
        Name of the column to use as bar in each group
    stack: string,
        Name of the column to use as stack in each bar
    value: string,
        Name of the column to use as value for each stack
    colors: list, optional
        Colors to use for each stack. If empty, default colorcycle is used.
    group_width_factor: float
        Factor of group width. When 1, groups are touching.
    bar_width_factor: float
        Factor of bar width. When 1, bars are touching
    group_label_offset: float
        y-locations offset of group labels in points. Depends on font size and padding.
    bar_label_rotation: float
        Rotation of the bar labels.
    """
    
    # Fill in None values
    if group_order is None:
        group_order = data[group].unique()
    if bar_order is None:
        bar_order = data[bar].unique()
    if stack_order is None:
        stack_order = data[stack].unique()
    if bar_labels is None:
        bar_labels = bar_order
    if group_labels is None:
        group_labels = group_order
    if stack_labels is None:
        stack_labels = stack_order
    if stack_title is None:
        stack_title = stack
    if colors is None:
        colors = mpl.cm.get_cmap(mpl.rcParams['image.cmap'])
    
    # Compute bar and group positions
    num_groups = len(group_order)
    num_bars = len(bar_order)
    num_stacks = len(stack_order)
    bar_width = group_width_factor / num_bars;
    group_offset = np.arange(num_groups) + 1
    bar_positions = np.array([group_offset + bar_idx * bar_width for bar_idx in range(num_bars)])
    group_centers = np.array([group_offset[group_idx] + (num_bars - 1) / 2  * bar_width 
                              for group_idx in range(num_groups)])
    empty_bars_left = 0
    
    # Plot the stacks in the given order
    for group_idx in range(num_groups):
        group_data = data[data[group] == group_order[group_idx]]
        for bar_idx in range(num_bars):
            bar_data = group_data[group_data[bar] == bar_order[bar_idx]]
            value_so_far = 0
            for stack_idx in range(num_stacks):
                stack_data = bar_data[bar_data[stack] == stack_order[stack_idx]]
                if not stack_data.empty:
                    plt.bar(bar_positions[bar_idx][group_idx], 
                        stack_data[value],
                        bottom=value_so_far,
                        width=bar_width_factor * bar_width,
                        edgecolor=edge_color,
                        color=colors[stack_idx % len(colors)])
                    value_so_far = value_so_far + stack_data[value].tolist()[0]
            if group_idx == 0 and value_so_far == 0:
                empty_bars_left += 1
        
        # Add group labels
        ax=plt.gca()
        
        # We are removing empty bars on the left
        if group_idx == 0:
            new_offset = group_offset[group_idx] + bar_width * empty_bars_left
            new_center = new_offset + (num_bars - 1 - empty_bars_left) / 2 * bar_width
            group_centers[group_idx] = new_center
            
        plt.annotate(group_labels[group_idx], 
                     xy=(group_centers[group_idx], 0),
                     xycoords=ax.get_xaxis_transform(),
                     xytext=(0, group_label_offset),
                     textcoords='offset points',
                     ha="center", va="center",
                     fontsize=mpl.rcParams['xtick.labelsize'])        
        
    # Add ticks, labels, legend...
    tick_positions = bar_positions.T.flatten()[empty_bars_left:]
    tick_labels =  np.tile(bar_labels, num_groups)[empty_bars_left:]
    plt.xticks(tick_positions, tick_labels, rotation=bar_label_rotation)
    plt.legend(stack_labels, title=stack_title, loc='upper right')