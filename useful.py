import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
from drormd import plot
import itertools
import sys
import re
from collections import OrderedDict 
from scipy.stats import spearmanr
from itertools import combinations
from scipy.stats import ks_2samp

#This is how you would import it:
    # Correctly expand the path to your Desktop/dror directory

#script_dir = os.path.expanduser("~/Desktop/dror")
#sys.path.append(script_dir)

# #Now try to import

#import useful
#importlib.reload(useful)

################ KS STATISTICS ####################

def collect_metric_data_by_group(data, condition_list):
    """Collect and flatten all values for each metric across the specified conditions."""
    metric_data = {}
    for cond in condition_list:
        if cond not in data:
            continue
        for rep_key in data[cond]:
            df = data[cond][rep_key]
            for metric in df.columns:
                if metric not in metric_data:
                    metric_data[metric] = []
                metric_data[metric].extend(df[metric].dropna().tolist())
    return {k: pd.Series(v) for k, v in metric_data.items()}

def compute_ks_by_metric(data, condition_groups, display=True, hline=0.3):
    """
    Computes the KS statistic for each metric between two groups of conditions.

    Args:
        data (dict): nested dictionary of the form data[condition][replicate] = DataFrame
        condition_groups (dict): e.g., {'GroupA': ['cond1', 'cond3'], 'GroupB': ['cond2']}
        display (bool): whether to show a bar plot of KS statistics

    Returns:
        ks_dict (dict): {metric_name: ks_statistic}
    """
    if len(condition_groups) != 2:
        raise ValueError("You must provide exactly two groups for KS comparison.")

    group_names = list(condition_groups.keys())
    group1_name, group2_name = group_names[0], group_names[1]
    group1_data = collect_metric_data_by_group(data, condition_groups[group1_name])
    group2_data = collect_metric_data_by_group(data, condition_groups[group2_name])

    ks_dict = {}

    # Compute KS for all overlapping metrics
    common_metrics = set(group1_data).intersection(group2_data)

    for metric in sorted(common_metrics):
        ks_stat, _ = ks_2samp(group1_data[metric], group2_data[metric])
        ks_dict[metric] = ks_stat

    if display:
        # Sort metrics by KS stat for clearer plotting
        sorted_metrics = sorted(ks_dict.items(), key=lambda x: x[1], reverse=True)
        labels = [k for k, _ in sorted_metrics]
        values = [v for _, v in sorted_metrics]

        plt.figure(figsize=(max(8, len(labels) * 0.6), 5))
        plt.axhline(hline, color='red', linestyle='--')
        plt.bar(labels, values, color='steelblue')
        plt.ylabel("KS Statistic")
        plt.xticks(rotation=45, ha='right')
        plt.title(f"KS Statistic between {group1_name} and {group2_name}")
        plt.tight_layout()
        plt.show()

    return ks_dict



######### CORRELATION ######################



def plot_bar_from_dict(data_dict, title="", ylabel="", xlabel="", sort=True, figsize=(10, 5), color='steelblue'):
    """
    Plots a bar graph from a dictionary of {label: value} pairs.

    Args:
        data_dict (dict): Keys are bar labels, values are bar heights.
        title (str): Title of the plot.
        ylabel (str): Label for the y-axis.
        xlabel (str): Label for the x-axis.
        sort (bool): Whether to sort bars by value descending.
        figsize (tuple): Figure size.
        color (str or list): Bar color(s).
    """
    if sort:
        items = sorted(data_dict.items(), key=lambda x: x[1], reverse=True)
    else:
        items = list(data_dict.items())

    labels, values = zip(*items)

    plt.figure(figsize=figsize)
    plt.bar(labels, values, color=color)
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()




def collect_metric_data(data, conditions, reps):
    """Collects and concatenates metric values across specified conditions and replicates."""
    combined_data = {}

    for cond in conditions:
        for r in range(1, reps+1):
            rep_key = f'{cond}_r{r}'
            if rep_key not in data[cond]:
                continue  # skip missing replicates
            df = data[cond][rep_key]
            for metric in df.columns:
                if metric not in combined_data:
                    combined_data[metric] = []
                combined_data[metric].extend(df[metric].dropna().tolist())

    return {metric: pd.Series(values) for metric, values in combined_data.items()}


def compute_spearman_matrix(metric_series_dict):
    """Compute the Spearman correlation matrix between all pairs of metrics."""
    metrics = list(metric_series_dict.keys())
    corr_matrix = pd.DataFrame(index=metrics, columns=metrics, dtype=float)

    for m1, m2 in combinations(metrics, 2):
        s1, s2 = metric_series_dict[m1], metric_series_dict[m2]
        length = min(len(s1), len(s2))
        corr, _ = spearmanr(s1[:length], s2[:length])
        corr_matrix.loc[m1, m2] = corr
        corr_matrix.loc[m2, m1] = corr

    # Diagonal is 1.0 by definition
    for m in metrics:
        corr_matrix.loc[m, m] = 1.0

    return corr_matrix


def make_correlation_dict(corr_matrix):
    """Convert the correlation matrix to a dictionary of metric pairs -> correlation."""
    corr_dict = {}
    for i in range(len(corr_matrix)):
        for j in range(i+1, len(corr_matrix)):
            m1, m2 = corr_matrix.index[i], corr_matrix.columns[j]
            corr_dict[(m1, m2)] = corr_matrix.loc[m1, m2]
    return corr_dict


def plot_corr_matrix(corr_matrix, figsize=(10, 8), cmap="coolwarm", title="Spearman Correlation Matrix"):
    """Plot a heatmap of the correlation matrix."""
    plt.figure(figsize=figsize)
    sns.heatmap(corr_matrix.astype(float), annot=True, cmap=cmap, fmt=".2f", square=True)
    plt.title(title)
    plt.tight_layout()
    plt.show()


def analyze_metric_correlations(data, conditions, reps, display=True):
    """Main function to collect data, compute correlations, and optionally plot."""
    metric_data = collect_metric_data(data, conditions, reps)
    corr_matrix = compute_spearman_matrix(metric_data)
    if display:
        plot_corr_matrix(corr_matrix)
    return make_correlation_dict(corr_matrix)

def classify_metric(metric, optional=[]):
    return 'side_chain' if 'chi' in metric or metric in optional else 'backbone'

def extract_metric_sets(corr_dict, cutoff, optional=[]):
    categorized = {
        'bb_bb': {},
        'bb_sc': {},
        'sc_sc': {}
    }

    for (m1, m2), corr in corr_dict.items():
        if abs(corr) < cutoff:
            continue
        type1 = classify_metric(m1, optional=optional)
        type2 = classify_metric(m2, optional=optional)

        if type1 == 'backbone' and type2 == 'backbone':
            categorized['bb_bb'][(m1, m2)] = corr
        elif type1 == 'side_chain' and type2 == 'side_chain':
            categorized['sc_sc'][(m1, m2)] = corr
        else:
            # Always store as (backbone, side_chain)
            if type1 == 'backbone':
                categorized['bb_sc'][(m1, m2)] = corr
            else:
                categorized['bb_sc'][(m2, m1)] = corr

    return categorized

def build_square_corr_matrix(pair_dict):
    """Turns a symmetric pairwise correlation dict into a square matrix DataFrame."""
    metrics = sorted({m for pair in pair_dict for m in pair})
    matrix = pd.DataFrame(index=metrics, columns=metrics, dtype=float)
    for m1, m2 in pair_dict:
        val = pair_dict[(m1, m2)]
        matrix.loc[m1, m2] = val
        matrix.loc[m2, m1] = val
    for m in metrics:
        matrix.loc[m, m] = 1.0
    return matrix

def build_rectangular_corr_matrix(pair_dict):
    """Builds a rectangular matrix: rows = backbone, cols = side chain."""
    bb = sorted({m1 for (m1, m2) in pair_dict})
    sc = sorted({m2 for (m1, m2) in pair_dict})
    matrix = pd.DataFrame(index=bb, columns=sc, dtype=float)
    for (m1, m2), val in pair_dict.items():
        matrix.loc[m1, m2] = val
    return matrix

def plot_labeled_heatmap(matrix, title, square=False):
    matrix = matrix.astype(float).copy()
    # Annotation DataFrame of formatted strings
    annot_df = matrix.map(lambda x: f"{x:.2f}" if pd.notna(x) else "")
    print(annot_df.values)
    fig, ax = plt.subplots(figsize=(max(8, len(matrix.columns)*0.6), max(6, len(matrix.index)*0.5)))
    sns.heatmap(
        matrix,
        annot=annot_df.values,        # DataFrame of string labels!
        fmt="s",               # string format
        cmap='coolwarm',
        vmin=-1, vmax=1,
        cbar_kws={'label': 'Spearman Correlation'},
        square=square,
        ax=ax
    )
    ax.set_title(title)
    plt.tight_layout()
    plt.show()



def plot_separated_correlation_matrices(corr_dict, cutoff=0.5, optional=[]):
    """
    Create 3 correlation heatmaps (bb-bb, sc-bb, sc-sc) using a specified cutoff.
    Side chain = contains 'chi'; backbone = everything else.
    """
    categorized = extract_metric_sets(corr_dict, cutoff,optional=optional)

    # Plot BB–BB
    if categorized['bb_bb']:
        mat_bb = build_square_corr_matrix(categorized['bb_bb'])
        plot_labeled_heatmap(mat_bb, f"Backbone vs Backbone (|ρ| > {cutoff})", square=True)
    else:
        print(f"No Backbone–Backbone metric pairs above cutoff {cutoff}")

    # Plot SC–SC
    if categorized['sc_sc']:
        mat_sc = build_square_corr_matrix(categorized['sc_sc'])
        plot_labeled_heatmap(mat_sc, f"Sidechain vs Sidechain (|ρ| > {cutoff})", square=True)
    else:
        print(f"No Side Chain–Side Chain metric pairs above cutoff {cutoff}")

    # Plot BB–SC (rectangular matrix)
    if categorized['bb_sc']:
        mat_mix = build_rectangular_corr_matrix(categorized['bb_sc'])
        plot_labeled_heatmap(mat_mix, f"Backbone vs Sidechain (|ρ| > {cutoff})", square=False)
    else:
        print(f"No Backbone–Side Chain metric pairs above cutoff {cutoff}")







##################### MAKE METRICS FILE ############################

def print_chi1(resnum, last='CG'):
    print(f"""    '{resnum}-chi1': ['{{{resnum}}} and name N',
               '{{{resnum}}} and name CA',
               '{{{resnum}}} and name CB',
               '{{{resnum}}} and name {last}'],""")


def print_chi2(resnum, last='CD1', third='CG'):
    print(f"""    '{resnum}-chi2': ['{{{resnum}}} and name CA',
               '{{{resnum}}} and name CB',
               '{{{resnum}}} and name {third}',
               '{{{resnum}}} and name {last}'],""")


def write_diheds(resnums, chis):
    """
    prints out dihedral dictionary entries for metrics.py
    resnums (list): something like ['A100', 'D401']
    so it contains the atom type and the number
    chis (list): list of chi angles i.e. [1, 2, ...]
    where you want to compute the chi[i] chi angle for
    resnum[i]
    """
    print('dihedrals : {')
    for i,resnum in enumerate(resnums):
        if 'P' in resnum:
            continue
        if chis[i] == 1:
            if 'I' in resnum or 'V' in resnum:
                print_chi1(resnum, last='CG1')
            elif 'T' in resnum or 'D' in resnum:
                print_chi1(resnum, last='OG1')
            elif 'S' in resnum:
                print_chi1(resnum, last='OG')
            else:
                print_chi1(resnum)
        elif chis[i] == 2:
            if 'M' in resnum:
                print_chi2(resnum, last='SD')
            elif 'N' in resnum:
                print_chi2(resnum, last='ND2')
            elif 'I' in resnum:
                print_chi2(resnum, last='CD', third='CG1')
            else:
                print_chi2(resnum)
    print('},')



########################## SIM ANALYSIS #########################


def add_min_column(csv_path, column1, column2, new_column_name, output_path=None):
    """
    Given a CSV file and two column names, this function creates a new column containing
    the minimum value of the two specified columns for each row and saves the modified CSV.

    Parameters:
        csv_path (str): Path to the input CSV file.
        column1 (str): Name of the first column.
        column2 (str): Name of the second column.
        new_column_name (str): Name of the new column that will store the minimum values.
        output_path (str, optional): Path to save the modified CSV. If None, overwrites the input file.

    Returns:
        pd.DataFrame: The modified DataFrame with the new column.
    """
    # Load CSV
    df = pd.read_csv(csv_path)

    # Ensure columns exist
    if column1 not in df.columns or column2 not in df.columns:
        raise ValueError(f"Columns '{column1}' and/or '{column2}' not found in CSV.")

    # Compute minimum value row-wise
    df[new_column_name] = df[[column1, column2]].min(axis=1)

    # Save the modified CSV
    output_path = output_path or csv_path  # If no output path is given, overwrite the input file
    df.to_csv(output_path, index=False)

    return df


def make_lists(data, cond, reps, start=0, finish=1000, offset=0, skip=False, met=None, filt=False):
    """
    Extracts a combined list of metric values from multiple simulation replicates for a given condition.

    Parameters:
        data (dict): Nested dictionary containing simulation data organized as data[condition][replicate] = DataFrame.
        cond (str): Condition name (e.g., 'WT', 'mutantA') used to access data[cond].
        reps (int): Number of replicates for the condition.
        start (int): Start index for slicing each metric trace. Default is 0.
        finish (int): End index for slicing each metric trace. Default is 1000.
        offset (int): Offset added to start/finish indices. Useful for burn-in trimming.
        skip (bool): Whether to skip replicate 10 (useful for known problematic runs). Default is False.
        met (str): The metric to extract (e.g., 'chi1', 'distance', etc.). If None, nothing is extracted.
        filt (tuple or bool): Optional filtering condition as a tuple: (column, operator, threshold),
                              e.g., ('chi1', '<', 90). If False, no filtering is applied.

    Returns:
        list: Flattened list of values across all replicates for the specified metric, optionally filtered and windowed.
              Dihedral angles are wrapped to the [0, 360) or [-180, 180] range depending on content.
    """
    if met:
        met_dist = []
        for i in range(1, reps + 1):
            if i == 10 and skip == True:
                continue
            r = f'{cond}_r{i}'
            data_r = data[cond][r]
            if filt:
                if filt[1] == '<':
                    data_r_2 = data_r[data_r[filt[0]] < filt[2]].reset_index(drop=True)
                elif filt[1] == '>':
                    data_r_2 = data_r[data_r[filt[0]] > filt[2]].reset_index(drop=True)
            if 'chi' in met:
                if not filt:
                    for i in range(start + offset, finish + offset):
                        if i < len(data_r[met]):
                            val = data_r[met][i]
                            if val < 0:
                                val += 360
                            if val > 360:
                                val -= 360
                            met_dist.append(val)
                else:
                    for i in range(len(data_r_2[met])):
                        if i < len(data_r_2[met]):
                            val = data_r_2[met][i]
                            if val < 0:
                                val += 360
                            if val > 360:
                                val -= 360
                            met_dist.append(val)
            else:
                if not filt:
                    met_dist += list(data_r[met].values)[start + offset:finish + offset]
                else:
                    met_dist += list(data_r_2[met].values)
        return met_dist




def sliding_mean(data_array, window):
    new_list = []
    for i in range(data_array.shape[0]):
        if i <= window:
            indices = range(0,2*i+1)
        else:
            indices = range(max(i - window, 0),
                            min(i + window + 1, len(data_array)))
        avg = 0
        for j in indices:
            avg += data_array[j]
        avg /= float(len(indices))
        new_list.append(avg)
    return np.array(new_list)


colors_lst = ['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']




def get_values(data, condition, nreps, metric):
    reps = ['r'+str(i) for i in range(1,nreps+1)] 
    traces = []
    for rep in reps:
        trace = data[condition][condition+'_'+rep][metric].values
        if 'chi' in metric:
            print(type(trace))
            trace_new = []
            for val in trace:
                if val < 0:
                    val += 360
                if val > 230:
                    val -= 360
                trace_new.append(val)
            traces.append(np.array(trace_new))
        else:
            traces.append(trace)
    return traces


####################### PLOTTING ###########################

def plot_time_traces(conditions, data, colors=colors_lst, leg_title='', only_mets=[], skip_mets=[],bw=1):
    """
    Plots smoothed time traces and kernel density estimates (KDEs) for a given metric across multiple conditions.

    Parameters:
        conditions (list of str): List of condition names (e.g., ['WT', 'mutantA']) corresponding to keys in `data`.
        data (dict): Nested dictionary organized as data[condition][replicate][metric] = time series array.
        colors (list of str): List of colors to use for each condition. Defaults to a predefined color palette.
        skip_mets (list of str): list of metrics not to plot
        only_mets (list of str): the only metrics to be plotted
        leg_title (str): Title for the legend box.

    Returns:
        None. Generates matplotlib plots for each metric across the specified conditions.

    Behavior:
        - For each metric in the first replicate of the first condition, the function:
            - Gathers all corresponding metric traces across replicates and conditions.
            - Applies special angle wrapping for metrics containing 'chi' to ensure cyclic continuity.
            - Plots time traces and KDE distributions using helper functions from `drormd.plot`.
        - Distinguishes metrics as either dihedral angles or distance-like values to label axes accordingly.
    """
    for metric in data[conditions[0]][f'{conditions[0]}_r1']:
        if only_mets:
            if metric not in only_mets:
                continue
        if skip_mets:
            if metric in skip_mets:
                continue
        trace, kde = plot.setup_trace_and_kde(metric, '')
        for i, condition in enumerate(conditions):
            A = [
            data[condition][rep][metric].dropna().to_numpy()       # <-- NaNs removed
            for rep in data[condition]
            if metric in data[condition][rep]                      # column exists
           and data[condition][rep][metric].notna().any()      # at least one finite value
            ]
            if not A:
                continue
            if 'chi' in metric:
                plot.dihedral_range(A, 0)
                ylabel = 'degree'
            else:
                ylabel = 'distance'
            plot.add_time_trace(trace, A, f'{condition}', colors=colors[i], lw=0, lw_smooth=1)
            plot.add_kde(kde, A, f'{condition}', flip=True, color=colors[i], burnin=0,bandwidth=bw)
        trace.set_xlim(0)
        trace.set_ylabel(ylabel)
        legend = trace.legend(title=leg_title, loc='best')
        plt.show()




def plot_traces_sing(data, conditions, reps, metric, ylim=[],
                    yticks=[], xlim=[], hlines=[], xticks=[]):
    """
    Plot raw and smoothed time-series traces for a given metric across multiple
    replicates of one or more experimental conditions.

    For each element in *conditions*, the function creates a single figure and
    places one subplot per replicate (indexed from 1 to ``reps − 1``) in a
    vertical stack.  Within every subplot the raw trace is shown in light grey
    and a 30-frame sliding-window mean (produced by ``sliding_mean``) is
    over-plotted in dark grey.

    Optional axis-formatting arguments (*ylim*, *yticks*, *xlim*, *hlines*,
    *xticks*) are applied uniformly to every subplot.

    Parameters
    ----------
    data : dict
        Nested dictionary such that
        ``data[condition][replicate][metric]`` returns a 1-D NumPy array
        containing the trace for *metric*.
    conditions : list[str]
        Names of the experimental conditions to visualise (must be keys in
        *data*).
    reps : int
        Number of replicate trajectories stored per condition, assumed to be
        one-based (e.g. ``…_r1``, ``…_r2``, …).
    metric : str
        Metric key to extract and plot.
    ylim : list[float, float], optional
        ``[ymin, ymax]`` limits for all y-axes.
    yticks : list[float], optional
        Custom y-tick locations.
    xlim : list[int, int], optional
        ``[xmin, xmax]`` limits for the x-axes.
    hlines : list[float], optional
        y-values at which to draw horizontal reference lines (dotted, black).
    xticks : list[int], optional
        Custom x-tick locations applied to the bottom subplot only.

    Returns
    -------
    None
        The function draws the plots directly via Matplotlib and does not
        return a value.
    """
    for condition in conditions:
        print(condition)
        traces = get_values(data, condition, reps, metric)
        plt.figure(1, figsize=(5,12))
        for i in range(1,reps):
            plt.subplot(6,1,i)
            plt.plot(sliding_mean(traces[i-1],30), color='gray',linewidth=3)
            plt.plot(traces[i-1], alpha=0.5, color='gray', linewidth=1)
            if ylim:
                plt.ylim(ylim[0], ylim[1])
            if yticks:
                plt.yticks(yticks,size=20)
            if xlim:
                plt.xlim(xlim[0], xlim[1])
            if i != reps:
                plt.xticks([])
            else:
                plt.xticks(size=20)

            if hlines:
                for hline in hlines:
                    plt.axhline(hline,linestyle=':', color='k')
            if xticks:
                plt.xticks(xticks)
        plt.show()  


def plot_smooth_frequency(data, labels=None, title="Smoothed Frequency Distribution",
                          vert=None, leg=True, xaxis='Value', colors=None, xlim=None,xticks=[],
                          alphas=[1]*6, dashed=[], thin=[]):
    """
    Plots smoothed frequency distributions (Kernel Density Estimates) for multiple datasets.

    Parameters:
        data (list of array-like): List of numerical datasets to visualize via KDE.
        labels (list of str, optional): Names for each dataset, used in the legend. Defaults to 'Dataset 1', etc.
        title (str): Title of the plot.
        vert (list of float, optional): List of x-values at which to draw vertical dashed reference lines.
        leg (bool): If True, includes a legend on the plot.
        xaxis (str): Label for the x-axis.
        colors (list of str, optional): Colors for each dataset line. Defaults to a preset palette if not given.
        xlim (list of float, optional): Specifies the x-axis limits as [xmin, xmax].

    Returns:
        None. Displays a matplotlib plot with smoothed frequency distributions.
    """
    if labels is None:
        labels = [f"Dataset {i+1}" for i in range(len(data))]

    if colors is None or len(colors) < len(data):
        default_palette = ['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666']
        colors = (colors or []) + default_palette[len(colors or []):len(data)]

    # Compute KDEs
    kdes = [gaussian_kde(d) for d in data]

    # Determine x range
    all_min = min(np.min(d) for d in data)
    all_max = max(np.max(d) for d in data)
    dist = all_max - all_min
    x_values = np.linspace(all_min - dist * 0.05, all_max + dist * 0.05, 1000)

        
    # Plot
    plt.figure(figsize=(10, 5))
    for i, kde in enumerate(kdes):
        y = kde(x_values)
        if i in dashed:
            linestyle = '--'
        else:
            linestyle = '-'
        if i in thin:
            linewidth = 3
        else:
            linewidth = 6

        plt.plot(x_values, y, label=labels[i], linewidth=linewidth, color=colors[i], alpha=alphas[i],
                 linestyle=linestyle)

    if vert:
        for val in vert:
            plt.axvline(val, linestyle='dashed', color='red')

    plt.xlabel(xaxis, size=25)
    plt.ylabel("Frequency", size=30)
    if xticks:
        plt.xticks(xticks)
    plt.xticks(size=30)
    plt.yticks([])
    plt.title(title)
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)

    if xlim:
        plt.xlim(xlim[0], xlim[1])

    if leg:
        plt.legend()

    plt.show()

###################### MAKE CONDITIONS FROM METRICS ########################


restypes = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
            'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

nums = ['1', '2', '3', '4', '5', '6', '7', '8', '9']

metrics_dir = sys.argv[1]

def extract_unique_keys(nested, preserve_order=True):
    """Return all substrings inside {...} anywhere in *nested*."""
    
    def walk(obj):
        """Yield every string found at any depth inside obj."""
        if isinstance(obj, str):
            yield obj
        elif isinstance(obj, dict):
            for v in obj.values():
                yield from walk(v)
        elif isinstance(obj, (list, tuple, set)):
            for item in obj:
                yield from walk(item)
        # ignore non-iterable, non-string leaf nodes (ints, floats, etc.)

    # Collect every string in one pass, then run the regex once
    all_text = " ".join(walk(nested))
    hits = re.findall(r"\{([^}]*)\}", all_text)

    if preserve_order:
        return list(OrderedDict.fromkeys(hits))  # de-duplicate, keep order
    else:
        return sorted(set(hits))     # Fast, but order not guaranteed


def get_helices(keys, TM=False):
    helices_used = set()
    
    for key in keys:
        if 'x' in key:
            helix = key.split('x')[0]
            if helix.isdigit():
                    helices_used.add(int(helix))
    
        elif 'TM' == key[:2] or 'H' == key[0]:
            if key[0] != 'H':
                try:
                    helices_used.add(int(key[2:]))
                except ValueError:
                    continue
            elif key[0] == 'H':
                helices_used.add(int(key[1:]))

    return helices_used


def get_ref_numbers(helices_used):
    # Prompt for segid
    segid = input("\n Enter the segid for the protein: ").strip()

    # Prompt for reference residue numbers (.50) for each helix
    helix_refs = {}
    print('\n')
    for helix in sorted(helices_used):
        ref_res = int(input(f"Enter the residue number for {helix}.50: "))
        helix_refs[helix] = {'ref_res': ref_res, 'segid': segid}

    return helix_refs, segid


def process_0x00(keys, helix_refs, segid=False, ref=False):
    dict_0x00 = {}
    for key in keys:
        helix, gpcrdb_pos = key.split('x')
        if ref:
            gpcrdb_pos = gpcrdb_pos[:-4]
        if helix.isdigit():
            helix_num = int(helix)
            gpcrdb_pos = int(gpcrdb_pos)
            offset = gpcrdb_pos - 50
            ref_res = helix_refs[helix_num]['ref_res']
            resid = ref_res + offset

            if segid:
                dict_0x00[key] = f"protein and segid {segid} and resid {resid}"
            else:
                dict_0x00[key] = f"protein and resid {resid}"
    return dict_0x00


def process_TMX(keys, helix_refs, segid=False, ref=False):
    dict_TMX = {}

    for key in keys:
        if key[0] == 'H':
            helix_num = int(key[1])
        else:
            try:
                helix_num = int(key[2])
            except ValueError:
                print(f"Can't process {key}, figure out on your own")
                continue

        ref_res = helix_refs[helix_num]['ref_res']
        segid = helix_refs[helix_num]['segid']

        # Prompt for GPCRDB range
        gpcrdb_range = input(f"Enter the GPCRDB range for {key} (e.g., 3.44-3.56): ").strip()
        start_gpcrdb, end_gpcrdb = gpcrdb_range.split('-')
        start_gpcrdb = int(start_gpcrdb.split('.')[1])
        end_gpcrdb = int(end_gpcrdb.split('.')[1])

        # Calculate residue range
        start_res = ref_res + (start_gpcrdb - 50)
        end_res = ref_res + (end_gpcrdb - 50)
        if ref:
            key += '_ref'
        if segid:
            dict_TMX[key] = f"protein and segid {segid} and resid {start_res} to {end_res}"
        else:
            dict_TMX[key] = f"protein and resid {start_res} to {end_res}"

    return dict_TMX



def process_X000(keys, segid=False, ask_segid=True):
    dict_X000 = {}
    if ask_segid:
        segid = input("\n Enter the segid for the protein: ").strip()

    if not keys:
        return {}
    for key in keys:
        num = key[1:]
        if segid:
            dict_X000[key] = f"protein and segid {segid} and resid {num}"
        else:
            dict_X000[key] = f"protein and resid {num}"
    return dict_X000

                
def process_BW(keys_0x00, keys_TMX=[], ref=False):
    # Identify the helices present
    if not keys_0x00:
        return {}
    if ref:
        print('Enter the following for the reference system')
        
    helices_used = get_helices(keys_0x00 + keys_TMX)

    helix_refs, segid = get_ref_numbers(helices_used)
    
    res_dict = {}
    
    res_dict = {**res_dict, **process_0x00(keys_0x00, helix_refs, segid, ref=ref)}
    
    if keys_TMX:
        res_dict = {**res_dict, **process_TMX(keys_TMX, helix_refs, segid, ref=ref)}

    return res_dict, segid


def make_cond_dict(metric_dict):
    keys = extract_unique_keys(metric_dict)
    print(keys)

    keys_0x00 = []
    keys_X000 = []
    keys_TMX = []
    keys_ref = []
    unknown = []
    for key in keys:
        if 'x' in key and 'ref' not in key:
            keys_0x00.append(key)
        elif 'H8' in key:
            keys_TMX.append(key)
        elif 'ref' in key and 'selection' not in key:
            keys_ref.append(key)
        elif key[0] in restypes and key[1] in nums:
            keys_X000.append(key)
        elif 'TM' == key[:2] and key[2] in nums:
            keys_TMX.append(key)
        else:
            print(f"Couldn't process key {key}, add on your own")
            unknown.append(key)

    ask_segid = True
    if keys_0x00 or keys_TMX:
        cond_dict_BW, segid = process_BW(keys_0x00, keys_TMX=keys_TMX)
        ask_segid = False
    else:
        cond_dict_BW = {}
    if keys_X000:
        if not ask_segid:
            cond_dict_X000 = process_X000(keys_X000, segid=segid)
        else:
            cond_dict_X000 = process_X000(keys_X000, ask_segid=ask_segid)
    else:
        cond_dict_X000 = {}

    cond_dict_ref = {}
    ref_same = ''

    if keys_ref:
        ref_same = input(f"\nShould the reference selections be the same as the original? Y/N : ").strip()
    
    if ref_same == 'Y':
        for key in keys_ref:
            cond_dict_ref[key] = cond_dict_BW[key[:-4]]
    elif ref_same == 'N':
        cond_dict_ref, _ = process_BW(keys_ref, ref=True)
    
    total_dict = {**cond_dict_BW, **cond_dict_X000, **cond_dict_ref}
    
    for key in unknown:
        total_dict[key] = ''

    return total_dict
