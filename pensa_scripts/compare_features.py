import argparse
import numpy as np

from pensa.features import \
    read_structure_features, \
    sort_features, \
    select_common_features, \
    get_multivar_res_timeseries
from pensa.statesinfo import \
    get_discrete_states
from pensa.comparison import \
    relative_entropy_analysis, \
    kolmogorov_smirnov_analysis, \
    ssi_ensemble_analysis, \
    mean_difference_analysis, \
    residue_visualization, \
    distances_visualization

def normalize_histidine_name(feat_dict):
    return {
        key: [f.replace("HSE", "HIS").replace("HSD", "HIS").replace("HSP", "HIS") for f in val]
        for key, val in feat_dict.items()
    }



import re

def match_residue_numbers(feat_a, feat_b):
    def replace_resnum(target, source):
        # Replace the last number in target with the one from source
        resnum_source = re.search(r'\d+$', source)
        if resnum_source:
            return re.sub(r'\d+$', resnum_source.group(), target)
        else:
            return target  # fallback: no change

    matched_feat_b = {}
    for key in feat_a:
        if key in feat_b:
            matched_feat_b[key] = [
                replace_resnum(b, a) for a, b in zip(feat_a[key], feat_b[key])
            ]
        else:
            matched_feat_b[key] = feat_b.get(key, [])
    print(matched_feat_b)
    return matched_feat_b



def workflow_torsions_jsd(args, feat_a, feat_b, data_a, data_b, tors='bb'):

    # Use only features that are common to both ensembles
    if args.only_common_sctors and tors == 'sc':
        select_a, select_b = select_common_features(
            feat_a[tors + '-torsions'],
            feat_b[tors + '-torsions']
        )
    else:
        select_a = np.arange(len(feat_a[tors + '-torsions']))
        select_b = np.arange(len(feat_b[tors + '-torsions']))

    features_a = list(np.array(feat_a[tors + '-torsions'])[select_a])
    features_b = list(np.array(feat_b[tors + '-torsions'])[select_b])

    for i, (a, b) in enumerate(zip(features_a, features_b)):
        if a != b:
            print(f"Mismatch at index {i}: A = {a}, B = {b}")

    # Relative Entropy analysis with BB torsions
    relen = relative_entropy_analysis(
        list(np.array(feat_a[tors + '-torsions'])[select_a]),
        list(np.array(feat_b[tors + '-torsions'])[select_b]),
        data_a[tors + '-torsions'][:, select_a],
        data_b[tors + '-torsions'][:, select_b],
        bin_width=None, bin_num=10, verbose=False,
        override_name_check=args.override_name_check
    )
    names, jsd, kld_ab, kld_ba = relen

    # Save all results (per feature) in CSV files
    np.savetxt(
        args.out_results + '_' + tors + '-torsions_relative-entropy.csv', np.array(relen).T,
        fmt='%s', delimiter=', ', header='Name, JSD(A, B), KLD(A, B), KLD(B, A)'
    )

    # Save the Jensen-Shannon distance as "B factor" in a PDB file
    residue_visualization(
        names, jsd, args.ref_file_a,
        args.out_plots + "_" + tors + "-torsions_jsd.pdf",
        args.out_vispdb + "_" + tors + "-torsions_jsd.pdb",
        y_label='max. JS dist. of ' + tors.upper() + ' torsions'
    )

    # Print the features with the highest values
    print(tors.upper() + " torsions with the strongest deviations (JSD):")
    sf = sort_features(names, jsd)
    for f in sf[:args.print_num]:
        print(f[0], f[1])

    return names, jsd


def workflow_torsions_kss(args, feat_a, feat_b, data_a, data_b, tors='bb'):

    # Use only features that are common to both ensembles
    if args.only_common_sctors and tors == 'sc':
        select_a, select_b = select_common_features(
            feat_a[tors + '-torsions'],
            feat_b[tors + '-torsions']
        )
    else:
        select_a = np.arange(len(feat_a[tors + '-torsions']))
        select_b = np.arange(len(feat_b[tors + '-torsions']))

    # Kolmogorov-Smirnov analysis with BB torsions
    ksana = kolmogorov_smirnov_analysis(
        list(np.array(feat_a[tors + '-torsions'])[select_a]),
        list(np.array(feat_b[tors + '-torsions'])[select_b]),
        data_a[tors + '-torsions'][:, select_a],
        data_b[tors + '-torsions'][:, select_b],
        verbose=False,
        override_name_check=args.override_name_check
    )
    names, kss, ksp = ksana

    # Save all results (per feature) in CSV files
    np.savetxt(
        args.out_results + '_' + tors + '-torsions_kolmogorov-smirnov.csv',
        np.array(ksana).T, fmt='%s', delimiter=', ', header='Name, KSS(A, B), p-value'
    )

    # Save the Kolmogorov-Smirnov statistic as "B factor" in a PDB file
    residue_visualization(
        names, kss, args.ref_file_a,
        args.out_plots + "_" + tors + "-torsions_kss.pdf",
        args.out_vispdb + "_" + tors + "-torsions_kss.pdb",
        y_label='max. KS stat. of ' + tors.upper() + ' torsions'
    )

    # Print the features with the highest values
    print(tors.upper() + " torsions with the strongest deviations (KSS):")
    sf = sort_features(names, kss)
    for f in sf[:args.print_num]:
        print(f[0], f[1])

    return names, kss


def workflow_torsions_ssi(args, feat_a, feat_b, data_a, data_b, tors='bb'):

    # Use only features that are common to both ensembles
    if args.only_common_sctors:
        select_a, select_b = select_common_features(
            feat_a[tors + '-torsions'],
            feat_b[tors + '-torsions']
        )
    else:
        select_a = np.arange(len(feat_a[tors + '-torsions']))
        select_b = np.arange(len(feat_b[tors + '-torsions']))

    multivar_res_feat_a, multivar_res_data_a = get_multivar_res_timeseries(
        feat_a, data_a, tors + '-torsions'
    )
    multivar_res_feat_b, multivar_res_data_b = get_multivar_res_timeseries(
        feat_b, data_b, tors + '-torsions'
    )

    # Find the state boundaries
    discrete_states_ab = get_discrete_states(
        multivar_res_data_a[tors + '-torsions'], multivar_res_data_b[tors + '-torsions']
    )

    # Run the main analysis
    ana = ssi_ensemble_analysis(
        multivar_res_feat_a[tors + '-torsions'], multivar_res_feat_b[tors + '-torsions'],
        multivar_res_data_a[tors + '-torsions'], multivar_res_data_b[tors + '-torsions'],
        discrete_states_ab, verbose=False, 
        override_name_check=args.override_name_check
    )
    resnames, ssi = ana

    # Save all results (per feature) in CSV files
    np.savetxt(
        args.out_results + '_' + tors + '-torsions_state-specific-information.csv',
        np.array(ana).T, fmt='%s', delimiter=', ', header='Name, SSI(A, B)'
    )

    # Save the state-specific information as "B factor" in a PDB file
    residue_visualization(
        resnames, ssi, args.ref_file_a,
        args.out_plots + "_" + tors + "-torsions_ssi.pdf",
        args.out_vispdb + "_" + tors + "-torsions_ssi.pdb",
        y_label='SSI of ' + tors.upper() + ' torsions'
    )

    # Print the features with the highest values
    print(tors.upper() + " torsions with the strongest deviations (SSI):")
    sf = sort_features(resnames, ssi)
    for f in sf[:args.print_num]:
        print(f[0], f[1])

    return resnames, ssi


# -------------#
# --- MAIN --- #
# -------------#

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--ref_file_a", type=str,
                        default='traj/rhodopsin_arrbound_receptor.gro')
    parser.add_argument("--trj_file_a", type=str,
                        default='traj/rhodopsin_arrbound_receptor.xtc')
    parser.add_argument("--ref_file_b", type=str,
                        default='traj/rhodopsin_gibound_receptor.gro')
    parser.add_argument("--trj_file_b", type=str,
                        default='traj/rhodopsin_gibound_receptor.xtc')
    parser.add_argument("--out_plots", type=str,
                        default='plots/rhodopsin_receptor')
    parser.add_argument("--out_vispdb", type=str,
                        default='vispdb/rhodopsin_receptor')
    parser.add_argument("--out_results", type=str,
                        default='results/rhodopsin_receptor')
    parser.add_argument("--start_frame", type=int, default=0)
    parser.add_argument("--print_num", type=int, default=12)
    parser.add_argument("--override_name_check",
                        dest="override_name_check", default=False, action="store_true")
    parser.add_argument("--only_common_sctors",
                        dest="only_common_sctors", default=False, action="store_true")
    args = parser.parse_args()

    # -- FEATURES -- #

    # Load Features
    feat_a_0, data_a = read_structure_features(
        args.ref_file_a, args.trj_file_a, args.start_frame
    )
    feat_b_0, data_b = read_structure_features(
        args.ref_file_b, args.trj_file_b, args.start_frame
    )

    feat_a = normalize_histidine_name(feat_a_0)
    feat_b = normalize_histidine_name(feat_b_0)

    # Report dimensions
    print('Feature dimensions from', args.trj_file_a)
    for k in data_a.keys():
        print(k, data_a[k].shape)
    print('Feature dimensions from', args.trj_file_b)
    for k in data_b.keys():
        print(k, data_b[k].shape)

    # -- TORSIONS -- #

    print('BACKBONE TORSIONS')

    names, jsd = workflow_torsions_jsd(
        args, feat_a, feat_b, data_a, data_b, tors='bb'
    )
    names, kss = workflow_torsions_kss(
        args, feat_a, feat_b, data_a, data_b, tors='bb'
    )
    names, ssi = workflow_torsions_ssi(
        args, feat_a, feat_b, data_a, data_b, tors='bb'
    )

    print('SIDECHAIN TORSIONS')

    names, jsd = workflow_torsions_jsd(
        args, feat_a, feat_b, data_a, data_b, tors='sc'
    )
    names, kss = workflow_torsions_kss(
        args, feat_a, feat_b, data_a, data_b, tors='sc'
    )
    names, ssi = workflow_torsions_ssi(
        args, feat_a, feat_b, data_a, data_b, tors='sc'
    )

    # -- BACKBONE C-ALPHA DISTANCES --

    print('BACKBONE C-ALPHA DISTANCES')
    print(args.override_name_check)
    if args.override_name_check:
        print('matching')
        feat_b = match_residue_numbers(feat_a, feat_b)

    # Relative entropy analysis for C-alpha distances
    relen = relative_entropy_analysis(
        feat_a['bb-distances'], feat_b['bb-distances'],
        data_a['bb-distances'], data_b['bb-distances'],
        bin_width=0.01, verbose=False
    )
    names, jsd, kld_ab, kld_ba = relen

    # Save all results (per feature) in a CSV file
    np.savetxt(args.out_results + '_bb-distances_relative-entropy.csv', np.array(relen).T,
               fmt='%s', delimiter=', ', header='Name, JSD(A, B), KLD(A, B), KLD(B, A)')

    # Print the features with the highest values
    print("Backbone C-alpha distances with the strongest deviations:")
    sf = sort_features(names, jsd)
    for f in sf[:args.print_num]:
        print(f[0], f[1])

    # Visualize the deviations in a matrix plot
    print(args.out_plots)
    matrix = distances_visualization(names, jsd,
                                     args.out_plots + "_bb-distances-distributions_jsd.pdf",
                                     vmin=0.0, vmax=1.0)

    # Difference-of-the-mean analysis for C-alpha distances
    meanda = mean_difference_analysis(feat_a['bb-distances'], feat_b['bb-distances'],
                                      data_a['bb-distances'], data_b['bb-distances'],
                                      verbose=False)
    names, avg, diff = meanda

    # Save all results (per feature) in a CSV file
    np.savetxt(args.out_results + '_bb-distances_difference-of-mean.csv', np.array(meanda).T,
               fmt='%s', delimiter=', ', header='Name, average, difference')

    # Sort the distances by their differences
    print("Backbone C-alpha distances with the strongest differences of their mean value:")
    sf = sort_features(names, diff)
    for f in sf[:args.print_num]:
        print(f[0], f[1])

    # Visualize the deviations in a matrix plot
    matrix = distances_visualization(
        names, diff, args.out_plots + "_bb-diststances_difference-of-mean.pdf"
    )
