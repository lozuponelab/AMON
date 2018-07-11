# TODO: make remove measured cos w/o reaction option
# TODO: make heatmap for significance of pathway enrichments
# TODO: test running with cos measured

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from scipy.stats import hypergeom
from os import path
from statsmodels.sandbox.stats.multicomp import multipletests
import numpy as np
from biom import load_table
import seaborn as sns

from microMetabPred.parse_KEGG import parse_ko, parse_rn, parse_co, parse_pathway, get_kegg_record_dict


sns.set()


def p_adjust(pvalues, method='fdr_bh'):
    res = multipletests(pvalues, method=method)
    return np.array(res[1], dtype=float)


def read_in_ids(file_loc):
    """
    Read in kos from whitespace separated list (.txt), tsv with KOs as row headers (.tsv/.csv) or biom table (.biom).
    """
    if file_loc.endswith('.txt'):
        return [i.strip() for i in open(file_loc).read().split()]
    elif file_loc.endswith('.tsv') or file_loc.endswith('.csv'):
        return list(pd.read_table(file_loc, sep=None, index_col=0).columns)
    elif file_loc.endswith('.biom'):
        return list(load_table(file_loc).ids(axis='observation'))
    else:
        raise ValueError('Input file %s does not have a valid file ending.')


def get_rns_from_kos(list_of_kos, ko_dict):
    reaction_set = list()
    for ko in list_of_kos:
        try:
            ko_record = ko_dict[ko]
            if 'RN' in ko_record['DBLINKS']:
                reaction_set += ko_record['DBLINKS']['RN']
        except KeyError:
            pass
    return reaction_set


def get_products_from_rns(list_of_rns, rn_dict):
    return set([co for rn in list_of_rns for co in rn_dict[rn]['EQUATION'][1]])


def make_compound_origin_table(cos_produced, other_cos_produced=None, cos_measured=None):
    table = pd.DataFrame(index=['microbes', 'host', 'detected'])
    if other_cos_produced is None:
        other_cos_produced = []
    if cos_measured is None:
        cos_measured = []
    for co in set(cos_produced) | set(other_cos_produced) | set(cos_measured):
        table[co] = [co in cos_produced, co in other_cos_produced, co in cos_measured]
    table = table.transpose()
    # get rid of any all false columns
    table = table[table.columns[table.sum().astype(bool)]]
    return table


def make_venn(bac_cos, host_cos=None, measured_cos=None, output_loc=None):
    if host_cos is not None and measured_cos is None:
        _ = venn2((set(bac_cos), set(host_cos)),
                  ("Compounds predicted\nproduced by bacteria", "Compounds predicted\nproduced by host"),)
    elif host_cos is None and measured_cos is not None:
        _ = venn2((set(bac_cos), set(host_cos)),
                  ("Compounds predicted\nproduced by bacteria", "Compounds measured"))
    else:
        _ = venn3((set(measured_cos), set(bac_cos), set(host_cos)),
                  ("Compounds measured", "Compounds predicted\nproduced by bacteria",
                      "Compounds predicted\nproduced by host"))
    if output_loc is not None:
        plt.savefig(output_loc, bbox_inches='tight', dpi=300)


def get_pathways_from_cos(co_dict):
    pathway_list = list()
    for co_record in co_dict.values():
        if 'PATHWAY' in co_record:
            pathway_list += [pathway[0] for pathway in co_record['PATHWAY']]
    return set(pathway_list)


def get_pathway_to_co_dict(pathway_dict, no_drug=True, no_glycan=True):
    pathway_to_co_dict = {pathway_record['NAME']: [compound[0] for compound in pathway_record['COMPOUND']]
                          for pathway_record in pathway_dict.values() if 'COMPOUND' in pathway_record}
    if no_drug:
        pathway_to_co_dict = {pathway: [co for co in cos if not co.startswith('D')]
                              for pathway, cos in pathway_to_co_dict.items()}
    if no_glycan:
        pathway_to_co_dict = {pathway: [co for co in cos if not co.startswith('G')]
                              for pathway, cos in pathway_to_co_dict.items()}
    return pathway_to_co_dict


def calculate_enrichment(cos, co_pathway_dict, min_pathway_size=10):
    all_cos = set([co for co_list in co_pathway_dict.values() for co in co_list])
    pathway_names = list()
    pathway_data = list()
    for pathway, pathway_cos in co_pathway_dict.items():
        pathway_present = set(pathway_cos)
        if len(pathway_present) > min_pathway_size:
            overlap = set(cos) & pathway_present
            prob = hypergeom.sf(len(overlap), len(all_cos), len(pathway_present), len(set(cos)))
            pathway_names.append(pathway)
            pathway_data.append([len(pathway_present), len(overlap), prob])
    enrichment_table = pd.DataFrame(pathway_data, index=pathway_names,
                                    columns=["pathway size", "overlap", "probability"])
    enrichment_table['adjusted probability'] = p_adjust(enrichment_table.probability)
    return enrichment_table.sort_values('adjusted probability')


def make_enrichment_clustermap(microbe_enrichment_p, host_enrichment_p, output_dir, min_p=.2, log=False):
    enriched_pathways = microbe_enrichment_p.loc[microbe_enrichment_p < min_p].index
    enriched_pathways = enriched_pathways | host_enrichment_p.loc[host_enrichment_p < min_p].index
    enrichment_p_df = pd.concat((microbe_enrichment_p, host_enrichment_p), axis=1, sort=True)
    enrichment_p_df = enrichment_p_df.loc[enriched_pathways]
    enrichment_p_df.columns = ("Predicted Bacterial Only", "Predicted Human Only")
    if log:
        enrichment_p_df = np.log(enrichment_p_df)
    g = sns.clustermap(enrichment_p_df, col_cluster=False, figsize=(2, 12), cmap="Blues_r", method="average")
    _ = plt.setp(g.ax_heatmap.get_xticklabels(), rotation=340, fontsize=12, ha="left")
    _ = plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=12)
    plt.savefig(path.join(output_dir, 'enrichment_heatmap.png'), dpi=500, bbox_inches='tight')


def main(kos_loc, output_dir, compounds_loc=None, other_kos_loc=None, detected_only=False, rxn_compounds_only=False,
         ko_file_loc=None, rn_file_loc=None, co_file_loc=None, pathway_file_loc=None):
    # read in all kos and get records
    kos = read_in_ids(kos_loc)
    all_kos = kos
    if other_kos_loc is not None:
        other_kos = read_in_ids(other_kos_loc)
        all_kos = all_kos + other_kos
    else:
        other_kos = None
    ko_dict = get_kegg_record_dict(set(all_kos), parse_ko, ko_file_loc)

    # get all reactions from kos
    rns = get_rns_from_kos(kos, ko_dict)
    all_rns = rns
    if other_kos is not None:
        other_rns = get_rns_from_kos(other_kos, ko_dict)
        all_rns = all_rns + other_rns
    else:
        other_rns = None

    # Get reactions from KEGG and pull kos produced
    rn_dict = get_kegg_record_dict(set(all_rns), parse_rn, rn_file_loc)
    cos_produced = get_products_from_rns(rns, rn_dict)
    if other_rns is not None:
        other_cos_produced = get_products_from_rns(other_rns, rn_dict)
    else:
        other_cos_produced = None

    # read in compounds that were measured if available
    if compounds_loc is not None:
        cos_measured = read_in_ids(compounds_loc)
    else:
        cos_measured = None
    origin_table = make_compound_origin_table(cos_produced, other_cos_produced, cos_measured)
    origin_table.to_csv(path.join(output_dir, 'origin_table.tsv'), sep='\t')

    # Get set of compounds
    if detected_only:
        if other_kos_loc is None:
            all_cos = set(cos_measured) | cos_produced
        else:
            all_cos = set(cos_measured) | cos_produced | other_cos_produced
    else:
        if other_cos_produced is None:
            all_cos = cos_produced
        else:
            all_cos = cos_produced | other_cos_produced

    # Get compound data from kegg
    co_dict = get_kegg_record_dict(all_cos, parse_co, co_file_loc)

    # remove compounds without reactions if required
    if rxn_compounds_only:
        cos_with_rxn = list()
        for compound, record in co_dict.items():
            if 'REACTION' in record:
                cos_with_rxn.append(compound)
        cos_measured = set(cos_measured) & set(cos_with_rxn)

    # Make venn diagram
    if compounds_loc is not None or other_kos_loc is not None:
        make_venn(cos_produced, other_cos_produced, cos_measured, path.join(output_dir, 'venn.png'))

    # Filter compounds down to only cos measured for cos produced and other cos produced
    if detected_only:
        cos_produced = set(cos_produced) & set(cos_measured)
        if other_cos_produced is not None:
            other_cos_produced = set(other_cos_produced) & set(cos_measured)

    # find compounds unique to microbes and to host if host included
    if other_cos_produced is not None:
        original_cos_produced = cos_produced
        cos_produced = cos_produced - other_cos_produced
        other_cos_produced = other_cos_produced - original_cos_produced

    # Get pathway info from pathways in compounds
    all_pathways = get_pathways_from_cos(co_dict)
    pathway_dict = get_kegg_record_dict(all_pathways, parse_pathway, pathway_file_loc)
    pathway_to_compound_dict = get_pathway_to_co_dict(pathway_dict, no_glycan=False)

    # calculate enrichment
    pathway_enrichment_df = calculate_enrichment(cos_produced, pathway_to_compound_dict)
    pathway_enrichment_df.to_csv(path.join(output_dir, 'bacteria_enrichment.tsv'), sep='\t')
    if other_kos_loc is not None:
        other_pathway_enrichment_df = calculate_enrichment(other_cos_produced, pathway_to_compound_dict)
        other_pathway_enrichment_df.to_csv(path.join(output_dir, 'host_enrichment.tsv'), sep='\t')
    else:
        other_pathway_enrichment_df = None

    # plot enrichment
    if other_pathway_enrichment_df is not None:
        make_enrichment_clustermap(pathway_enrichment_df['adjusted probability'],
                                   other_pathway_enrichment_df['adjusted probability'], output_dir)
