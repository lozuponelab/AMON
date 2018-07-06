# TODO: make remove measured cos w/o reaction option
# TODO: make read in biom table of kos or cos an option
# TODO: accept kegg files to parse as input
# TODO: rewrite get dict methods to all be one method that takes a parser and file loc

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from scipy.stats import hypergeom
from os import path
from statsmodels.sandbox.stats.multicomp import multipletests
import numpy as np

from microMetabPred.parse_KEGG import get_from_kegg_api, parse_ko, parse_rn, parse_co, parse_pathway, \
                                      get_from_kegg_flat_file


def p_adjust(pvalues, method='fdr_bh'):
    res = multipletests(pvalues, method=method)
    return np.array(res[1], dtype=float)


def parse_whitespace_sep(file_loc):
    return [i.strip() for i in open(file_loc).readlines()]


def get_ko_dict(list_of_kos, ko_file_loc=None):
    if ko_file_loc is None:
        ko_raw_records = get_from_kegg_api(list_of_kos)
        ko_records = [parse_ko(i) for i in ko_raw_records]
    else:
        ko_records = get_from_kegg_flat_file(ko_file_loc, list_of_kos, parse_ko)
    print("%s ko records acquired" % len(ko_records))
    return {ko_record['ENTRY']: ko_record for ko_record in ko_records}


def get_rns_from_kos(list_of_kos, ko_dict):
    reaction_set = list()
    for ko in list_of_kos:
        try:
            ko_record = ko_dict[ko]
            if 'RN' in ko_record['DBLINKS']:
                reaction_set += ko_record['DBLINKS']['RN']
        except KeyError:
            # print('KO %s does not exist in KEGG?' % ko)
            pass
    return reaction_set


def get_reaction_dict(list_of_reactions, rn_file_loc=None):
    if rn_file_loc is None:
        rn_raw_records = get_from_kegg_api(list_of_reactions)
        rn_records = [parse_rn(i) for i in rn_raw_records]
    else:
        rn_records = get_from_kegg_flat_file(rn_file_loc, list_of_reactions, parse_rn)
    print("%s rn records acquired" % len(rn_records))
    return {rn_record['ENTRY']: rn_record for rn_record in rn_records}


def get_products_from_rns(list_of_rns, rn_dict):
    return set([co for rn in list_of_rns for co in rn_dict[rn]['EQUATION'][1]])


def make_compound_origin_table(cos_produced, other_cos_produced=None, cos_measured=None):
    table = pd.DataFrame(index=['microbe', 'host', 'detected'])
    if other_cos_produced is None:
        other_cos_produced = []
    if cos_measured is None:
        cos_measured = []
    for co in set(cos_produced) ^ set(other_cos_produced) ^ set(cos_measured):
        table[co] = [co in cos_produced, co in other_cos_produced, co in cos_measured]
    table = table.transpose()
    # get rid of any all false columns
    table = table[table.columns[table.sum().astype(bool)]]
    return table


def make_venn(bac_cos, host_cos=None, measured_cos=None, output_loc=None):
    if host_cos is not None and measured_cos is None:
        venn = venn2((set(bac_cos), set(host_cos)),
                     ("Compounds predicted produced by bacteria", "Compounds predicted produced by host"))
    if host_cos is None and measured_cos is not None:
        venn = venn2((set(bac_cos), set(host_cos)),
                     ("Compounds predicted produced by bacteria", "Compounds measured"))
    else:
        venn = venn3((set(measured_cos), set(bac_cos), set(host_cos)),
                     ("Compounds measured", "Compounds predicted produced by bacteria",
                      "Compounds predicted produced by host"))
    if output_loc is not None:
        plt.savefig(output_loc)


def get_compound_dict(list_of_compounds, co_file_loc):
    if co_file_loc is None:
        co_raw_records = get_from_kegg_api(list_of_compounds)
        co_records = [parse_co(i) for i in co_raw_records]
    else:
        co_records = get_from_kegg_flat_file(co_file_loc, list_of_compounds, parse_co)
    print("%s co records acquired" % len(co_records))
    return {co_record['ENTRY']: co_record for co_record in co_records}


def get_pathways_from_cos(co_dict):
    pathway_list = list()
    for co_record in co_dict.values():
        if 'PATHWAY' in co_record:
            pathway_list += [pathway[0] for pathway in co_record['PATHWAY']]
    return pathway_list


def get_pathway_dict(list_of_pathways, pathway_file_loc):
    if pathway_file_loc is None:
        co_raw_records = get_from_kegg_api(list_of_pathways)
        pathway_records = [parse_co(i) for i in co_raw_records]
    else:
        pathway_records = get_from_kegg_flat_file(pathway_file_loc, list_of_pathways, parse_pathway)
    print("%s pathway records acquired" % len(pathway_records))
    return {pathway_record['ENTRY']: pathway_record for pathway_record in pathway_records}


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


def calculate_enrichment(cos, co_pathway_dict, min_pathway_size=3):
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
    enrichment_table = pd.DataFrame(pathway_data, index=pathway_names, columns=["pathway size", "overlap", "probability"])
    enrichment_table['adjusted probability'] = p_adjust(enrichment_table.probability)
    return enrichment_table.sort_values('adjusted probability')


def main(kos_loc, output_dir, compounds_loc=None, other_kos_loc=None, detected_only=False,
         ko_file_loc=None, rn_file_loc=None, co_file_loc=None, pathway_file_loc=None):
    # read in all kos and get records
    kos = parse_whitespace_sep(kos_loc)
    all_kos = kos
    if other_kos_loc is not None:
        other_kos = parse_whitespace_sep(other_kos_loc)
        all_kos = all_kos + other_kos
    ko_dict = get_ko_dict(set(all_kos), ko_file_loc)

    # get all reactions from kos
    rns = get_rns_from_kos(kos, ko_dict)
    all_rns = rns
    if other_kos_loc is not None:
        other_rns = get_rns_from_kos(other_kos, ko_dict)
        all_rns = all_rns + rns
    rn_dict = get_reaction_dict(set(all_rns), rn_file_loc)
    cos_produced = get_products_from_rns(rns, rn_dict)
    if other_kos_loc is not None:
        other_cos_produced = get_products_from_rns(other_rns, rn_dict)
    else:
        other_cos_produced = None

    # read in compounds that were measured if available
    if compounds_loc is not None:
        cos_measured = parse_whitespace_sep(compounds_loc)
    else:
        cos_measured = None
    origin_table = make_compound_origin_table(cos_produced, other_cos_produced, cos_measured)
    origin_table.to_csv(path.join(output_dir, 'origin_table.tsv'), sep='\t')
    if compounds_loc is not None or other_kos_loc is not None:
        make_venn(cos_produced, other_cos_produced, cos_measured, path.join(output_dir, 'venn.png'))

    # calculate enrichment
    if detected_only:
        cos_produced = cos_produced & cos_measured
        if other_cos_produced is not None:
            other_cos_produced = other_cos_produced & cos_measured
        all_cos = cos_measured
    else:
        if other_cos_produced is None:
            all_cos = cos_produced
        else:
            all_cos = cos_produced & other_cos_produced
    co_dict = get_compound_dict(all_cos, co_file_loc)
    all_pathways = get_pathways_from_cos(co_dict)
    pathway_dict = get_pathway_dict(all_pathways, pathway_file_loc)
    pathway_to_compound_dict = get_pathway_to_co_dict(pathway_dict, no_glycan=False)
    enrichment_table = calculate_enrichment(cos_produced, pathway_to_compound_dict)
    enrichment_table.to_csv(path.join(output_dir, 'bacteria_enrichment.tsv'), sep='\t')
    if other_cos_produced is not None:
        other_cos_enrichment_table = calculate_enrichment(other_cos_produced, pathway_to_compound_dict)
        other_cos_enrichment_table.to_csv(path.join(output_dir, 'host_enrichment'), sep='\t')
