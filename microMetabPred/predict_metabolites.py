# TODO: make remove measured cos w/o reaction option

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from scipy.stats import hypergeom
from os import path

from microMetabPred.parse_KEGG import get_from_kegg_api, parse_ko, parse_rn, parse_co, parse_pathway


def parse_whitespace_sep(file_loc):
    return set(i.strip() for i in open(file_loc).readlines())


def get_cos_from_kos(kos):
    ko_raw_records = get_from_kegg_api(kos)
    print("kos acquired")
    ko_records = [parse_ko(i) for i in ko_raw_records]
    reaction_set = set()
    for ko_record in ko_records:
        if 'RN' in ko_record:
            reaction_set.update(ko_record['DBLINKS']['RN'])
    rn_raw_records = get_from_kegg_api(reaction_set)
    print("rns acquired")
    rn_records = [parse_rn(i) for i in rn_raw_records]
    compounds_generated = set([co for rn_record in rn_records for co in rn_record['EQUATION'][1]])
    return compounds_generated


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


def get_co_pathway_dict(cos, no_drug=True, no_glycan=True):
    co_raw_records = get_from_kegg_api(cos)
    co_records = [parse_co(i) for i in co_raw_records]
    pathways = set([pathway[0] for co_record in co_records for pathway in co_record['PATHWAY']])
    pathway_raw_records = get_from_kegg_api([i.replace('map', 'ko') for i in pathways])
    pathway_records = [parse_pathway(i) for i in pathway_raw_records]
    co_pathway_dict = {pathway_record['NAME']: [compound[0] for compound in pathway_record['COMPOUND']]
                       for pathway_record in pathway_records if 'COMPOUND' in pathway_record}
    if no_drug:
        co_pathway_dict = {pathway: [co for co in cos if not co.startswith('D')]
                           for pathway, cos in co_pathway_dict.items()}
    if no_glycan:
        co_pathway_dict = {pathway: [co for co in cos if not co.startswith('G')]
                           for pathway, cos in co_pathway_dict.items()}
    return co_pathway_dict


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
    return pd.DataFrame(pathway_data, index=pathway_names, columns=["pathway size", "overlap", "probability"])


def main(kos_loc, output_dir, compounds_loc=None, other_kos_loc=None, detected_only=False):
    kos = parse_whitespace_sep(kos_loc)
    cos_produced = get_cos_from_kos(kos)
    if other_kos_loc is not None:
        other_kos = parse_whitespace_sep(other_kos_loc)
        other_cos_produced = get_cos_from_kos(other_kos)
    else:
        other_cos_produced = None
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
    co_pathway_dict = get_co_pathway_dict(all_cos)
    enrichment_table = calculate_enrichment(cos_produced, co_pathway_dict)
    enrichment_table.to_csv(path.join(output_dir, 'bacteria_enrichment.tsv'), sep='\t')
    if other_cos_produced is not None:
        other_cos_enrichment_table = calculate_enrichment(other_cos_produced, co_pathway_dict)
        other_cos_enrichment_table.to_csv(path.join(output_dir, 'host_enrichment'), sep='\t')
