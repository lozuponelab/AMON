import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles
from scipy.stats import hypergeom
from os import path, makedirs
from statsmodels.sandbox.stats.multicomp import multipletests
import numpy as np
from biom import load_table
import seaborn as sns
import json
from collections import defaultdict, OrderedDict
from datetime import datetime
from warnings import warn

from KEGG_parser.parsers import parse_ko, parse_rn, parse_co, parse_pathway
from KEGG_parser.downloader import get_kegg_record_dict

sns.set()

# TODO: take multiple files


class Logger(OrderedDict):
    """"""
    def __init__(self, output):
        super(Logger, self).__init__()
        self.output_file = output
        self['start time'] = datetime.now()

    def output_log(self):
        with open(self.output_file, 'w') as f:
            self['finish time'] = datetime.now()
            self['elapsed time'] = self['finish time'] - self['start time']
            for key, value in self.items():
                f.write(key + ': ' + str(value) + '\n')


def p_adjust(pvalues, method='fdr_bh'):
    res = multipletests(pvalues, method=method)
    return np.array(res[1], dtype=float)


def read_in_ids(file_loc, keep_separated=False, samples_are_columns=False, name=None):
    """
    Read in kos from whitespace separated list (.txt), tsv with KOs as row headers (.tsv/.csv) or biom table (.biom).
    """
    if file_loc.endswith('.txt'):
        if name is None:
            raise ValueError('Name must be given if giving .txt list')
        return {name: set([i.strip() for i in open(file_loc).read().split()])}
    elif (file_loc.endswith('.tsv') or file_loc.endswith('.csv')) and keep_separated:
        genome_table = pd.read_table(file_loc, sep=None, index_col=0, engine='python')
        samples_dict = dict()
        if samples_are_columns:
            genome_table = genome_table.transpose()
        for sample in genome_table.index:
            samples_dict[sample] = set(genome_table.columns[genome_table.loc[sample].astype(bool)])
        return samples_dict
    elif file_loc.endswith('.tsv') or file_loc.endswith('.csv'):
        if name is None:
            raise ValueError('Name must be given if giving .tsv or .csv and not separating')
        return {name: set(pd.read_table(file_loc, sep=None, index_col=0, engine='python').columns)}
    elif file_loc.endswith('.biom') and keep_separated:
        id_table = load_table(file_loc)
        samples_dict = dict()
        for data, sample, _ in id_table.iter(axis='sample'):
            samples_dict[sample] = set(id_table.ids(axis='observation')[data.astype(bool)])
        return samples_dict
    elif file_loc.endswith('.biom'):
        if name is None:
            raise ValueError('Name must be given if giving .biom and not separating')
        id_table = load_table(file_loc)
        # first remove KO's which aren't present in any samples
        ids_to_keep = id_table.ids(axis='observation')[id_table.sum(axis='observation') > 0]
        id_table.filter(ids_to_keep, axis='observation', inplace=True)
        return {name: set(id_table.ids(axis='observation'))}
    else:
        raise ValueError('Input file %s does not have a parsable file ending.')


def get_rns_from_kos(dict_of_kos: dict, ko_dict: dict):
    sample_rns = dict()
    for sample, list_of_kos in dict_of_kos.items():
        reaction_set = list()
        for ko in list_of_kos:
            try:
                ko_record = ko_dict[ko]
                if 'REACTION' in ko_record.keys():
                    rxn_ids = [rxn[0] for rxn in ko_record['REACTION']]
                    reaction_set += rxn_ids
            except KeyError:
                pass
        sample_rns[sample] = reaction_set
    return sample_rns


def get_products_from_rns(dict_of_rns: dict, rn_dict: dict):
    return {sample: set([co for rn in list_of_rns for co in rn_dict[rn]['EQUATION'][1]])
            for sample, list_of_rns in dict_of_rns.items()}


def reverse_dict_of_lists(dict_of_lists):
    reversed_dict = defaultdict(list)
    for key, list_ in dict_of_lists.items():
        for item in list_:
            reversed_dict[item].append(key)
    return reversed_dict


def make_compound_origin_table(sample_cos_produced: dict, cos_measured=None):
    columns = list(sample_cos_produced.keys())
    rows = list()
    cos_to_samples_dict = reverse_dict_of_lists(sample_cos_produced)
    for co, samples in cos_to_samples_dict.items():
        rows.append([sample in samples for sample in columns])
    table = pd.DataFrame(rows, index=cos_to_samples_dict.keys(), columns=columns)
    if cos_measured is not None:
        table['detected'] = [co in cos_measured for co in table.index]
    return table


def merge_dicts_of_lists(*dicts):
    merged_dicts = defaultdict(list)
    for dict_ in dicts:
        for key, list_ in dict_.items():
            merged_dicts[key] += list_
    return merged_dicts


def make_kegg_mapper_input(sample_ids, detected_ids=None, origin_colors=('blue', 'green', 'yellow'),
                           detected_color='orange'):
    samples = list(sample_ids.keys())
    microbe_ids = sample_ids[samples[0]]
    if len(samples) == 2:
        host_ids = sample_ids[samples[1]]
    else:
        host_ids = ()
    if detected_ids is None:
        detected_ids = ()
    ids = list()
    colors = list()
    for id_ in set(microbe_ids) | set(host_ids) | set(detected_ids):
        # save id
        ids.append(id_)
        # check where id is present
        microbe_present = id_ in microbe_ids
        host_present = id_ in host_ids
        detected_present = id_ in detected_ids
        origin_color = None
        detect_color = None
        if microbe_present and host_present:
            origin_color = origin_colors[1]
        elif microbe_present:
            origin_color = origin_colors[0]
        elif host_present:
            origin_color = origin_colors[2]
        else:
            pass
        if detected_present:
            detect_color = detected_color
        color = ''
        if origin_color is not None:
            color += origin_color
        if detect_color is not None:
            color += ',%s' % detect_color
        colors.append(color)
    df = pd.Series(colors, index=ids)
    return df


def make_venn(sample_cos_produced, measured_cos=None, output_loc=None, name1='gene_set_1', name2='gene_set_2'):
    samples = list(sample_cos_produced.keys())
    bac_cos = sample_cos_produced[samples[0]]
    if len(samples) == 2:
        host_cos = sample_cos_produced[samples[1]]
    else:
        host_cos = None
    if host_cos is None and measured_cos is None:
        raise ValueError("Must give host_cos or measured_cos to make venn diagram")
    if host_cos is not None and measured_cos is None:
        _ = venn2((set(bac_cos), set(host_cos)),
                  ("Compounds predicted\nproduced by %s" % name1.replace('_', ' '),
                   "Compounds predicted\nproduced by %s" % name2.replace('_', ' ')),
                  set_colors=('white',)*2)
        _ = venn2_circles((set(bac_cos), set(host_cos)), linestyle='solid')
    elif host_cos is None and measured_cos is not None:
        _ = venn2((set(bac_cos), set(measured_cos)),
                  ("Compounds predicted\nproduced by %s" % name1.replace('_', ' '), "Compounds measured"),
                  set_colors=('white',)*2)
        _ = venn2_circles((set(bac_cos), set(measured_cos)), linestyle='solid')
    else:
        _ = venn3((set(measured_cos), set(bac_cos), set(host_cos)),
                  ("Compounds measured", "Compounds predicted\nproduced by %s" % name1.replace('_', ' '),
                   "Compounds predicted\nproduced by %s" % name2.replace('_', ' ')),
                  set_colors=('white',)*3)
        _ = venn3_circles((set(measured_cos), set(bac_cos), set(host_cos)), linestyle='solid')
    if output_loc is not None:
        plt.savefig(output_loc, bbox_inches='tight', dpi=300)
    else:
        plt.show()


def get_pathways_from_cos(co_dict):
    pathway_list = list()
    for co_record in co_dict.values():
        if 'PATHWAY' in co_record:
            pathway_list += [pathway[0] for pathway in co_record['PATHWAY']]
    return set(pathway_list)


def get_unique_from_dict_of_lists(dict_of_lists):
    unique_dict_of_lists = dict()
    for key, list_ in dict_of_lists.items():
        all_other_values = set([value for other_key, other_list in dict_of_lists.items() for value in other_list
                                if other_key != key])
        unique_dict_of_lists[key] = set(list_) - all_other_values
    return unique_dict_of_lists


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
    # if 0 rows in enrichment table, return None
    # otherwise there's a zero division error during p adj
    if len(enrichment_table.index)==0:
        warn("No pathways were identified from the KOs provided."
             "Please verify that your KOs are valid (e.g., formatted as K02041).")
        return None
    
    
    enrichment_table['adjusted probability'] = p_adjust(enrichment_table.probability)
    if np.any((enrichment_table['adjusted probability'] < .05) & (enrichment_table['overlap'] == 0)):
        return None
    else:
        return enrichment_table.sort_values('adjusted probability')


def make_enrichment_clustermap(pathway_enrichment_dfs: dict, key, output_loc, min_p=.1, log=False):
    enrichment_p_df = pd.DataFrame.from_dict({sample: pathway_enrichment_df[key] for sample, pathway_enrichment_df in
                                              pathway_enrichment_dfs.items()})
    enrichment_p_df = enrichment_p_df.loc[enrichment_p_df.index[(enrichment_p_df<min_p).sum(axis=1) > 0]]
    enrichment_p_df = enrichment_p_df[enrichment_p_df.columns[(enrichment_p_df<min_p).sum(axis=0) > 0]]
    if log:
        enrichment_p_df = np.log(enrichment_p_df)
    g = sns.clustermap(enrichment_p_df, col_cluster=False, figsize=(2, 12), cmap="Blues_r", method="average")
    _ = plt.setp(g.ax_heatmap.get_xticklabels(), rotation=340, fontsize=12, ha="left")
    _ = plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=12)
    plt.savefig(output_loc, dpi=500, bbox_inches='tight')


def main(kos_loc, output_dir, other_kos_loc=None, compounds_loc=None, name1='gene_set_1', name2='gene_set_2',
         keep_separated=False, samples_are_columns=False, detected_only=False, rxn_compounds_only=False,
         unique_only=True, ko_file_loc=None, rn_file_loc=None, co_file_loc=None, pathway_file_loc=None,
         write_json=False, try_async=False):
    # create output dir to throw error quick
    makedirs(output_dir)
    logger = Logger(path.join(output_dir, "AMON_log.txt"))

    # read in all kos and get records
    sample_kos = read_in_ids(kos_loc, keep_separated=keep_separated,
                             samples_are_columns=samples_are_columns, name=name1)
    logger['kos_loc'] = path.abspath(kos_loc)
    if other_kos_loc is not None:
        sample_kos.update(read_in_ids(other_kos_loc, keep_separated=keep_separated,
                                      samples_are_columns=samples_are_columns, name=name2))
        logger['other_kos_loc'] = path.abspath(other_kos_loc)
    all_kos = set([value for values in sample_kos.values() for value in values])
    logger['Number of samples'] = len(sample_kos)
    logger['Total number of KOs'] = len(all_kos)

    ko_dict = get_kegg_record_dict(set(all_kos), parse_ko, ko_file_loc, try_async=try_async)
    if write_json:
        open(path.join(output_dir, 'ko_dict.json'), 'w').write(json.dumps(ko_dict))
        logger['KO json location'] = path.abspath(path.join(output_dir, 'ko_dict.json'))

    # get all reactions from kos
    sample_rns = get_rns_from_kos(sample_kos, ko_dict)
    all_rns = set([value for values in sample_rns.values() for value in values])
    logger['Total number of reactions'] = len(all_rns)

    # get reactions from kegg
    rn_dict = get_kegg_record_dict(set(all_rns), parse_rn, rn_file_loc, try_async=try_async)
    if write_json:
        open(path.join(output_dir, 'rn_dict.json'), 'w').write(json.dumps(rn_dict))
        logger['RN json location'] = path.abspath(path.join(output_dir, 'rn_dict.json'))

    # Get reactions from KEGG and pull cos produced
    sample_cos_produced = get_products_from_rns(sample_rns, rn_dict)

    # read in compounds that were measured if available
    if compounds_loc is not None:
        cos_measured = list(read_in_ids(compounds_loc, name='Compounds', keep_separated=False).values())[0]
        logger['compounds_loc'] = path.abspath(compounds_loc)
    else:
        cos_measured = None

    # make compound origin table
    origin_table = make_compound_origin_table(sample_cos_produced, cos_measured)

    # get rid of any all false columns
    origin_table = origin_table[origin_table.columns[origin_table.sum().astype(bool)]]
    origin_table.to_csv(path.join(output_dir, 'origin_table.tsv'), sep='\t')
    logger['Origin table location'] = path.abspath(path.join(output_dir, 'origin_table.tsv'))

    # make kegg mapper input if 2 or fewer samples
    if len(sample_cos_produced) <= 2:
        kegg_mapper_input = make_kegg_mapper_input(merge_dicts_of_lists(sample_kos, sample_cos_produced), cos_measured)
        kegg_mapper_input.to_csv(path.join(output_dir, 'kegg_mapper.tsv'), sep='\t')
        logger['KEGG mapper location'] = path.abspath(path.join(output_dir, 'kegg_mapper.tsv'))

    # Get full set of compounds
    all_cos_produced = set([value for values in sample_cos_produced.values() for value in values])
    logger['Number of cos produced across samples'] = len(all_cos_produced)
    if detected_only:
        all_cos_produced = set(all_cos_produced) | set(cos_measured)
        logger['Number of cos produced and detected'] = len(all_cos_produced)

    # Get compound data from kegg
    co_dict = get_kegg_record_dict(all_cos_produced, parse_co, co_file_loc, try_async=try_async)
    if write_json:
        open(path.join(output_dir, 'co_dict.json'), 'w').write(json.dumps(co_dict))

    # remove compounds without reactions if required
    if rxn_compounds_only:
        cos_with_rxn = list()
        for compound, record in co_dict.items():
            if 'REACTION' in record:
                cos_with_rxn.append(compound)
        cos_measured = set(cos_measured) & set(cos_with_rxn)

    # Make venn diagram
    if (compounds_loc is not None or len(sample_cos_produced) > 1) and len(sample_cos_produced) <= 2:
        make_venn(sample_cos_produced, cos_measured, path.join(output_dir, 'venn.png'))

    # Filter compounds down to only cos measured for cos produced and other cos produced
    if detected_only:
        sample_cos_produced = {sample: set(cos_produced) & set(cos_measured) for sample, cos_produced
                               in sample_cos_produced.items()}

    # find compounds unique to microbes and to host if host included
    if unique_only:
        sample_cos_produced = get_unique_from_dict_of_lists(sample_cos_produced)

    # Get pathway info from pathways in compounds
    all_pathways = [pathway.replace('map', 'ko') for pathway in get_pathways_from_cos(co_dict)]
    pathway_dict = get_kegg_record_dict(all_pathways, parse_pathway, pathway_file_loc, try_async=try_async)
    pathway_to_compound_dict = get_pathway_to_co_dict(pathway_dict, no_glycan=False)

    # calculate enrichment
    pathway_enrichment_dfs = dict()
    for sample, cos_produced in sample_cos_produced.items():
        pathway_enrichment_df = calculate_enrichment(cos_produced, pathway_to_compound_dict)
        if pathway_enrichment_df is not None:
            pathway_enrichment_df.to_csv(path.join(output_dir, '%s_compound_pathway_enrichment.tsv' % sample), sep='\t')
            logger['%s pathway enrichment'] = path.abspath(path.join(output_dir,
                                                                     '%s_compound_pathway_enrichment.tsv' % sample))
            pathway_enrichment_dfs[sample] = pathway_enrichment_df

    if len(pathway_enrichment_dfs) > 0:
        make_enrichment_clustermap(pathway_enrichment_dfs, 'adjusted probability',
                                   path.join(output_dir, 'enrichment_heatmap.png'))
        logger['Enrichment clustermap location'] = path.abspath(path.join(output_dir, 'enrichment_heatmap.png'))

    logger.output_log()
