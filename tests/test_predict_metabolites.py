import pytest
from numpy.testing import assert_allclose
import pandas as pd
from biom.table import Table
from os.path import isfile
import numpy as np

from AMON.predict_metabolites import p_adjust, read_in_ids, make_compound_origin_table, get_rns_from_kos, \
                                     get_products_from_rns, get_pathways_from_cos, get_pathway_to_co_dict, \
                                     make_venn, calculate_enrichment, make_enrichment_clustermap, \
                                     make_kegg_mapper_input, reverse_dict_of_lists, merge_dicts_of_lists,\
                                     get_unique_from_dict_of_lists


@pytest.fixture()
def list_of_pvalues():
    return .001, .05, .5


def test_p_adjust(list_of_pvalues):
    p_adj = p_adjust(list_of_pvalues)
    assert len(p_adj) == 3
    assert_allclose(p_adj, (.003, .075, .5))


@pytest.fixture()
def list_of_kos():
    return ['K00001', 'K00002', 'K00003']


@pytest.fixture()
def ids_txt(tmpdir_factory, list_of_kos):
    fn = tmpdir_factory.mktemp("data").join("ko_list.txt")
    open(str(fn), 'w').write('%s\n' % '\n'.join(list_of_kos))
    return str(fn)


def test_read_in_ids_txt(ids_txt, list_of_kos):
    ids = read_in_ids(ids_txt, name='Sample1')
    assert len(ids) == 1
    assert ids['Sample1'] == set(list_of_kos)


@pytest.fixture()
def ids_csv(tmpdir_factory, list_of_kos):
    fn = tmpdir_factory.mktemp("data").join("ko_list.csv")
    columns = list_of_kos
    index = ('Sample1', 'Sample2')
    data = np.array([[1, 1, 0],
                     [1, 0, 1]])
    df = pd.DataFrame(data, index=index, columns=columns)
    df.to_csv(str(fn))
    return str(fn)


def test_read_in_ids_csv(ids_csv, list_of_kos):
    ids = read_in_ids(ids_csv, name='Sample1')
    assert len(ids) == 1
    assert set(ids['Sample1']) == set(list_of_kos)


def test_read_in_ids_csv_keep_separated(ids_csv):
    sample_dict = read_in_ids(ids_csv, keep_separated=True)
    assert len(sample_dict) == 2
    assert 'Sample1' in sample_dict
    assert set(sample_dict['Sample1']) == {'K00001', 'K00002'}
    assert 'Sample2' in sample_dict
    assert set(sample_dict['Sample2']) == {'K00001', 'K00003'}


@pytest.fixture()
def ids_tsv(tmpdir_factory, list_of_kos):
    fn = tmpdir_factory.mktemp("data").join("ko_list.tsv")
    columns = list_of_kos
    index = ('Sample1', 'Sample2')
    data = [[1, 1, 0],
            [1, 0, 1]]
    df = pd.DataFrame(data, index=index, columns=columns)
    df.to_csv(str(fn), sep='\t')
    return str(fn)


def test_read_in_ids_tsv(ids_tsv, list_of_kos):
    ids = read_in_ids(ids_tsv, name='Sample1')
    assert len(ids) == 1
    assert set(ids['Sample1']) == set(list_of_kos)


def test_read_in_ids_tsv_keep_separated(ids_tsv):
    sample_dict = read_in_ids(ids_tsv, keep_separated=True)
    assert 2 == len(sample_dict)
    assert 'Sample1' in sample_dict
    assert set(sample_dict['Sample1']) == {'K00001', 'K00002'}
    assert 'Sample2' in sample_dict
    assert set(sample_dict['Sample2']) == {'K00001', 'K00003'}


@pytest.fixture()
def ids_biom(tmpdir_factory, list_of_kos):
    fn = tmpdir_factory.mktemp("data").join("ko_list.biom")
    observations = list_of_kos + ['K99999']
    samples = ('Sample1', 'Sample2')
    data = np.array([[1, 1, 0, 0],
                     [1, 0, 1, 0]])
    table = Table(data.transpose(), observations, samples)
    table.to_json('testing', open(str(fn), 'w'))
    return str(fn)


def test_read_in_ids_biom(ids_biom, list_of_kos):
    ids = read_in_ids(ids_biom, keep_separated=False, name='Sample1')
    assert len(ids) == 1
    assert set(ids['Sample1']) == set(list_of_kos)


def test_read_in_ids_biom_keep_separated(ids_biom):
    sample_dict = read_in_ids(ids_biom, keep_separated=True)
    assert len(sample_dict) == 2
    assert 'Sample1' in sample_dict
    assert set(sample_dict['Sample1']) == {'K00001', 'K00002'}
    assert 'Sample2' in sample_dict
    assert set(sample_dict['Sample2']) == {'K00001', 'K00003'}


def test_read_in_ids_bad_ending():
    with pytest.raises(ValueError):
        ids = read_in_ids('fake.xkcd')


@pytest.fixture()
def ko_dict():
    ko1 = {'ENTRY': 'K00001', 'REACTION': ['R00000', 'R00001']}
    ko2 = {'ENTRY': 'K00002', 'DBLINKS': {'COG': ['COG0000']}}
    return {'K00001': ko1, 'K00002': ko2}


@pytest.fixture()
def dict_of_kos(list_of_kos):
    return {'Sample1': {'K00001', 'K00002'},
            'Sample2': {'K00001'}
    }


def test_get_rns_from_kos(dict_of_kos, ko_dict):
    rns = get_rns_from_kos(dict_of_kos, ko_dict)
    assert len(rns) == 2
    assert len(rns['Sample1']) == 2
    assert len(rns['Sample2']) == 2


@pytest.fixture()
def list_of_rns():
    return ['R00000', 'R00001']


@pytest.fixture()
def rn_dict():
    rn1 = {'ENTRY': 'R00000', 'EQUATION': (('C00001', 'C00002'), ('C00003', 'C00004'))}
    rn2 = {'ENTRY': 'R00001', 'EQUATION': (('C00006', 'C00005'), ('C00003', 'C00007'))}
    return {'R00000': rn1, 'R00001': rn2}


@pytest.fixture()
def dict_of_rns(list_of_rns):
    return {'Sample1': list_of_rns}


def test_get_products_from_rns(dict_of_rns, rn_dict):
    products = get_products_from_rns(dict_of_rns, rn_dict)
    assert len(products) == 1
    assert len(products['Sample1']) == 3


def test_reverse_dict_of_lists():
    dict_of_lists = {'Sample1': ['C00001', 'C00002'],
                     'Sample2': ['C00001', 'C00003']}
    reversed_dict_of_lists = reverse_dict_of_lists(dict_of_lists)
    assert len(reversed_dict_of_lists) == 3
    assert len(reversed_dict_of_lists['C00001']) == 2
    assert len(reversed_dict_of_lists['C00002']) == 1
    assert len(reversed_dict_of_lists['C00003']) == 1


@pytest.fixture()
def list_of_cos():
    return ['C00002', 'C00003', 'C00005']


@pytest.fixture()
def list_of_other_cos():
    return ['C00001', 'C00002', 'C00004']


@pytest.fixture()
def dict_of_cos(list_of_cos, list_of_other_cos):
    return {'Sample1': list_of_cos,
            'Sample2': list_of_other_cos}


@pytest.fixture()
def list_of_measured_cos():
    return ['C00001', 'C00003', 'C00006']


def test_make_compound_origin_table(dict_of_cos, list_of_measured_cos):
    # two samples
    table1 = make_compound_origin_table(dict_of_cos)
    assert table1.shape == (5, 2)
    assert tuple(table1.sum(axis=0)) == (3, 3)
    # two samples and measured compounds
    table2 = make_compound_origin_table(dict_of_cos, list_of_measured_cos)
    assert table2.shape == (5, 3)
    assert tuple(table2.sum(axis=0)) == (3, 3, 2)


def test_merge_dicts_of_lists(dict_of_kos, dict_of_cos):
    merged_dicts = merge_dicts_of_lists(dict_of_kos, dict_of_cos)
    assert len(merged_dicts) == 2
    assert set(merged_dicts['Sample1']) == {'K00001', 'K00002', 'C00002', 'C00003', 'C00005'}
    assert set(merged_dicts['Sample2']) == {'K00001', 'C00001', 'C00002', 'C00004'}


def test_make_kegg_mapper_input(list_of_cos, dict_of_cos, list_of_measured_cos):
    # bacteria only
    kegg_mapper_table1 = make_kegg_mapper_input({'Sample1': list_of_cos})
    # host but not measured
    kegg_mapper_table2 = make_kegg_mapper_input(dict_of_cos)
    # measured but no host
    kegg_mapper_table3 = make_kegg_mapper_input({'Sample1': list_of_cos}, detected_ids=list_of_measured_cos)
    # all
    kegg_mapper_table4 = make_kegg_mapper_input(dict_of_cos, list_of_measured_cos)
    assert kegg_mapper_table4.shape == (6,)
    assert kegg_mapper_table4['C00001'] == 'yellow,orange'
    assert kegg_mapper_table4['C00002'] == 'green'
    assert kegg_mapper_table4['C00004'] == 'yellow'
    assert kegg_mapper_table4['C00005'] == 'blue'


def test_get_unique_from_dict_of_lists(dict_of_cos):
    unique_dict_of_cos = get_unique_from_dict_of_lists(dict_of_cos)
    assert len(unique_dict_of_cos) == 2
    assert set(unique_dict_of_cos['Sample1']) == {'C00003', 'C00005'}
    assert set(unique_dict_of_cos['Sample2']) == {'C00001', 'C00004'}


def test_make_venn(dict_of_cos, list_of_measured_cos, tmpdir):
    with pytest.raises(ValueError):
        make_venn({'Sample1': list_of_measured_cos})
    p = tmpdir.mkdir('test_venn')
    b_h_venn_path = str(p.join('b_h_venn.png'))
    make_venn(dict_of_cos, output_loc=b_h_venn_path)
    assert isfile(b_h_venn_path)
    b_m_venn_path = str(p.join('b_m_venn.png'))
    make_venn({'Sample1': list_of_measured_cos}, measured_cos=list_of_measured_cos, output_loc=b_m_venn_path)
    assert isfile(b_m_venn_path)
    b_h_m_venn_path = str(p.join('b_m_h_venn.png'))
    make_venn(dict_of_cos, list_of_measured_cos, output_loc=b_h_m_venn_path)
    assert isfile(b_h_m_venn_path)


@pytest.fixture()
def co_dict():
    co1 = {'ENTRY': 'C00000', 'PATHWAY': (('ko00001', 'Name1'), ('ko00003', 'Name3'))}
    co2 = {'ENTRY': 'C00001', 'PATHWAY': (('ko00006', 'Name6'), ('ko00003', 'Name3'))}
    return {'C00000': co1, 'C00001': co2}


def test_get_pathways_from_cos(co_dict):
    pathways = get_pathways_from_cos(co_dict)
    assert len(pathways) == 3
    assert pathways == {'ko00001', 'ko00003', 'ko00006'}


@pytest.fixture()
def pathway_dict():
    pathway1 = {'ENTRY': 'ko00000', 'NAME': 'One fake pathway', 'COMPOUND': (('C00001', ''), ('C00002', ''))}
    pathway2 = {'ENTRY': 'ko00001', 'NAME': 'Another fake pathway',
                'COMPOUND': (('C00001', 'A compound'), ('C00005', ''), ('C00003', ''), ('C00007', ''))}
    return {'ko00000': pathway1, 'ko00001': pathway2}


def test_get_pathway_to_co_dict(pathway_dict):
    pathway_to_co_dict = get_pathway_to_co_dict(pathway_dict)
    assert len(pathway_to_co_dict) == 2
    assert tuple(pathway_to_co_dict['One fake pathway']) == ('C00001', 'C00002')
    assert tuple(pathway_to_co_dict['Another fake pathway']) == ('C00001', 'C00005', 'C00003', 'C00007')


@pytest.fixture()
def pathway_co_dict():
    return {'one fake pathway': ('C00002', 'C00004', 'C00006'),
            'two fake pathway': ('C00001', 'C00002', 'C00005'),
            'three fake pathway': ('C00001', 'C00004', 'C00006', 'C00007', 'C00008')}


def test_calculate_enrichment(list_of_cos, pathway_co_dict):
    enrichment = calculate_enrichment(list_of_cos, pathway_co_dict, min_pathway_size=0)
    assert enrichment.shape == (3, 4)
    assert tuple(enrichment.index) == ('two fake pathway', 'one fake pathway', 'three fake pathway')


@pytest.fixture()
def p_values():
    return pd.Series((.0000005, .05, .5), index=('two fake pathway', 'one fake pathway', 'three fake pathway'))


@pytest.fixture()
def host_p_values():
    return pd.Series((.001, .2, .7), index=('one fake pathway', 'two fake pathway', 'three fake pathway'))


@pytest.fixture()
def enrichment_dfs():
    data1 = [[1, .001],
            [1, .2],
            [10, .7],
            [3, .0005]]
    df1 = pd.DataFrame(data1, index=('one fake pathway', 'two fake pathway', 'three fake pathway', 'four fake pathway'),
                       columns=('overlap', 'p-value'))
    data2 = [[1, .0000005],
            [1, .05],
            [10, .5],
            [1, .002]]
    df2 = pd.DataFrame(data2, index=('one fake pathway', 'two fake pathway', 'three fake pathway', 'four fake pathway'),
                       columns=('overlap', 'p-value'))
    return {'Sample1': df1,
            'Sample2': df2}


def test_make_enrichment_clustermap(enrichment_dfs, tmpdir):
    p = tmpdir.mkdir('test_venn')
    enrichment_path = str(p.join('enrichment_heatmap.png'))
    make_enrichment_clustermap(enrichment_dfs, 'p-value', output_loc=enrichment_path)
    assert isfile(enrichment_path)
