import pytest
from numpy.testing import assert_allclose
import pandas as pd
from biom.table import Table
from os.path import isfile

from AMON.predict_metabolites import p_adjust, read_in_ids, make_compound_origin_table, get_rns_from_kos, \
                                     get_products_from_rns, get_pathways_from_cos, get_pathway_to_co_dict, \
                                     make_venn, calculate_enrichment, make_enrichment_clustermap, \
                                     make_kegg_mapper_input


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
    ids = read_in_ids(ids_txt)
    assert len(ids) == 3
    assert tuple(ids) == tuple(list_of_kos)


@pytest.fixture()
def ids_csv(tmpdir_factory, list_of_kos):
    fn = tmpdir_factory.mktemp("data").join("ko_list.csv")
    columns = list_of_kos
    index = ('Sample1', 'Sample2')
    data = [[1, 1, 1],
            [1, 1, 1]]
    df = pd.DataFrame(data, index=index, columns=columns)
    df.to_csv(str(fn))
    return str(fn)


def test_read_in_ids_csv(ids_csv, list_of_kos):
    ids = read_in_ids(ids_csv)
    assert len(ids) == 3
    assert tuple(ids) == tuple(list_of_kos)


@pytest.fixture()
def ids_tsv(tmpdir_factory, list_of_kos):
    fn = tmpdir_factory.mktemp("data").join("ko_list.tsv")
    columns = list_of_kos
    index = ('Sample1', 'Sample2')
    data = [[1, 1, 1],
            [1, 1, 1]]
    df = pd.DataFrame(data, index=index, columns=columns)
    df.to_csv(str(fn), sep='\t')
    return str(fn)


def test_read_in_ids_tsv(ids_tsv, list_of_kos):
    ids = read_in_ids(ids_tsv)
    assert len(ids) == 3
    assert tuple(ids) == tuple(list_of_kos)


@pytest.fixture()
def ids_biom(tmpdir_factory, list_of_kos):
    fn = tmpdir_factory.mktemp("data").join("ko_list.biom")
    observations = list_of_kos
    samples = ('Sample1', 'Sample2')
    data = [[1, 1, 1],
            [1, 1, 1]]
    table = Table(data, observations, samples)
    table.to_json('testing', open(str(fn), 'w'))
    return str(fn)


def test_read_in_ids_biom(ids_biom, list_of_kos):
    ids = read_in_ids(ids_biom)
    assert len(ids) == 3
    assert tuple(ids) == tuple(list_of_kos)


def test_read_in_ids_bad_ending():
    with pytest.raises(ValueError):
        ids = read_in_ids('fake.xkcd')


@pytest.fixture()
def ko_dict():
    ko1 = {'ENTRY': 'K00001', 'DBLINKS': {'RN': ['R00000', 'R00001']}}
    ko2 = {'ENTRY': 'K00002', 'DBLINKS': {'COG': ['COG0000']}}
    return {'K00001': ko1, 'K00002': ko2}


def test_get_rns_from_kos(list_of_kos, ko_dict):
    rns = get_rns_from_kos(list_of_kos, ko_dict)
    assert len(rns) == 2


@pytest.fixture()
def list_of_rns():
    return ['R00000', 'R00001']


@pytest.fixture()
def rn_dict():
    rn1 = {'ENTRY': 'R00000', 'EQUATION': (('C00001', 'C00002'), ('C00003', 'C00004'))}
    rn2 = {'ENTRY': 'R00001', 'EQUATION': (('C00006', 'C00005'), ('C00003', 'C00007'))}
    return {'R00000': rn1, 'R00001': rn2}


def test_get_products_from_rns(list_of_rns, rn_dict):
    products = get_products_from_rns(list_of_rns, rn_dict)
    assert len(products) == 3


@pytest.fixture()
def list_of_cos():
    return ['C00002', 'C00003', 'C00005']


@pytest.fixture()
def list_of_other_cos():
    return ['C00001', 'C00002', 'C00004']


@pytest.fixture()
def list_of_measured_cos():
    return ['C00001', 'C00003', 'C00006']


@pytest.fixture()
def test_make_compound_origin_table(list_of_cos, list_of_other_cos, list_of_measured_cos):
    # bacteria only
    table1 = make_compound_origin_table(list_of_cos)
    assert table1.shape == (3, 3)
    assert tuple(table1.sum(axis=0)) == (3, 0, 0)
    # host but not measured
    table2 = make_compound_origin_table(list_of_cos, other_cos_produced=list_of_other_cos)
    assert table2.shape == (5, 3)
    assert tuple(table2.sum(axis=0)) == (3, 3, 0)
    # measured but no host
    table3 = make_compound_origin_table(list_of_cos, cos_measured=list_of_measured_cos)
    assert table3.shape == (5, 3)
    assert tuple(table3.sum(axis=0)) == (3, 0, 3)
    # all
    table4 = make_compound_origin_table(list_of_cos, list_of_other_cos, list_of_measured_cos)
    assert table4.shape == (6, 3)
    assert tuple(table4.sum(axis=0)) == (3, 3, 3)
    return table1, table2, table3, table4


def test_make_kegg_mapper_input(test_make_compound_origin_table):
    origin_table1, origin_table2, origin_table3, origin_table4 = test_make_compound_origin_table
    # bacteria only
    kegg_mapper_table1 = make_kegg_mapper_input(origin_table1)
    # host but not measured
    kegg_mapper_table2 = make_kegg_mapper_input(origin_table2)
    # measured but no host
    kegg_mapper_table3 = make_kegg_mapper_input(origin_table3)
    # all
    kegg_mapper_table4 = make_kegg_mapper_input(origin_table4)
    assert kegg_mapper_table4.shape == (6,)
    assert kegg_mapper_table4['C00001'] == 'yellow,orange'
    assert kegg_mapper_table4['C00002'] == 'green'
    assert kegg_mapper_table4['C00005'] == 'blue'


def test_make_venn(list_of_cos, list_of_other_cos, list_of_measured_cos, tmpdir):
    with pytest.raises(ValueError):
        make_venn(list_of_cos)
    p = tmpdir.mkdir('test_venn')
    b_h_venn_path = str(p.join('b_h_venn.png'))
    make_venn(list_of_cos, host_cos=list_of_other_cos, output_loc=b_h_venn_path)
    assert isfile(b_h_venn_path)
    b_m_venn_path = str(p.join('b_m_venn.png'))
    make_venn(list_of_cos, measured_cos=list_of_measured_cos, output_loc=b_m_venn_path)
    assert isfile(b_m_venn_path)
    b_h_m_venn_path = str(p.join('b_m_h_venn.png'))
    make_venn(list_of_cos, list_of_other_cos, list_of_measured_cos, output_loc=b_h_m_venn_path)
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


def test_make_enrichment_clustermap(p_values, host_p_values, tmpdir):
    p = tmpdir.mkdir('test_venn')
    enrichment_path = str(p.join('enrichment_heatmap.png'))
    make_enrichment_clustermap(p_values, host_p_values, output_loc=enrichment_path)
    assert isfile(enrichment_path)
