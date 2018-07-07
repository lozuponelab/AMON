import pytest

from microMetabPred.predict_metabolites import  make_compound_origin_table, get_rns_from_kos, get_products_from_rns


@pytest.fixture()
def list_of_kos():
    return ['K00001', 'K00002', 'K00003']


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
    return ['C00003', 'C00002', 'C00005']


@pytest.fixture()
def list_of_other_cos():
    return ['C00001', 'C00002', 'C00004']


@pytest.fixture()
def list_of_measured_cos():
    return ['C00001', 'C00006', 'C00003']


def test_make_compound_origin_table(list_of_cos, list_of_other_cos, list_of_measured_cos):
    # bacteria only
    table1 = make_compound_origin_table(list_of_cos)
    assert table1.shape == (3, 1)
    assert tuple(table1.sum(axis=0)) == (3,)
    # host but not measured
    table2 = make_compound_origin_table(list_of_cos, other_cos_produced=list_of_other_cos)
    assert table2.shape == (5, 2)
    assert tuple(table2.sum(axis=0)) == (3, 3)
    # measured but no host
    table3 = make_compound_origin_table(list_of_cos, cos_measured=list_of_measured_cos)
    assert table3.shape == (5, 2)
    assert tuple(table3.sum(axis=0)) == (3, 3)
    # all
    table4 = make_compound_origin_table(list_of_cos, list_of_other_cos, list_of_measured_cos)
    assert table4.shape == (6, 3)
    assert tuple(table4.sum(axis=0)) == (3, 3, 3)
