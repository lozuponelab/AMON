import pytest

from microMetabPred.predict_metabolites import get_cos_from_kos,  make_compound_origin_table


@pytest.fixture()
def list_of_kos():
    return ['K00001', 'K00002']


def test_get_cos_from_kos(list_of_kos):
    cos = get_cos_from_kos(list_of_kos)
    assert type(cos) == set
    assert len(cos) > 0


@pytest.fixture()
def list_of_cos():
    return ['C00003', 'C00002']


@pytest.fixture()
def list_of_other_cos():
    return []


@pytest.fixture()
def list_of_measured_cos():
    return []


def test_make_compound_origin_table(list_of_cos):
    table = make_compound_origin_table(list_of_cos)
    assert table.shape == (2, 1)
