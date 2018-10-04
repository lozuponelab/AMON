import pytest
import asyncio

from KEGG_parser.parsers import parse_ko, parse_rn, parse_co, parse_pathway, parse_organism
from KEGG_parser.downloader import get_from_kegg_flat_file, get_from_kegg_api, get_kegg_record_dict, \
    get_kegg_link_from_api


@pytest.fixture()
def list_of_kos():
    return ['K00001', 'K00002']


@pytest.fixture()
def list_of_rxns():
    return ['R02124', 'R06927']


@pytest.fixture()
def loop():
    return asyncio.get_event_loop()


def test_get_from_kegg_kos(loop, list_of_kos):
    raw_entries = get_from_kegg_api(loop, list_of_kos, parse_ko)
    assert len(raw_entries) == 2


def test_get_from_kegg_rxns(loop, list_of_rxns):
    raw_entries = get_from_kegg_api(loop, list_of_rxns, parse_rn)
    assert len(raw_entries) == 2


@pytest.fixture()
def ko_raw_record():
    return "ENTRY       K00001                      KO\n" \
           "NAME        E0.0.0.0\n" \
           "DEFINITION  a fake gene\n" \
           "PATHWAY     ko00000 a fake pathway\n" \
           "DISEASE     H00000 A bad one\n" \
           "CLASS       Metabolism; Carbohydrate Metabolism; Glycolysis / Gluconeogenesis[PATH:ko00010]\n" \
           "DBLINKS     RN: R00000\n" \
           "            COG: COG0000\n" \
           "GENES       HSA: hsa00000\n" \
           "REFERENCE\n" \
           "  AUTHORS   Fake G.\n" \
           "  TITLE     Not Real\n" \
           "  JOURNAL   Nurture (2001)\n" \
           "  SEQUENCE  [fke:FK_0000]"


def test_parse_ko(ko_raw_record):
    ko_record = parse_ko(ko_raw_record)
    assert len(ko_record) == 8
    assert tuple(ko_record['DBLINKS']['RN']) == tuple(['R00000'])
    assert tuple(ko_record['DBLINKS']['COG']) == tuple(['COG0000'])


@pytest.fixture()
def rn_raw_record():
    return "ENTRY       R00000                      Reaction\n" \
           "NAME        a fake reaction\n" \
           "DEFINITION  reactant 1 + reactant 2 <=> product 1 + product 2\n" \
           "EQUATION    C00000(m+n) + 3 C00001 <=> C00002(m) + (n+1) C00003\n" \
           "RPAIR       RP05682  C00001_C00011 leave\n" \
           "            RP07361  C00011_C01010 leave" \
           "            RP07385  C00014_C01010 leave" \
           "            RP13530  C00011_C01010 leave" \
           "            RP13531  C00014_C01010 leave\n" \
           "ENZYME      0.0.0.0\n" \
           "DBLINKS     HSA: hsa00000\n" \
           "REFERENCE\n" \
           "  AUTHORS   Fake G.\n" \
           "  TITLE     Not Real\n" \
           "  JOURNAL   Nurture (2001)\n"


def test_parse_rn(rn_raw_record):
    rn_record = parse_rn(rn_raw_record)
    assert len(rn_record) == 7
    assert tuple(rn_record['EQUATION'][0]) == tuple(['C00000', 'C00001'])
    assert tuple(rn_record['EQUATION'][1]) == tuple(['C00002', 'C00003'])


@pytest.fixture()
def co_raw_record():
    return "ENTRY       C00000                      Compound\n" \
           "NAME        a fake compound;\n" \
           "            another name for this compound\n" \
           "FORMULA     C6H12O6\n" \
           "EXACT_MASS  1.000000\n" \
           "MOL_WEIGHT  1.000001\n" \
           "REMARK      Same as G00000\n" \
           "REACTION    R00906 R00922 R00935 R00936 R00937 R00977 R00994 R01000\n" \
           "            R01033 R01034 R01036 R01061 R01088 R01093 R01094 R01130\n" \
           "PATHWAY     map00000 Fake pathway\n" \
           "DBLINKS     ChEBI: 00000\n" \
           "ATOM        1\n" \
           "            1   O0  O    22.1250  -16.2017\n" \
           "BOND        0\n" \
           "REFERENCE\n" \
           "  AUTHORS   Fake G.\n" \
           "  TITLE     Not Real\n" \
           "  JOURNAL   Nurture (2001)\n" \
           "SEQUENCE    [fke:FK_0000]\n"


def test_parse_co(co_raw_record):
    co_record = parse_co(co_raw_record)
    assert len(co_record) == 9
    assert len(co_record['REACTION']) == 16
    assert co_record['PATHWAY'][0] == ('map00000', 'Fake pathway')


@pytest.fixture()
def gl_raw_record():
    return "ENTRY       G00000                      Glycan\n" \
           "COMPOSITION (GlcA)2 (GlcNAc)2 (S)2\n" \
           "MASS        1000\n" \
           "REACTION    R00000 R00001\n" \
           "PATHWAY     map00000  Fake pathway\n" \
           "            map01100  Metabolic pathways\n" \
           "MODULE      M00000  A fake module\n" \
           "NODE        6" \
           "            1   GlcA     13.6  -2.1\n" \
           "            2   GlcNAc    4.6  -2.1\n" \
           "            3   S        -1.4   1.9\n" \
           "            4   GlcA     -4.4  -2.1\n" \
           "            5   S       -10.4   1.9\n" \
           "            6   GlcNAc  -13.4  -2.1\n" \
           "EDGE        5\n" \
           "            1     2:a1    1:4  \n" \
           "            2     3       2:6  \n" \
           "            3     4:b1    2:4  \n" \
           "            4     5       4:2  \n" \
           "            5     6:a1    4:4  \n"


def test_parse_gl(gl_raw_record):
    co_record = parse_co(gl_raw_record)
    assert len(co_record) == 6
    assert len(co_record['REACTION']) == 2
    assert co_record['PATHWAY'][0] == ('map00000', 'Fake pathway')


@pytest.fixture()
def pathway_raw_record():
    return "ENTRY       ko00000                      Pathway\n" \
           "NAME        a fake pathway\n" \
           "DESCRIPTION A long description of what this pathway does\n" \
           "CLASS       09100 Metabolism; 09101 Carbohydrate metabolism\n" \
           "PATHWAY_MAP ko00000 a fake pathway\n" \
           "COMPOUND    C00000 A fake compound\n" \
           "            C99999 Another fake compound\n" \
           "DBLINKS     GO: 0019852 0019572\n" \
           "ORGANISM    Homo sapiens (human) [GN:hsa]\n" \
           "REFERENCE\n" \
           "  AUTHORS   Michal G.\n" \
           "  TITLE     Biochemical Pathways\n" \
           "  JOURNAL   Wiley (1999)\n"


def test_parse_pathway(pathway_raw_record):
    pathway_record = parse_pathway(pathway_raw_record)
    assert len(pathway_record) == 7
    assert set([compound[0] for compound in pathway_record['COMPOUND']]) == {'C00000', 'C99999'}


@pytest.fixture()
def organism_raw_record():
    return "ENTRY       000000            CDS       H.sapiens\n" \
           "NAME        SBS6969, ASBS1234\n" \
           "DEFINITION  a fake guman gene\n" \
           "ORTHOLOGY   K00000 Some fake orthology\n" \
           "PATHWAY     hsa00000  Specific Stuff\n" \
           "            hsa99999  The Place\n" \
           "CLASS       Things; Stuff; Specific Stuff [PATH:hsa00000]\n" \
           "            This; Must Be; The Place [PATH:hsa99999]\n" \
           "POSITION    1p1.1\n" \
           "MOTIF       Pfam: ASDF, JKL1, HYDE\n" \
           "            PROSITE: A_PROSITE_ID ANOTHER_PROSITE_ID\n" \
           "DBLINKS     NCBI-GI: 00000000\n" \
           "            HGNC: 00000\n" \
           "AASEQ       12\n" \
           "            QRHMIRHTGDGPYKCQECGKAFDRPSLFRIHERTHTGEKPHECKQCGKAFISFTNFQSHM\n" \
           "            IRHTGDGPYKCKVCGRAFIFPSYVRKHERTHTGEKPYECNKCGKTFSSSSNVRTHERTHT\n" \
           "            GEKPYECKECGKAFISLPSVRRHMIKHTGDGPYKCQVCGRAFDCPSSFQIHERTHTGEKP\n" \
           "NTSEQ       1234\n" \
           "            agaactcacactggtgagaaaccctatgcatgtccggaatgtgggaaagccttcatttct\n" \
           "            ctcccaagtgttcgaagacacatgattaagcacactggagatggaccatataaatgtcag\n" \
           "            gaatgtgggaaagcctttgatcgcccaagtttatttcagatacatgaaagaactcacact\n" \
           "            ggagagaaaccctatgaatgtcaggaatgtgcaaaagctttcatttctcttccaagtttt\n"


def test_parse_organism(organism_raw_record):
    organism_record = parse_organism(organism_raw_record)
    assert len(organism_record) == 9
    assert len(organism_record['ORTHOLOGY']) == 2
    assert organism_record['ORTHOLOGY'][0] == 'K00000'


@pytest.fixture()
def ko_flat_file(tmpdir_factory, ko_raw_record):
    fn = tmpdir_factory.mktemp("data").join("ko")
    fn.write('%s///\n' % '///\n'.join([ko_raw_record] * 3))
    return str(fn)


def test_get_from_kegg_flat_file(ko_flat_file):
    organism_records = get_from_kegg_flat_file(ko_flat_file)
    assert len(organism_records) == 3


def test_get_kegg_record_dict(list_of_kos, ko_flat_file):
    ko_dict_web = get_kegg_record_dict(list_of_kos, parse_ko)
    assert len(ko_dict_web) == 2
    ko_dict_local = get_kegg_record_dict(list_of_kos, parse_ko, ko_flat_file)
    assert len(ko_dict_local) == 1


def test_get_kegg_link_from_api():
    link_dict = get_kegg_link_from_api('ko', 'hsa')
    assert tuple(link_dict.keys())[0].startswith('K')
