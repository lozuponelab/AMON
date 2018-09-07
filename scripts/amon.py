import argparse

from AMON.predict_metabolites import main

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--gene_set', help="KEGG KO's from bacterial community or organism of interest in the "
                                                 "form of a white space separated list, a tsv or csv with KO ids as "
                                                 "column names or a biom file with KO ids as observations",
                        required=True)
    parser.add_argument('-o', '--output_dir', help="directory to store output", required=True)
    parser.add_argument('--detected_compounds', help="list of compounds detected via metabolomics")
    parser.add_argument('--other_gene_set', help="white space separated list of KEGG KO's from the host, another "
                                                 "organism or other environment")

    parser.add_argument('--detected_only', help="only use detected compounds in enrichment analysis",
                        action='store_true', default=False)
    parser.add_argument('--rn_compound_only', help="only use compounds with associated reactions", action='store_true',
                        default=False)

    parser.add_argument('--ko_file_loc', help='Location of ko file from KEGG FTP download')
    parser.add_argument('--rn_file_loc', help='Location of reaction file from KEGG FTP download')
    parser.add_argument('--co_file_loc', help='Location of compound file from KEGG FTP download')
    parser.add_argument('--pathway_file_loc', help='Location of pathway file from KEGG FTP download')
    parser.add_argument('--save_entries', help='Save json file of KEGG entries at all levels used in analysis for '
                                               'deeper analysis', action='store_true', default=False)

    parser.add_argument('--verbose', help="verbose output", action='store_true', default=False)

    args = parser.parse_args()
    kos_loc = args.gene_set
    output_dir = args.output_dir
    detected_compounds = args.detected_compounds
    other_kos_loc = args.other_gene_set
    detected_compounds_only = args.detected_only
    rn_compounds_only = args.rn_compound_only

    ko_file_loc = args.ko_file_loc
    rn_file_loc = args.rn_file_loc
    co_file_loc = args.co_file_loc
    pathway_file_loc = args.pathway_file_loc
    write_json = args.save_entries

    verbose = args.verbose

    if detected_compounds_only and detected_compounds is None:
        raise ValueError('Cannot have detected compounds only and not provide detected compounds')

    main(kos_loc, output_dir, detected_compounds, other_kos_loc, detected_compounds_only, rn_compounds_only,
         ko_file_loc=ko_file_loc, rn_file_loc=rn_file_loc, co_file_loc=co_file_loc, pathway_file_loc=pathway_file_loc,
         write_json=write_json, verbose=verbose)
