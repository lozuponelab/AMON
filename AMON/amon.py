#!/usr/bin/env python

import argparse
from AMON.get_kegg_files import get_kegg_files
from AMON.predict_metabolites import main as amon_main

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Primary inputs
    parser.add_argument('-i', '--gene_set', help="KEGG KO's from bacterial community or organism of interest in the "
                                                 "form of a white space separated list, a tsv or csv with KO ids as "
                                                 "column names or a biom file with KO ids as observations",
                        required=True)
    parser.add_argument('-o', '--output_dir', help="directory to store output", required=True)
    parser.add_argument('--detected_compounds', help="list of compounds detected via metabolomics")
    parser.add_argument('--other_gene_set', help="white space separated list of KEGG KO's from the host, another "
                                                 "organism or other environment")
    # Names for gene sets
    parser.add_argument('--gene_set_name', help="Name to use for first gene set (should have no spaces, underscore "
                                                "separated)")
    parser.add_argument('--other_gene_set_name', help="Name to use for second gene set (should have no spaces, "
                                                      "underscore separated)")
    # Options
    parser.add_argument('--keep_separated', help='If input in biom or tabular format keep samples separate for '
                                                 'analysis', action='store_true', default=False)
    parser.add_argument('--samples_are_columns', help='If data is in tabular format, by default genes are columns and '
                                                      'samples rows, to indicate that samples are columns and genes '
                                                      'are rows use this flag', action='store_true', default=False)
    # Filters
    parser.add_argument('--detected_only', help="only use detected compounds in enrichment analysis",
                        action='store_true', default=False)
    parser.add_argument('--rn_compound_only', help="only use compounds with associated reactions", action='store_true',
                        default=False)
    parser.add_argument('--unique_only', help='only use compounds that are unique to a sample in enrichment',
                        action='store_true', default=False)
    # Local KEGG files
    parser.add_argument('--save_entries', help='Save json file of KEGG entries at all levels used in analysis for '
                                               'deeper analysis', action='store_true', default=False)
    parser.add_argument("--force-download-kegg", help="Re-download KEGG flat files from cloud",
                        action="store_true", default=False)
    
    args = parser.parse_args()

    name1 = args.gene_set_name or "gene_set_1"
    name2 = args.other_gene_set_name or "gene_set_2"

    if args.detected_only and args.detected_compounds is None:
        raise ValueError(
            "Cannot use --detected_only without --detected_compounds"
        )

    # Resolve KEGG flat files (cloud + cache)
    kegg_paths = get_kegg_files(force_download=args.force_download_kegg)

    amon_main(
        args.gene_set,
        args.output_dir,
        args.other_gene_set,
        args.detected_compounds,
        name1,
        name2,
        args.keep_separated,
        args.samples_are_columns,
        args.detected_only,
        args.rn_compound_only,
        args.unique_only,
        ko_file_loc=kegg_paths["ko"],
        rn_file_loc=kegg_paths["reaction"],
        co_file_loc=kegg_paths["compound"],
        pathway_file_loc=kegg_paths["pathway"],
        enzyme_file_loc=kegg_paths["enzyme"],
        write_json=args.save_entries,
    )


if __name__ == "__main__":
    main()