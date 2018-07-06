import argparse

from microMetabPred.predict_metabolites import main

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input', help="white space separated list of KEGG KO's from bacterial community",
                        required=True)
    parser.add_argument('-o', '--output_dir', help="directory to store output", required=True)
    parser.add_argument('--detected_compounds', help="list of compounds detected via metabolomics")
    parser.add_argument('--host_kos', help="white space separated list of KEGG KO's from host or other environment")
    parser.add_argument('--detected_only', help="only use detected metabolites in enrichment analysis",
                        action='store_true', default=False)

    args = parser.parse_args()
    kos_loc = args.input
    output_dir = args.output_dir
    detected_compounds = args.detected_compounds
    other_kos_loc = args.host_kos
    detected_only = args.detected_only

    if detected_only and detected_compounds is None:
        raise ValueError('Cannot have detected compounds only and not provide detected compounds')

    main(kos_loc, output_dir, detected_compounds, other_kos_loc, detected_only)
