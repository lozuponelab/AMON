#!/usr/bin/env python

import argparse
from KEGG_parser.parsers import parse_organism
from KEGG_parser.downloader import (
    get_from_kegg_flat_file,
    get_kegg_link_from_api,
)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', help="KEGG organism identifier or KEGG organism flat file",
                        required=True)
    parser.add_argument('-o', '--output', help="Output file of new line separated list of KOs from genome",
                        required=True)
    parser.add_argument('--from_flat_file', help="Indicates that input is a flat flile to be parsered directly",
                        action='store_true', default=False)

    args = parser.parse_args()

    if args.from_flat_file:
        org_records = get_from_kegg_flat_file(
            args.input, parser=parse_organism
        )
        kos = [
            rec["ORTHOLOGY"][0]
            for rec in org_records
            if "ORTHOLOGY" in rec
        ]
    else:
        kos = get_kegg_link_from_api("ko", args.input).keys()

    with open(args.output, "w") as f:
        f.write("\n".join(kos) + "\n")


if __name__ == "__main__":
    main()