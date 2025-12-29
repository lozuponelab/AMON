#!/usr/bin/env python

import argparse

from KEGG_parser.parsers import parse_organism
from KEGG_parser.downloader import get_from_kegg_flat_file, get_kegg_link_from_api

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', help="KEGG organism identifier or KEGG organism flat file",
                        required=True)
    parser.add_argument('-o', '--output', help="Output file of new line separated list of KOs from genome",
                        required=True)
    parser.add_argument('--from_flat_file', help="Indicates that input is a flat flile to be parsered directly",
                        action='store_true', default=False)

    args = parser.parse_args()

    kegg_org_file_loc = args.input
    output_file = args.output

    if args.from_flat_file:
        org_records = get_from_kegg_flat_file(kegg_org_file_loc, parser=parse_organism)
        # for org_record in org_records:
        kos = [org_record['ORTHOLOGY'][0] for org_record in org_records if 'ORTHOLOGY' in org_record]
    else:
        link_dict = get_kegg_link_from_api('ko', args.input)
        kos = link_dict.keys()

    with open(output_file, 'w') as f:
        f.write('%s\n' % '\n'.join(kos))
