#! /usr/bin/env python
import csv
import sys
import argparse
import random
from collections import defaultdict


def main():
    p = argparse.ArgumentParser()
    p.add_argument('summary_csv')
    p.add_argument('metadata_csv')
    p.add_argument('-o', '--output-subsampled', required=True)
    p.add_argument('-n', '--choose-num', default=3)
    p.add_argument('--seed', type=int, default=1, help='random number seed')
    args = p.parse_args()

    summary_by_accession = {}
    for row in csv.DictReader(open(args.summary_csv, newline='')):
        acc = row['accession']
        summary_by_accession[acc] = row

    fieldnames = list(row.keys())
    print(f"got {len(summary_by_accession)} summary entries.")

    metadata_by_accession = {}
    for row in csv.DictReader(open(args.metadata_csv, newline='')):
        acc = row['accession']
        metadata_by_accession[acc] = row

    print(f"got {len(metadata_by_accession)} metadata entries.")

    # produce subsets by column
    subset_options = defaultdict(set)
    col = 'biome3'
    for acc in summary_by_accession:
        val = metadata_by_accession[acc][col]
        subset_options[val].add(acc)

    print(f"got {len(subset_options)} subsets by column '{col}'")
    print(f"choosing {args.choose_num} rows from each subset.")

    print(f"setting random number seed to: {args.seed}")
    random.seed(args.seed)

    with open(args.output_subsampled, 'w', newline='') as fp:
        fieldnames.append('subsampled_on')
        w = csv.DictWriter(fp, fieldnames=fieldnames)
        w.writeheader()

        for subset_name, members in subset_options.items():
            subsample = random.choices(list(members), k=args.choose_num)
            if len(set(subsample)) < args.choose_num:
                print(f"WARNING: only got {len(set(subsample))} accs for {subset_name}")

            for acc in set(subsample):
                row = summary_by_accession[acc]
                new_row = dict(row)
                new_row['subsampled_on'] = subset_name

                w.writerow(new_row)
            



if __name__ == '__main__':
    sys.exit(main())
