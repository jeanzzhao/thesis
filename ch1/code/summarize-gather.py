#! /usr/bin/env python
import os.path
import json
import csv
import sys
import argparse

import pandas
import sourmash
from sourmash import sourmash_args

class GatherInfo:
    headers = ["accession", "ref_f_unweighted", "ref_f_weighted", "n_matches"]

    def __init__(self, metag_acc, gather_csv):
        self.metag_acc = metag_acc
        self.gather_csv = gather_csv

    def calc(self):
        self.calc_ref_based_kmer_stuff(self.gather_csv)

    def get_row(self):
        xx = [self.metag_acc]
        for x in self.headers[1:]:
            val = getattr(self, x)
            xx.append(val)
        return xx

    def calc_ref_based_kmer_stuff(self, gather_csv):
        df = pandas.read_csv(gather_csv)
        row = df.tail(1).squeeze()
        sum_weighted_found = row['sum_weighted_found'] 
        total_weighted_hashes = row['total_weighted_hashes']
        # oops, this is the same as: print(df['f_unique_weighted'].sum())
        print('XXX', sum_weighted_found, total_weighted_hashes)
        print(f"for {self.metag_acc}:")
        print(f"total ref k-mers found (abund): {sum_weighted_found / total_weighted_hashes * 100:.1f}")
        print(f"total ref k-mers found (flat): {df['f_unique_to_query'].sum() * 100:.1f}")
        self.n_matches = len(df)
        self.ref_f_unweighted = df['f_unique_to_query'].sum()
        self.ref_f_weighted = sum_weighted_found / total_weighted_hashes


def main(argv):
    p = argparse.ArgumentParser(argv)
    p.add_argument('csvs', nargs='+')
    p.add_argument('-o', '--output-csv', help='output summary CSV here',
                   required=True)
    
    args = p.parse_args()

    results = []
    for gather_csv in args.csvs:
        acc = os.path.basename(gather_csv).split('.')[0]

        info = GatherInfo(acc, gather_csv)
        info.calc()
        results.append(info.get_row())

    with open(args.output_csv, 'w', newline='') as fp:
        w = csv.writer(fp)
        w.writerow(GatherInfo.headers)
        for rr in results:
            w.writerow(rr)


if __name__ == '__main__':
    sys.exit(main(sys.argv))
