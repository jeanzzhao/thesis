#! /usr/bin/env python
import sys
import argparse
import os
from pickle import dump, load
from pprint import pprint
import collections
import csv

from jsonapi_client import Session as APISession
from jsonapi_client import Modifier
import requests
import pandas as pd
# Data transformation
from functools import reduce
from collections import defaultdict

TEST=True

WORT_PATH="/group/ctbrowngrp/irber/data/wort-data/wort-sra/sigs/{acc}.sig"


def truncate_biome(biome, position):
    biome = biome.split(':')
    assert len(biome) >= position, biome
    biome = biome[1:position]
    return ':'.join(biome)


def get_runs_from_samples(samples_data):
    xx = []
    for item in samples_data:
        for item2 in item['data']:
            xx.append(item2['relationships']['runs']['links']['related'])
    return xx


def read_pickle(filename):
    if os.path.exists(filename):
        print(f"reading from '{filename}'")
        with open(filename, 'rb') as fp:
            x = load(fp)
            return x
        print('done!')

    return None


def main():
    p = argparse.ArgumentParser()
    p.add_argument('-o', '--output-spreadsheet')
    p.add_argument('--save-sig-paths')
    args = p.parse_args()

    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_colwidth', None)

    ##
    ## First get a list of all the biomes.
    ##

    # does it exist? if so, read. if not, grab, then save.
    biome_filename = '1-biomes.pickle'
    biome_json = read_pickle(biome_filename)

    all_biomes = []
    for result in biome_json:
        for record in result['data']:
            biome_name = record['id']
            sample_count = record['attributes']['samples-count']
            all_biomes.append((biome_name, sample_count))

    # sort and filter
    all_biomes.sort(key=lambda x: x[1])
    print("all:", len(all_biomes))
    all_biomes = [ (x, y) for (x, y) in all_biomes if y > 0 ]
    print("0-filtered:", len(all_biomes))
    #all_biomes = [ (x, y) for (x, y) in all_biomes if y >= 10 ]
    #print("10-filtered:", len(all_biomes))
    #all_biomes = [ (x, y) for (x, y) in all_biomes if y >= 50 ]
    #print("50-filtered:", len(all_biomes))

    # select only biomes that have four parts
    all_biomes = [ (x, y) for (x, y) in all_biomes if x.count(':') == 3 ]
    print("hierarchy filtered:", len(all_biomes))

    ##
    ## Then, for each biome, get a list of associated samples.
    ##

    biome_samples_filename = '2-biome-samples.pickle'
    samples_by_biome = read_pickle(biome_samples_filename)

    ##
    ## Extract the list of runs from each biome (no web request needed)
    ##

    runs_by_biome = defaultdict(list)
    for biome_name, samples_vv in samples_by_biome.items():
        runs = get_runs_from_samples(samples_vv)
        runs_by_biome[biome_name].extend(runs)

    runs_by_sample_filename = '3c-runs_by_sample.pickle'
    runs_by_sample = read_pickle(runs_by_sample_filename)

    sub_runs_by_biome = {}
    biome_runs_counter = collections.Counter()
    biome_illumina_counter = collections.Counter()
    total = 0
    total_illumina = 0
    for biome_name, runlist in runs_by_biome.items():
        x = []
        for sample_url in runlist:
            if sample_url in runs_by_sample:
                for item in runs_by_sample[sample_url]:
                    data = item['data']
                    for item in data:
                        attr = item.get('attributes')
                        if attr['experiment-type'] == 'metagenomic' and attr['instrument-platform'] == 'ILLUMINA':
                            print((biome_name, attr['accession'], attr['experiment-type'], attr['instrument-platform'], attr['instrument-model']))
                    
                            x.append((biome_name, attr['accession'], attr['experiment-type'], attr['instrument-platform'], attr['instrument-model']))
                    total += 1

        if x:
            sub_runs_by_biome[biome_name] = x
            biome_runs_counter[biome_name] = len(runlist)
            biome_illumina_counter[biome_name] = len(x)
            total_illumina += len(x)

    print(f'found {total_illumina} total w/Illumina, of {total}')

    for biome_name, ill_run_count in biome_illumina_counter.most_common():
        run_count = biome_runs_counter[biome_name]
        star = '*' if run_count > 200 else ''
        print(f"{star}{ill_run_count} of {run_count} - {biome_name}")

    if args.output_spreadsheet:
        print(f"outputting spreadsheet to '{args.output_spreadsheet}'")

        sig_path_out = None
        if args.save_sig_paths:
            sig_path_out = open(args.save_sig_paths, 'wt')

        with open(args.output_spreadsheet, 'w', newline='') as fp:
            w = csv.writer(fp)

            w.writerow(['biome1', 'biome2', 'biome3',
                       'accession', 'experiment_type',
                        'platform', 'model', 'sigpath'])

            n_missing = 0
            total_count = 0
            for biome_name, run_tuples in sub_runs_by_biome.items():
                for rt in run_tuples:
                    total_count += 1
                    biome_name, acc, exptype, platform, model = rt
                    sigpath = WORT_PATH.format(acc=acc)
                    if not os.path.exists(sigpath):
                        sigpath = ''
                        n_missing += 1
                    else:
                        if sig_path_out:
                            sig_path_out.write(sigpath + "\n")

                    w.writerow([truncate_biome(biome_name, 2),
                                truncate_biome(biome_name, 3),
                                truncate_biome(biome_name, 4),
                                acc,
                                exptype,
                                platform,
                                model,
                                sigpath])

            print(f"of {total_count} rows, {n_missing} missing sig files.")
            if sig_path_out:
                print(f"wrote {total_count - n_missing} sig paths to '{args.save_sig_paths}'")

    sys.exit(0)


if __name__ == '__main__':
    sys.exit(main())
