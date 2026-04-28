#! /usr/bin/env python
# OG location: https://github.com/ctb/2024-summarize-assembly-mapping/blob/main/summarize-ref-assembly.py
import os.path
import json
import csv
import sys
import argparse
import yaml
import pprint
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

import pandas
import sourmash
from sourmash import sourmash_args
from sourmash import minhash
from itertools import chain, combinations


# https://docs.python.org/3/library/itertools.html
def powerset(iterable, *, start=2):
    "powerset([1,2,3]) → () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(start, len(s)+1))

def isect_all_sigs(siglist):
    get_name = lambda x: [ss.name for ss in x]
    
    # check: all names are distinct
    assert len(set(get_name(siglist))) == len(siglist)

    # make a powerset, all combinations
    orig_pset = [ list(x) for x in powerset(siglist, start=1) ]
    orig_names = [ get_name(x) for x in orig_pset ]

    # sort most combinations to fewest
    pset = sorted(orig_pset, key=lambda x: -len(x))

    # walk through, calculating intersections and tracking them.
    # no hash should appear in more than one set.
    subtract_me = set()
    isects = []
    names_to_hashes = {}
    for n, combo in enumerate(pset):
        combo = list(combo)
        combo_name = get_name(combo)
        ss = combo.pop()
        hashes = set(ss.minhash.hashes) - subtract_me

        while combo and hashes:
            ss = combo.pop()
            hashes.intersection_update(ss.minhash.hashes)

        # track hashes by combination name so we can re-order later.
        names_to_hashes[tuple(combo_name)] = hashes
        subtract_me.update(hashes)

    # use leftover ss to get a minhash :)
    blank_mh = ss.minhash.copy_and_clear()

    # return list in order of original powerset;
    # convert hashes to minhash objects
    combos = []
    names = []
    for name in orig_names:
        hashes = names_to_hashes[tuple(name)]
        names.append(name)
        mh = blank_mh.copy_and_clear().flatten()
        mh.add_many(hashes)
        combos.append(mh)

    return names, combos


# add:
# * size of metagenome (metag_weighted)
# we have:
# * assembly_refmap_isect_w
# * assembly_f_readmapped_w
# * assembly_f_weighted

class MetagenomeInfo:
    headers = ["accession", "assembly_f_unweighted", "assembly_f_weighted",
               "assembly_f_readmapped", "assembly_f_readmapped_w",
               "ref_f_unweighted", "ref_f_weighted", "f_reads_mapped",
               "assembly_refmap_isect_w",
               "yaml_n_bases", "yaml_n_reads", "yaml_kmers",
               "yaml_known_hashes", "yaml_unknown_hashes", "yaml_total_hashes",
               'w_isect_all',
               'w_isect_matches',
               'w_no_isect',
               'w_isect_mapassem',
               'w_mapassem_only',
               'w_map_only']
               
    def __init__(self, metag_acc, *, ksize=31, grist_dir, assembly_dir):
        self.metag_acc = metag_acc
        self.ksize = ksize
        self.grist_dir = grist_dir
        self.assembly_dir = assembly_dir

        self.metag_sig_path = os.path.join(grist_dir, "sigs",
                                           f"{self.metag_acc}.trim.sig.zip")
        self.metag_sig = sourmash_args.load_one_signature(self.metag_sig_path, ksize=self.ksize)
        self.gather_matches_mh = None
    

    def calc(self):
        "Calculate all the numbers."
        print(f'for accession {self.metag_acc}:')
        self.load_gather_matches(self.grist_dir)
        self.calc_assembly_stuff(self.assembly_dir)
        self.calc_ref_based_kmer_stuff(self.grist_dir)
        self.calc_mapping_stuff(self.grist_dir)

    def get_row(self):
        """Return a CSV row, starting with metagenome accession and
        containing the rest of the attributes in self.headers."""
        xx = [self.metag_acc]
        for x in self.headers[1:]:
            val = getattr(self, x)
            xx.append(val)
        return xx

    def calc_assembly_stuff(self, assembly_dir):
        sigfile = os.path.join(assembly_dir,
                               f"{self.metag_acc}.megahit.fa.gz.sig")
        assert os.path.exists(sigfile), sigfile

        self.assembly_sig = sourmash_args.load_one_signature(sigfile,
                                                             ksize=self.ksize)

        # percent of flat k-mers accounted for by assembly
        self.assembly_f_unweighted = self.metag_sig.contained_by(self.assembly_sig)
        print(f"assembly/unweighted: {self.assembly_f_unweighted*100:.1f}% (assembly_f_unweighted)")

        # abundance weighted version:
        assembly_mh = self.assembly_sig.minhash.flatten()
        metag_mh = self.metag_sig.minhash
        intersect = assembly_mh.intersection(metag_mh.flatten()).inflate(metag_mh)

        # now sum:
        total_weighted_sum = metag_mh.sum_abundances
        intersect_weighted_sum = intersect.sum_abundances
        print(f"assembly/weighted: {intersect_weighted_sum / total_weighted_sum * 100:.1f}% (assembly_f_weighted)")
        self.assembly_f_weighted = intersect_weighted_sum / total_weighted_sum

        sig_mapped = os.path.join(assembly_dir,
                                  f'{self.metag_acc}.x.ma.fq.gz.sig')
        ma_sig = sourmash_args.load_one_signature(sig_mapped, ksize=self.ksize)
        self.assembly_f_readmapped = self.metag_sig.contained_by(ma_sig)
        print(f"% k-mers in reads mapped to assembly (flat): {self.assembly_f_readmapped*100:.1f}% (assembly_f_readmapped)")

        # calculate weighted version.
        self.assembly_f_readmapped_w = self.metag_sig.contained_by_weighted(ma_sig)
        print(f"% k-mers in reads mapped to assembly (weighted): {self.assembly_f_readmapped_w*100:.1f}% (assembly_f_readmapped_w)")

        gather_sig = sourmash.SourmashSignature(self.gather_matches_mh,
                                                name=self.metag_sig.name+'.gather')
        siglist = [self.metag_sig, ma_sig, self.assembly_sig, gather_sig]

        names, isects = isect_all_sigs(siglist)

        sizes = [ len(x) for x in isects ]
        #pprint.pprint(list(zip(range(0, len(names)), names, sizes)))

        # restore:
        metag_weighted_mh = self.metag_sig.minhash
        isects = [ mh.inflate(metag_weighted_mh) for mh in isects ]
        sizes = [ mh.sum_abundances for mh in isects ]

        assert sum(sizes) == metag_weighted_mh.sum_abundances, (sum(sizes), metag_weighted_mh.sum_abundances)

        # [(0, ['DRR014782'], 300177),
        self.w_no_isect = sizes[0]
        # (1, ['DRR014782.ma'], 0),
        # sizes[1]: map assem only, always zero
        # (2, ['DRR014782.assembly'], 766),
        # sizes[2]: novel assem k-mers only, close to 0
        # (3, ['DRR014782.gather'], 1317562),
        # sizes[3]: stuff from ref genomes that is not in metagenome, ignore
        # (4, ['DRR014782', 'DRR014782.ma'], 107373),
        self.w_mapassem_only = sizes[4]
        # (5, ['DRR014782', 'DRR014782.assembly'], 3),
        self.w_assem_only = sizes[5] # metag+assembly only, close to zero
        # (6, ['DRR014782', 'DRR014782.gather'], 250384),
        self.w_map_only = sizes[6] # metag+gather only
        # (7, ['DRR014782.ma', 'DRR014782.assembly'], 0),
        # sizes[7]: ma and assembly only, zero
        # (8, ['DRR014782.ma', 'DRR014782.gather'], 0),
        # sizes[8]: ma and gather only, zero
        # (9, ['DRR014782.assembly', 'DRR014782.gather'], 947),
        # sizes[9]: assembly and gather only, close to zero
        # (10, ['DRR014782', 'DRR014782.ma', 'DRR014782.assembly'], 34747),
        # sizes[10]: metag + map assem + assembly
        self.w_isect_mapassem = sizes[10] # metag + mapassem + assembly
        # (11, ['DRR014782', 'DRR014782.ma', 'DRR014782.gather'], 33249),
        self.w_isect_matches = sizes[11] # metag + mapassem + gather
        # (12, ['DRR014782', 'DRR014782.assembly', 'DRR014782.gather'], 5),
        # sizes[12]: metag, assembly, gather, close to 0
        # (13, ['DRR014782.ma', 'DRR014782.assembly', 'DRR014782.gather'], 0),
        # sizes[13]: map assem, assembly, gather, no metagenome, 0
        # (14,
        # ['DRR014782', 'DRR014782.ma', 'DRR014782.assembly', 'DRR014782.gather'],
        # 229312)]

        self.w_isect_all = sizes[14]

    def calc_ref_based_kmer_stuff(self, grist_dir):
        gather_csv = os.path.join(grist_dir, "gather", f"{self.metag_acc}.gather.csv.gz")
        df = pandas.read_csv(gather_csv)

        self.ref_f_unweighted = df['f_unique_to_query'].sum()
        self.ref_f_weighted = df['f_unique_weighted'].sum()
        
        print(f"total ref k-mers found (abund): {self.ref_f_weighted * 100:.1f} (ref_f_weighted)")
        print(f"total ref k-mers found (flat): {self.ref_f_unweighted * 100:.1f} (ref_f_unweighted)")

        assembly_mh = self.assembly_sig.minhash
        isect_mh = minhash.flatten_and_intersect_scaled(self.gather_matches_mh,
                                                        assembly_mh)
        isect_weighted = self.metag_sig.minhash.contained_by_weighted(isect_mh)
        self.assembly_refmap_isect_w = isect_weighted
        print(f"overlap between ref/assembly, weighted {self.assembly_refmap_isect_w*100:.1f}% (assembly_refmap_isect_w)")

    def calc_mapping_stuff(self, grist_dir):
        leftover_csv = os.path.join(grist_dir, 'leftover',
                                    f"{self.metag_acc}.summary.csv")
        df = pandas.read_csv(leftover_csv)
        total_mapped_reads = df['n_mapped_reads'].sum()
        print(f"total mapped reads: {total_mapped_reads}")

        read_stats_file = os.path.join(grist_dir, 'trim', f"{self.metag_acc}.trim.json")
        with open(read_stats_file, 'rb') as fp:
            read_stats = json.load(fp)
        total_reads = read_stats['summary']['after_filtering']['total_reads']
        f_mapped = total_mapped_reads / total_reads
        print(f"fraction of reads that mapped: {f_mapped*100:.1f}% (f_mapped)")
        self.f_reads_mapped = f_mapped

        summary_file = os.path.join(grist_dir, f'{self.metag_acc}.info.yaml')
        with open(summary_file, 'rb') as fp:
            summary_d = yaml.load(fp, Loader=Loader)
        for key in summary_d:
            attr = f'yaml_{key}'
            setattr(self, attr, summary_d[key])

    def load_gather_matches(self, grist_dir):
        matches_zip = os.path.join(grist_dir,
                                   f'gather/{self.metag_acc}.matches.sig.zip')
        matches_mh = None
        for ss in sourmash.load_file_as_signatures(matches_zip):
            if matches_mh is None:
                matches_mh = ss.minhash.copy_and_clear()

            matches_mh += ss.minhash

        self.gather_matches_mh = matches_mh


def main(argv):
    p = argparse.ArgumentParser(argv)
    p.add_argument('accs', nargs='+')
    p.add_argument('-o', '--output-csv', help='output CSV here',
                   required=True)
    p.add_argument('-t', '--top-level-directory', required=True,
                   help="directory containing assembly/, atta/, and grist/outputs/, e.g. '/group/ctbrowngrp4/2025-zyzhao-assemloss/85/'")
    
    args = p.parse_args()

    tld = args.top_level_directory.rstrip('/') + '/'
    grist_dir = os.path.join(tld, 'grist/outputs')
    assembly_dir = os.path.join(tld, 'assembly/')

    # check all the paths
    for path in (tld, grist_dir, assembly_dir):
        assert os.path.isdir(path), f'{path} is not a directory or does not exist?'

    # gather info for each accession
    results = []
    for metag in args.accs:
        info = MetagenomeInfo(metag,
                              grist_dir=grist_dir,
                              assembly_dir=assembly_dir)
        info.calc()
        results.append(info.get_row())

    # write out CSV
    with open(args.output_csv, 'w', newline='') as fp:
        w = csv.writer(fp)
        w.writerow(MetagenomeInfo.headers)
        for rr in results:
            w.writerow(rr)


if __name__ == '__main__':
    sys.exit(main(sys.argv))
