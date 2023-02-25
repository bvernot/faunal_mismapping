#!/usr/bin/env python3

import argparse
import functools
import math
import random
import signal
import sys
import os.path
from collections import Counter
from collections import defaultdict
from collections import namedtuple

import pysam
import pandas as pd
import numpy as np

signal.signal(signal.SIGPIPE, signal.SIG_DFL)


def majority(bases, prop):
    """Select the most frequent allele, requiring at least `prop` of all
    reads to agree on it.

    Cases with two alleles at the same proportion are handled by
    implicit behaviour of the `most_common()` method, which assigns
    the order of frequency arbitrarily.
    """
    counts = Counter(bases).most_common()
    allele, n =  counts[0]

    # return the most frequent allele in case there are at most two
    # alleles present at a site
    if n / len(bases) >= prop and len(counts) <= 2:
        return allele
    else:
        return "N"


def flush(i, calls, out_fun, times_called):
    """Save genotype calls to a file (either a VCF or a pileup file)."""
    print(f"\r{i + 1} positions processed", end="", file=sys.stderr)
    calls = pd.DataFrame(
        calls, columns=["chrom", "pos", "ref", "coverage", "call"]
    ).query('ref != "N" & call != "N"')
    out_fun(calls, times_called)
    return times_called + 1


def get_ref_base(col):
    tuples = col.pileups[0].alignment.get_aligned_pairs(with_seq=True)
    for read_pos, ref_pos, ref_base in tuples:
        if ref_pos == col.reference_pos:
            return ref_base.upper()
        pass
    pass


site_results = dict()
config_results = dict()

def process_site(chrom, pos, sites, ref_bam, bases, third_sites, read_names):

    # bases = [b for b in bases] # if b in (a1, a2)]

    if len(bases) > 1:
        print('MORE BASES')
        pass

    a1, a2, category = sites[(chrom, pos)]
    ref_control, _, _, ref_human, a3, *flag = third_sites[(chrom,pos)]

    anc_base = [a1, a2]
    anc_base.remove(ref_human)
    if len(anc_base) != 1:
        print('WERID ERROR')
        pass
    anc_base = anc_base[0]

    der_base = ref_human

    a4 = [b for b in ['A', 'C', 'G', 'T'] if b not in [anc_base, der_base, a3]]
    if len(a4) != 1:
        print('WERID ERROR')
        pass
    a4 = a4[0]
    
    # print(bases, ref_human, a2, a3)
    # print(anc_base, der_base, a3, a4)

    config = (anc_base, der_base, a3, a4)
    if config not in config_results:
        config_results[config] = {'n_anc' : 0,
                                  'n_der' : 0,
                                  'n_a3' : 0,
                                  'n_a4' : 0}
        pass

    config_results[config]['n_anc'] += sum(b == anc_base for b in bases)
    config_results[config]['n_der'] += sum(b == der_base for b in bases)
    config_results[config]['n_a3'] += sum(b == a3 for b in bases)
    config_results[config]['n_a4'] += sum(b == a4 for b in bases)

    if len(bases) == len(read_names):
        for i in range(len(bases)):
            print('READ_NAMES', ['anc_base', 'der_base', 'a3', 'a4'][[anc_base, der_base, a3, a4].index(bases[i])], read_names[i])
            pass
        pass
    else:
        print('READ_NAMES', 'LENGTH_MISMATCH')
        pass

    # print(config_results)
    
    ##############
    ### SOMETHING TO SPLIT OUT Ti Tv
    ##############

    # print(ref, a1, a2, a3, ref_fasta)
    # if ref not in (a1, a2):
    #     print('Skipping site because ref not in a1/a2!', ref, a1, a2)
    #     return

    if category not in site_results:
        site_results[category] = {'n_sites_a12' : 0,
                                 'n_sites' : 0,
                                 'n_der' : 0,
                                 'n_anc' : 0,
                                 'ref_bam' : 0,
                                  # 'ref_nonmatch' : 0,
                                 'n_first' : 0,
                                 'n_second' : 0, 
                                 'n_third' : 0,
                                 'n_fourth' : 0,
                                 'ref_bam_not_a3' : 0}
        pass


    # ref_match = sum(b == ref for b in bases if b in (a1, a2))
    # ref_nonmatch = sum(b != ref for b in bases if b in (a1, a2))
    # third_nonmatch =  sum(b != ref for b in bases if b not in (a1, a2))
    ref_match = sum(b == ref_human for b in bases if b in (a1, a2)) #ref from bam
    ref_nonmatch = sum(b != ref_human for b in bases if b in (a1, a2))
    first_allele = sum(b == a1 for b in bases)# if b in (a1, a2))
    second_allele = sum(b == a2 for b in bases)# if b in (a1, a2))
    third_allele = sum(b == a3 for b in bases)# if b not in (a1, a2))
    fourth_allele = sum(b == b for b in bases if b not in (a1, a2, a3))
    ref_bam_match = sum(b == ref_bam for b in bases)
    # print(bases, ref, ref_match, ref_nonmatch)

    site_results[category]["n_sites_a12"] += ref_match + ref_nonmatch
    site_results[category]["n_der"] += ref_match
    site_results[category]["n_anc"] += ref_nonmatch

    site_results[category]["n_sites"] += ref_match + ref_nonmatch + third_allele + fourth_allele

    site_results[category]["n_first"] += first_allele
    site_results[category]["n_second"] += second_allele
    site_results[category]["n_third"] += third_allele
    site_results[category]["n_fourth"] += fourth_allele

    # site_results[category]["ref_nonmatch"] += ref_nonmatch
    site_results[category]["ref_bam"] += ref_bam_match
    site_results[category]["ref_bam_not_a3"] += ref_bam != a3

    # print(site_results)
    
    return
    

def call_bases(call_fun, out_fun, bam, mincov, minbq, minmq, minlen, chrom, sites, flush_interval, limit, random_bases, third, subs):
    """Sample bases in a given region of the genome based on the pileup
    of reads. If no coordinates were specified, sample from the whole BAM file.
    """
    calls = []
    i = 0
    last_flush = 0
    times_called = 0

    # subs = {}
    read_names = dict()
    
    for i, col in enumerate(bam.pileup(contig=chrom, compute_baq=False,
                                       min_base_quality=minbq,
                                       min_mapping_quality=minmq)):

        if limit is not None and i > limit: break

        # print(col.reference_name, col.reference_pos)
        # print(sites)


        debug_reads = True
        debug_reads = False

        
        
        for fragment_idx, pileup in enumerate(col.pileups):
            
            if pileup.alignment.query_name in read_names:
                continue

            if len(pileup.alignment.query_sequence) < minlen:
                if debug_reads: print('Filtering read due to minimum length', pileup.alignment.query_name, len(pileup.alignment.query_sequence))
                continue

            if debug_reads: print()
            
            ## save this read name so we don't process it again
            read_names[pileup.alignment.query_name] = 0
            
            ## use these reads to build up a substitution model
            if debug_reads: print(pileup.alignment.query_name)
            bases = col.get_query_sequences(add_indels=True)
            if debug_reads: print(bases)
            if debug_reads: print(col.get_query_positions())
            
            if debug_reads: print(pileup.alignment.get_aligned_pairs(with_seq=True, matches_only=False))
            if debug_reads: print(pileup.alignment.query_sequence, '         <-- query sequence (no gaps)')
            # if debug_reads: print(pileup.alignment.get_reference_sequence()) ## includes clipped bases?
            # if debug_reads: print(''.join(str(x) if x[2] is not None else '-' for x in pileup.alignment.get_aligned_pairs(with_seq=True, matches_only=True)))
            # if debug_reads: print(''.join(x[1] if x[1] is not None else '-' for x in pileup.alignment.get_aligned_pairs(with_seq=True, matches_only=True)))
            q_seq = ''.join(pileup.alignment.query_sequence[x[0]] if x[0] is not None else '-' for x in pileup.alignment.get_aligned_pairs(with_seq=True, matches_only=False))
            r_seq = ''.join(x[2] if x[1] is not None else '-' for x in pileup.alignment.get_aligned_pairs(with_seq=True, matches_only=False))
            q_pos = [x[0] for x in pileup.alignment.get_aligned_pairs(with_seq=True, matches_only=False)]
            r_pos = [x[1] for x in pileup.alignment.get_aligned_pairs(with_seq=True, matches_only=False)]
            
            if debug_reads: print(q_seq, '<-- query sequence (with gaps)')
            if debug_reads: print(r_seq, '<-- reference sequence (with gaps)')
            
            for pos in range(len(q_seq)):
                if r_pos[pos] is not None \
                   and q_pos[pos] is not None \
                   and (col.reference_name, r_pos[pos]+1) not in sites:
                    s = (r_seq[pos].upper(), q_seq[pos].upper())
                    if s not in subs: subs[s] = 0
                    subs[s] += 1
                    pass
                
                if debug_reads: print(q_seq)
                if debug_reads: print(r_seq)
                if debug_reads: print(' ' * pos + q_seq[pos], '<-- SNP' if r_pos[pos] is not None and (col.reference_name, r_pos[pos]+1) in sites else '')
                if debug_reads: print(' ' * pos + r_seq[pos])
                pass
            # if debug_reads: print(' ' * col.get_query_positions()[fragment_idx] + pileup.alignment.query_sequence[col.get_query_positions()[fragment_idx]])
            # if debug_reads: print(' ' * col.get_query_positions()[fragment_idx] + ''.join(x[2] if x[1] is not None else '-' for x in pileup.alignment.get_aligned_pairs(with_seq=True, matches_only=False))[col.get_query_positions()[fragment_idx]])
            
            ## print the SNP site(s) if any exist in this read
            s_pos = [p for p in r_pos if p is not None and (col.reference_name, p+1) in sites]
            if debug_reads: print(''.join(q_seq[i]+'<- SNP' if p is not None and (col.reference_name, p+1) in sites else ' ' for i,p in enumerate(r_pos)))
            if debug_reads: print(''.join(' ' if x[1] is not None and ((col.reference_name, x[1]+1,) not in sites) else str(x[2]) + ' HEY' for x in pileup.alignment.get_aligned_pairs(with_seq=True, matches_only=False)))
            if debug_reads: print('SNP pos:', s_pos)

            if debug_reads: print('substitution matrix:', subs)
            pass

        ###
        ## only process the sites that we requested
        if sites is not None and (col.reference_name, col.reference_pos + 1) not in sites:
            continue

        
        bases = col.get_query_sequences(add_indels=True)

        # filter out sites with no reads and sites with indels
        if bases and "*" not in bases and all(len(i) == 1 for i in bases):
            bases = [b.upper() for b in col.get_query_sequences()
                     if b in "ACGTacgt"]

            #print('~~~~~~~~~~')
            #print(bases)
            if minlen is not None:
                read_lens = [len(pileupread.alignment.get_reference_positions(full_length=True)) for pileupread in col.pileups]
                bases = [bases[i] for i,l in enumerate(read_lens) if l >= minlen]
                pass
            #print(read_lens)
            #print(bases)

            #### try to look at deamination also?
            # print(col.pileups[0].alignment.query_sequence)
            # print(''.join(x[2] if x[2] is not None else '-' for x in col.pileups[0].alignment.get_aligned_pairs(with_seq=True)))
            
            if random_bases:
                random.seed(1984)
                bases = random.choices(bases, k = math.floor(len(bases) / 2))
                pass

            if len(bases) >= mincov:

                ## for calculating stats
                process_site(col.reference_name,
                             col.reference_pos + 1,
                             sites,
                             get_ref_base(col),
                             bases, third,
                             [pileup.alignment.query_name for pileup in col.pileups])
                
                ## for outputing a file with every site listed
                calls.append((
                    col.reference_name,
                    col.reference_pos + 1,
                    get_ref_base(col),
                    len(bases),
                    call_fun(bases)
                ))
                pass
            pass

        if i-last_flush > flush_interval:
            times_called = flush(i, calls, out_fun, times_called)
            calls = []
            last_flush = i
            pass

        pass
    
    flush(i, calls, out_fun, times_called)
    return


def write_pileup(pileups, times_called, output, new_file, records_dict, report = True):

    if not report: return
    
    if new_file and times_called == 0: print("SITES\tchrom\tpos\tref\tpileup\tA\tC\tG\tT\tmatch", file=output)
    for i in pileups.itertuples():
        counts = Counter(i.call)
        counts_str = '\t'.join(str(counts[i]) for i in 'ACGT')
        if i.ref == i.call:
            match = 'ref_match'
        else:
            match = 'ref_mismatch'
            pass
        records_dict[match] += 1
        print(f"SITES\t{i.chrom}\t{i.pos}\t{i.ref}\t{''.join(i.call)}\t{counts_str}\t{match}", file=output)


def check_range(value):
    """Make sure that the required proportion of reads in agreement on a
    base is between 0.5 and 1.0.
    """
    value = float(value)
    if value < 0.5 or value > 1.0:
        raise argparse.ArgumentTypeError(
            "Required proportion for majority calling "
            f"needs to be between 0.5 and 1.0 (value given: {value}).")
    return value


def read_sites(sites_file, add_chr):
    sites = dict()
    with open(sites_file, "rt") as sf:
        for i, line in enumerate(sf):
            line = line.rstrip().split()
            chrom,pos = line[0:2]
            if i == 0 and not pos.isnumeric():
                print('Keeping header:', line, file=sys.stderr)
                sites['category_header'] = tuple(line[4:])
                continue
            
            a1_der, a2_anc = line[2:4]

            category = tuple(line[4:])
            
            if add_chr:
                chrom = 'chr' + chrom
                pass
            pos = int(pos)
            sites[(chrom,pos)] = (a1_der, a2_anc, category)
            pass
        pass
    return sites
    
def read_third_file(third_file):
    third_sites = dict()
    with open(third_file, "r") as tf:
        for i, line in enumerate(tf):
            line = line.rstrip().split()
            chrom, pos, ref_control, a1, a2, ref_fasta, a3, *flag = line
            pos = int(pos)
            third_sites[(chrom,pos)] = (ref_control, a1, a2, ref_fasta, a3, *flag)
    return third_sites

def check_dicts(third, sites):
    sites_no_header = dict(sites)
    del sites_no_header['category_header']
    set1 = set(third.keys())
    set2 = set(sites_no_header.keys())
    merged = set1 ^ set2
    assert merged == set()

    for key in third: #key == chrom, pos
        _, a1, a2, _, _, *_ = third[key]
        assert a1 == sites[key][0]
        assert a2 == sites[key][1]
        pass
    pass



def report_stats(args, sites, site_results):

    all_cats = sorted(site_results.keys())

    all_stats_keys = site_results[all_cats[0]].keys()
    print(all_stats_keys)

    sum_stats = {k: sum(d[k] for d in site_results.values() if k in d) for k in all_stats_keys}

    print(sum_stats)

    if args.by_cats:
        # if 'category_header' in sites:
        print("\n")
        print('CATSITES',
              '\t'.join(sites['category_header']),
              '\t'.join(all_stats_keys),
              sep='\t')

        for category in all_cats:
            print('CATSITES',
                  '\t'.join(category), 
                  '\t'.join(site_results[category][k] for k in all_stats_keys),
                  sep='\t')
            pass
        pass
    

    print('\n')
    print('ALLSITES', '\t'.join(all_stats_keys), sep='\t')
    print('ALLSITES',
          '\t'.join(str(sum_stats[k]) for k in all_stats_keys),
          sep='\t')

    pass



def calc_lik_for_params(all_cats, config_results, f_anc, f_der, f_a4, subs):

    lik = 0
    for config in all_cats:

        # print(config)
        config_stats = config_results[config]
        # print(config_stats)

        ## observation counts
        n_anc, n_der, n_a3, n_a4 = [config_stats[c] for c in ('n_anc', 'n_der', 'n_a3', 'n_a4')]
        # print(n_anc, n_der, n_a3, n_a4)

        ## alleles for each state
        anc, der, a3, a4 = config

        ## loop through observations
        bases = ['A', 'C', 'G', 'T']
        for obs in bases:

            ## "loop" through hidden state genotypes
            temp_lik = 0
            
            ## anc
            temp_lik += (f_anc + f_a4) * subs[(anc,obs)]
            
            ## der
            temp_lik += (f_der + f_a4) * subs[(der,obs)]
            
            ## a3
            temp_lik += (1 - f_anc - f_der - 3*f_a4) * subs[(a3,obs)]
            
            ## a4
            temp_lik += (f_a4) * subs[(a4,obs)]

            i = config.index(obs)
            n_obs = [n_anc, n_der, n_a3, n_a4][i]
            # print(i, ':', n_obs)
            
            lik += n_obs * math.log(temp_lik)

            pass

        pass

    return lik


def calc_lik(args, sites, config_results, subs):

    all_cats = sorted(config_results.keys())

    all_stats_keys = config_results[all_cats[0]].keys()

    print('lik', all_cats)
    print('lik', all_stats_keys)

    sum_stats = {k: sum(d[k] for d in config_results.values() if k in d) for k in all_stats_keys}

    print('lik', sum_stats)

    f_anc, f_der, f_a4 = (.1, .5, .05)

    bases = ['A', 'C', 'G', 'T']

    #subs = {(h,o) : 0.997 if h == o else 0.001 for h in bases for o in bases}
    print(subs)

    all_liks = []
    all_params = []

    print('LIKELIHOOD', 'f_anc', 'f_der', 'f_a4', 'lik', 'rel_f_anc', 'rel_f_der', sep='\t')

    for f_anc in np.arange(0,1.1,.01):
        for f_der in np.arange(0,1.1,.01):
            for f_a4 in np.arange(0,.3,.01):

                if f_anc + f_der + 3*f_a4 > 1:
                    continue
                
                lik = calc_lik_for_params(all_cats, config_results, f_anc, f_der, f_a4, subs)
                print('LIKELIHOOD', f_anc, f_der, f_a4, lik, f_anc/(f_anc+f_der), f_der/(f_anc+f_der), sep='\t')
                all_params += [(f_anc, f_der, f_a4)]
                all_liks += [lik]
                
                pass
            pass
        pass

    max_lik_i = np.argmax(all_liks)
    f_anc, f_der, f_a4 = all_params[max_lik_i]
    max_lik = all_liks[max_lik_i]
    print('MAXIMUM_LIKELIHOOD', f_anc, f_der, f_a4, max_lik, f_anc/(f_anc+f_der), f_der/(f_anc+f_der), sep='\t')

    return(max_lik, max_lik_i, all_liks, all_params)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Call alleles from a BAM file using various criteria")
    parser.add_argument("--bam", help="BAM file to sample from", required=True)
    parser.add_argument("--chrom", help="Chromosome to sample from")
    parser.add_argument("--strategy", help="How to 'genotype'?", choices=["random", "majority", "pileup", "none"], required=True)
    parser.add_argument("--proportion", help="Required proportion of the majority allele", type=check_range, default=0.5)
    parser.add_argument("--seed", help="Set seed for random allele sampling [random]")
    parser.add_argument("--mincov", help="Minimum coverage", type=int, default=1)
    parser.add_argument("--minbq", help="Minimum base quality", type=int, default=13)
    parser.add_argument("--minmq", help="Minimum read mapping quality", type=int, default=0)
    parser.add_argument("--minlen", help="Minimum read length", type=int, default=35)
    parser.add_argument("--add-chr", help="Add 'chr' to the beginining of every site's chromosome, to match bam file.", action='store_true')
    parser.add_argument("--sample-name", help="Sample name to put in a VCF header")
    parser.add_argument("--output", help="Output file name")
    parser.add_argument("--sites", help="Restrict output to this list of sites (chrom, 1 based pos)")
    parser.add_argument("--flush", help="Print to file every N bases", type=int, default=100000)
    parser.add_argument("--limit", help="Stop after processing N bases from bam", type=int, default=None)
    parser.add_argument("--third_file", help="Defines which nucleotide is the third allele")
    parser.add_argument("--random_bases", help="Restricts analysis to a random subset of bases.", type=bool)
    parser.add_argument("--by-cats", help="Print statistics for each SNP category", action=argparse.BooleanOptionalAction)


    args = parser.parse_args()

    # if args.strategy != "pileup" and not args.sample_name:
    #     parser.error(f"Sample name has to be specified when writing a VCF file")
    #     pass

    if args.output is None:
        new_file = True
        output = sys.stdout
    else:
        new_file = not os.path.isfile(args.output)
        output = open(args.output, "w" if new_file else "a")
        pass
    records_dict = defaultdict(int)

    bam = pysam.AlignmentFile(args.bam)

    if not bam.has_index():
        print("BAM file index is missing", file=sys.stderr)
        sys.exit(1)
        pass

    if args.strategy == "pileup":
        call_fun = lambda x: "".join(x)
        out_fun = functools.partial(write_pileup, output=output, new_file=new_file, records_dict=records_dict)

    elif args.strategy == "random":
        if args.seed:
            random.seed(args.seed)
        call_fun = lambda x: random.choice(x)
        #out_fun = functools.partial(write_vcf, output=args.output,
        #                            sample_name=args.sample_name)
        out_fun = functools.partial(write_pileup, output=output, new_file=new_file, records_dict=records_dict, report=False)

    elif args.strategy == "none":
        call_fun = lambda x: "".join(x)
        out_fun = functools.partial(write_pileup, output=output, new_file=new_file, records_dict=records_dict, report=False)
        pass

    # elif args.strategy == "majority":
    #     call_fun = functools.partial(majority, prop=args.proportion)
    #     out_fun = functools.partial(write_vcf, output=args.output,
    #                                 sample_name=args.sample_name)

    if args.third_file:
        third = read_third_file(args.third_file)


    if args.sites is not None:
        sites = read_sites(args.sites, args.add_chr)
    else:
        sites = None
        pass

    check_dicts(third, sites)

    subs_matrix = dict()
    call_bases(call_fun, out_fun, bam, args.mincov,
               args.minbq, args.minmq, args.minlen,
               args.chrom, sites, args.flush, args.limit,
               args.random_bases, third, subs = subs_matrix)

    subs_probs = dict()
    print('Substitution matrix + probabilities')
    print(subs_matrix)
    for b1 in list('ACGT'):
        s = sum(subs_matrix[(b1,b2)] for b2 in list('ACGT'))
        for b2 in list('ACGT'):
            subs_probs[(b1,b2)] = subs_matrix[(b1,b2)] / s 
            print('', b1, b2, subs_matrix[(b1,b2)], subs_probs[(b1,b2)], sep='\t')
            pass
        pass
            
    
    # print(records_dict, file=sys.stderr)

    # print(site_results)


    report_stats(args, sites, site_results)
    calc_lik(args, sites, config_results, subs_probs)

    exit()
    
    
    all_cats = sorted(site_results.keys())
    print(all_cats)
    all_cats = site_results.keys()
    print(all_cats)


    nsites = 0
    nder = 0
    ref_bam = 0
    non_ref = 0

    n_first = 0
    n_second = 0
    n_third = 0
    n_fourth = 0

    ref_bam_not_a3 = 0

    if args.by_cats:
        if 'category_header' in sites:
            print("\n")
            print('CATSITES', '\t'.join(sites['category_header']), 'nder', 'nsites',
            'ref_bam', 'non_ref', 'n_1st', 'n_2nd', 'n_3rd', 'n_4th', 'ref_bam_not_a3', sep='\t')


    all_cats = sorted(site_results.keys())
    #    all_cats = sorted(site_results.keys(), key = lambda x : tuple(x[0], int(x[1]), x[2] if int(x[2]) < 14 else '14+'))

    for category in all_cats:

        nsites += site_results[category]['n_sites']
        nder += site_results[category]['n_der']
        ref_bam += site_results[category]['ref_bam']
        non_ref += site_results[category]['ref_nonmatch']
        n_first += site_results[category]['n_first']
        n_second += site_results[category]['n_second']
        n_third += site_results[category]['n_third']
        n_fourth += site_results[category]['n_fourth']
        ref_bam_not_a3 += site_results[category]["ref_bam_not_a3"]

        if args.by_cats:
            print('CATSITES', '\t'.join(category), 
            site_results[category]['n_der'], 
            site_results[category]['n_sites'], 
            site_results[category]['ref_bam'],
            site_results[category]['ref_nonmatch'],
            site_results[category]['n_first'],
            site_results[category]['n_second'], 
            site_results[category]['n_third'], 
            site_results[category]['n_fourth'], 
            site_results[category]["ref_bam_not_a3"], sep='\t')

    print('\n')
    print('ALLSITES', 'nder', 'nsites','ref_bam', 'nonref', 'n_1st', 'n_2nd', 
    'n_3rd', 'n_4th', 'ref_bam_not_a3', sep='\t')

    print('ALLSITES', nder, nsites, ref_bam, non_ref, n_first, n_second, 
        n_third, n_fourth, ref_bam_not_a3, sep='\t')
    # print("Number of times that the reference of the bam file does not match the a3 allele given by the third_sites file:", site_results[category]['ref_bam_not_a3'],)

# def write_vcf(calls, output, sample_name):
#     filename = output + ".vcf"
#     new_file = not os.path.isfile(filename)
#     with open(filename, "w" if new_file else "a") as vcf:
#         if new_file:
#             print(
#                 "##fileformat=VCFv4.1\n"
#                 "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
#                 "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of high-quality bases\">\n"
#                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}".
#                 format(sample=sample_name), file=vcf
#             )
#         for i in calls.itertuples():
#             alt, gt = (".", 0) if i.ref == i.call else (i.call, 1)
#             print(f"{i.chrom}\t{i.pos}\t.\t{i.ref}\t{alt}\t.\t.\t.\tGT:DP\t{gt}:{i.coverage}", file=vcf)


# def write_pileup_og(pileups, output):
#     filename = output + ".txt"
#     new_file = not os.path.isfile(filename)
#     with open(filename, "w" if new_file else "a") as tsv:
#         if new_file: print("chrom\tpos\tref\tpileup\tA\tC\tG\tT", file=tsv)
#         for i in pileups.itertuples():
#             counts = Counter(i.call)
#             counts_str = '\t'.join(str(counts[i]) for i in 'ACGT')
#             print(f"{i.chrom}\t{i.pos}\t{i.ref}\t{''.join(i.call)}\t{counts_str}", file=tsv)
