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
# import pandas as pd
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


# def flush(i, calls, out_fun, times_called):
#     """Save genotype calls to a file (either a VCF or a pileup file)."""
#     print(f"\r{i + 1} positions processed", end="", file=sys.stderr)
#     calls = pd.DataFrame(
#         calls, columns=["chrom", "pos", "ref", "coverage", "call"]
#     ).query('ref != "N" & call != "N"')
#     out_fun(calls, times_called)
#     return times_called + 1


def get_ref_base(col):
    tuples = col.pileups[0].alignment.get_aligned_pairs(with_seq=True)
    for read_pos, ref_pos, ref_base in tuples:
        if ref_pos == col.reference_pos:
            return ref_base.upper()
        pass
    pass


#site_results = dict()
#config_results = dict()

def process_site(chrom, pos, sites, ref_bam, bases,
                 # third_sites,
                 results, read_names):

    # bases = [b for b in bases] # if b in (a1, a2)]

    if len(bases) > 1:
        if args.debug: print('Depth > 1 at position', chrom, pos)
        pass

    a1, a2, a3, ref_human, category = sites[(chrom, pos)]
    #print('sites[(chrom, pos)]', sites[(chrom, pos)])
    #ref_control, _, _, ref_human, a3, *flag = third_sites[(chrom,pos)]

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


    ####
    ## a "configuration" is some combination of anc, der, a3, a4
    ## i.e. A, T, C, G
    ## we keep stats separately for different configs because we use the substitution matrix (I think?)
    
    config = (anc_base, der_base, a3, a4)
    config_results = results['config']

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


    ####
    ## two cases where we need to loop over all reads overlapping this SNP
    ## 1. when reporting info for each read (for debugging, for example)
    ## 2. when checking simulated data to see what species it actually comes from

    if args.report_sim_truth is not None or args.report_reads:

        sim_results = results['sim_source']

        ## loop over all reads ('bases') overlapping this site
        for i in range(len(bases)):

            state = ['anc_base', 'der_base', 'a3', 'a4'][[anc_base, der_base, a3, a4].index(bases[i])]

            if args.report_reads:
                print('READ_NAMES', state, chrom, pos, bases[i], read_names[i], sep='\t')
                pass

            if args.report_sim_truth is not None:
                species = read_names[i].split('_')[0]
                if species not in sim_results:
                    sim_results[species] = dict()
                    sim_results[species]['N'] = 0
                    sim_results[species]['anc_base'] = 0
                    sim_results[species]['der_base'] = 0
                    sim_results[species]['a3'] = 0
                    sim_results[species]['a4'] = 0
                    pass
                sim_results[species][state] += 1
                sim_results[species]['N'] += 1
                # print(sim_results)
                pass
            
            pass
        pass


    site_results = results['sites']

    ####
    ## initialize site_results for this category (burden, etc)
    ## 

    # if False and category not in site_results:
    #     site_results[category] = {'n_reads_a12' : 0,
    #                              'n_reads' : 0,
    #                              'n_der' : 0,
    #                              'n_anc' : 0,
    #                              'ref_bam' : 0,
    #                               # 'ref_nonmatch' : 0,
    #                              'n_first' : 0,
    #                              'n_second' : 0, 
    #                              'n_third' : 0,
    #                              'n_fourth' : 0,
    #                              'ref_bam_not_a3' : 0}
    #     pass


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

    site_results[category]["n_reads_a12"] += ref_match + ref_nonmatch
    site_results[category]["n_der"] += ref_match
    site_results[category]["n_anc"] += ref_nonmatch

    site_results[category]["n_reads"] += ref_match + ref_nonmatch + third_allele + fourth_allele
    site_results[category]["n_snps"] += 1

    site_results[category]["n_first"] += first_allele
    site_results[category]["n_second"] += second_allele
    site_results[category]["n_third"] += third_allele
    site_results[category]["n_fourth"] += fourth_allele

    # site_results[category]["ref_nonmatch"] += ref_nonmatch
    site_results[category]["ref_bam"] += ref_bam_match
    site_results[category]["ref_bam_not_a3"] += ref_bam != a3

    # print(site_results)
    
    return


def compute_subs_counts(pileups, results, args, reference_name):

    subs = results['subs_matrix']
    read_names = results['read_names']
    debug_reads = True
    debug_reads = False
    
    for fragment_idx, pileup in enumerate(pileups):
        
        if pileup.alignment.query_name in read_names:
            continue
        
        if debug_reads: print()
        
        ## save this read name so we don't process it again
        # read_names[pileup.alignment.query_name] = 0
        
        ## use these reads to build up a substitution model
        if debug_reads: print(pileup.alignment.query_name)
        
        if debug_reads: print(pileup.alignment.get_aligned_pairs(with_seq=True, matches_only=False))
        if debug_reads: print(pileup.alignment.query_sequence, '         <-- query sequence (no gaps)')
        q_seq = ''.join(pileup.alignment.query_sequence[x[0]] if x[0] is not None else '-' for x in pileup.alignment.get_aligned_pairs(with_seq=True, matches_only=False))
        r_seq = ''.join(x[2] if x[1] is not None else '-' for x in pileup.alignment.get_aligned_pairs(with_seq=True, matches_only=False))
        q_pos = [x[0] for x in pileup.alignment.get_aligned_pairs(with_seq=True, matches_only=False)]
        r_pos = [x[1] for x in pileup.alignment.get_aligned_pairs(with_seq=True, matches_only=False)]
        
        if debug_reads: print(q_seq, '<-- query sequence (with gaps)')
        if debug_reads: print(r_seq, '<-- reference sequence (with gaps)')
            
        for pos in range(len(q_seq)):
            if r_pos[pos] is not None \
               and q_pos[pos] is not None \
               and (reference_name, r_pos[pos]+1) not in sites:
                s = (r_seq[pos].upper(), q_seq[pos].upper())
                if s not in subs: subs[s] = 0
                subs[s] += 1
                pass
            
            if debug_reads: print(q_seq)
            if debug_reads: print(r_seq)
            if debug_reads: print(' ' * pos + q_seq[pos], '<-- SNP' if r_pos[pos] is not None and (reference_name, r_pos[pos]+1) in sites else '')
            if debug_reads: print(' ' * pos + r_seq[pos])
            pass
        
        ## print the SNP site(s) if any exist in this read
        s_pos = [p for p in r_pos if p is not None and (reference_name, p+1) in sites]
        if debug_reads: print(''.join(q_seq[i]+'<- SNP' if p is not None and (reference_name, p+1) in sites else ' ' for i,p in enumerate(r_pos)))
        if debug_reads: print(''.join(' ' if x[1] is not None and ((reference_name, x[1]+1,) not in sites) else str(x[2]) + ' HEY' for x in pileup.alignment.get_aligned_pairs(with_seq=True, matches_only=False)))
        if debug_reads: print('SNP pos:', s_pos)
        
        if debug_reads: print('substitution matrix:', subs)
        pass
    return



# call_fun, out_fun, 
#def call_bases(bam, mincov, minbq, minmq, minlen, chrom, sites, flush_interval, limit, random_bases, third, subs):
# @profile
def call_bases(bam, sites,
               # third,
               results):
    """Sample bases in a given region of the genome based on the pileup
    of reads. If no coordinates were specified, sample from the whole BAM file.
    """
    calls = []
    i = 0
    # last_flush = 0
    times_called = 0

    subs = results['subs_matrix']
    results['read_names'] = dict()
    
    for i, col in enumerate(bam.pileup(contig=args.chrom, compute_baq=False,
                                       min_base_quality=args.minbq,
                                       min_mapping_quality=args.minmq)):

        if args.limit is not None and i > args.limit: break

        ## keep track of if this pileup column is actually in the sites we want to use
        col_in_sites = sites is not None and (col.reference_name, col.reference_pos + 1) in sites

        ######
        ## Do basic filtering (len, etc)
        ##### This is slow, b/c we're checking the same read over and over
        ##### We are also filtering "reads" based on bases at a given position, which doesn't really make sense..

        og_pileups = col.pileups
        if col_in_sites: og_bases = col.get_query_sequences(add_indels=True)
        # if col_in_sites: print(og_bases)
        # og_names = [p.alignment.query_name for p in og_pileups]

        pileups = []
        if col_in_sites: bases = []
        # read_names = []

        # if len(bases) >= args.mincov .....
        
        for i,pileup in enumerate(og_pileups):

            if  pileup.alignment.query_name in results['read_names'] \
                and results['read_names'][pileup.alignment.query_name] == 'fail':
                continue
            
            if len(pileup.alignment.query_sequence) < args.minlen:
                # if debug_reads: print('Filtering read due to minimum length', pileup.alignment.query_name, len(pileup.alignment.query_sequence))
                results['read_names'][pileup.alignment.query_name] = 'fail'
                continue

            ######
            ## compute substitutions matrix for all non-filtered reads, but only once per read
            
            if pileup.alignment.query_name not in results['read_names']:
                ### compute substitutions matrix for all non-filtered reads
                compute_subs_counts([pileup], results, args, col.reference_name)
                ## save this read name so we don't process it again
                results['read_names'][pileup.alignment.query_name] = 'pass'
                pass

            # filter out sites with no reads and sites with indels

            if col_in_sites and len(og_bases[i]) != 1:
                continue
            if col_in_sites and og_bases[i] == "*":
                continue
            if col_in_sites and og_bases[i] not in "ACGTacgt":
                continue

            pileups += [pileup]
            if col_in_sites: bases += [og_bases[i].upper()]
            # read_names += [read_names[i]]
            pass

        if len(pileups) == 0:
            continue

        #print(bases)
        #print(pileup)
        

        ######
        ## Now only process the sites that we requested
        if not col_in_sites:
            continue

        ######
        ### sometimes we want to just do random base sampling
        if not args.all_reads and len(bases) > 1:
            # random.seed(1984)
            # print(bases, len(bases))
            i = random.sample(range(len(bases)), 1)[0]
            # print(i)
            bases = [bases[i]]
            pileups = [pileups[i]]
            # read_names = [read_names[i]]
            pass
        
        ######
        ## Actually calculate site statistics
        process_site(col.reference_name,
                     col.reference_pos + 1,
                     sites,
                     get_ref_base(col),
                     bases,
                     # third,
                     results,
                     [pileup.alignment.query_name for pileup in pileups])
                
        # ## for outputing a file with every site listed
        # calls.append((
        #     col.reference_name,
        #     col.reference_pos + 1,
        #     get_ref_base(col),
        #     len(bases),
        #     call_fun(bases)
        # ))

        pass
    
    return


# def write_pileup(pileups, times_called, output, new_file, records_dict, report = True):

#     if not report: return
    
#     if new_file and times_called == 0: print("SITES\tchrom\tpos\tref\tpileup\tA\tC\tG\tT\tmatch", file=output)
#     for i in pileups.itertuples():
#         counts = Counter(i.call)
#         counts_str = '\t'.join(str(counts[i]) for i in 'ACGT')
#         if i.ref == i.call:
#             match = 'ref_match'
#         else:
#             match = 'ref_mismatch'
#             pass
#         records_dict[match] += 1
#         print(f"SITES\t{i.chrom}\t{i.pos}\t{i.ref}\t{''.join(i.call)}\t{counts_str}\t{match}", file=output)


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
    site_category_counts = defaultdict(int)
    
    with open(sites_file, "rt") as sf:
        for i, line in enumerate(sf):
            line = line.rstrip().split('\t')
            chrom,pos = line[0:2]
            # if i == 0 and not pos.isnumeric():
            if i == 0:
                sites['category_header'] = tuple(line[6:])

                print(' SNP file columns:  ', ', '.join(line), file=sys.stderr)
                print(' Categories from SNP file:  ', ', '.join(sites['category_header']), file=sys.stderr)
                if args.categories is not None:
                    category_indices = [i for i,cat in enumerate(sites['category_header']) if cat in args.categories]
                    sites['category_header'] = [sites['category_header'][i] for i in category_indices]
                    print(' Using only these categories (from --categories):  ', ', '.join(sites['category_header']), file=sys.stderr)
                    pass

                # print('SNP file header:', file=sys.stderr)
                # print('\t'.join(line), file=sys.stderr)
                # print('Categories from SNP file:', file=sys.stderr, sep='\t')
                # print('\t'.join(sites['category_header']), file=sys.stderr, sep='\t')

                continue
            
            a1_der, a2_anc = line[2:4]
            ref_fasta, a1_der, a2_anc, a3 = line[2:6]

            ####
            ## "category" is stuff like burden, n_3, etc
            ##
            
            category = tuple(line[6:])

            ## we're not always using all of the categories (--categories)
            ## in that case, filter down to the ones we need

            if args.categories is not None:
                category = tuple(category[i] for i in category_indices)
                pass

            # print('category', category)
            
            if add_chr:
                chrom = 'chr' + chrom
                pass
            pos = int(pos)
            sites[(chrom,pos)] = (a1_der, a2_anc, a3, ref_fasta, category)
            site_category_counts[category] += 1
            pass
        pass
    return (sites, site_category_counts)
    
def read_sites_og(sites_file, add_chr):
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



def report_stats(args, sites, site_category_counts, results):

    site_results = results['sites']

    all_cats = sorted(site_results.keys())
    print('all_cats', all_cats)

    all_stats_keys = site_results[all_cats[0]].keys()
    print('all_stats_keys', all_stats_keys)

    sum_stats = {k: sum(d[k] for d in site_results.values() if k in d) for k in all_stats_keys}

    print('sum_stats', sum_stats)

    if args.use_categories:
        # if 'category_header' in sites:
        print("\n")
        print('CATSITES',
              '\t'.join(sites['category_header']),
              'nsnps_in_cat',
              'cov_in_cat',
              'prop_snps_in_cat',
              'faunal_prop_b1',
              'faunal_prop_b1_w',
              # 'faunal_reads_b1',
              '\t'.join(all_stats_keys),
              sep='\t')

        if 'coverage' in args.strategy:
            # all_cats = [c for c,_ in sorted( ((c,int(c)) for c in all_cats),
            #                                  key = lambda _,x : x )]
            all_cats = [(str(x),) for x in range(-1,31)] + [('30+',), ('no-msa',), ('ensembl-fail',)]
            all_cats = [(str(x),) for x in range(-1,31)] + [('30+',), ('only-primate-msa',), ('ensembl-filter',)]
            baseline_cats = [(x,) for x in args.cov_baseline]
            pass

        ## TODO - THIS FAILS IF THERE'S LITERALLY NO DNA, I THINK?
        ## average coverage in "baseline" categories
        print([(site_category_counts[category], category) for category in baseline_cats])
        baseline_cov1 = sum([site_results[category]['n_reads'] / site_category_counts[category] for category in baseline_cats if site_category_counts[category] > 0]) / \
            sum(site_category_counts[category] > 0 for category in baseline_cats)
        ## weighted average coverage in "baseline" categories
        baseline_cov1_w = sum([site_results[category]['n_reads'] for category in baseline_cats if site_category_counts[category] > 0]) / \
            sum([site_category_counts[category] for category in baseline_cats if site_category_counts[category] > 0])
        print('baseline-cov', baseline_cov1)
        print('baseline-cov-w', baseline_cov1_w)

        n_reads_faunal = 0
        n_reads_faunal_w = 0
        n_reads_total = 0
        
        for category in all_cats:
            #print(category)
            #print(site_results[category])
            #print(all_stats_keys)

            if site_category_counts[category] > 0 and site_results[category]['n_reads'] > 0:
                cat_cov = site_results[category]['n_reads'] / site_category_counts[category]
                prop_faunal_cov   = max(0, (cat_cov - baseline_cov1)) / cat_cov
                prop_faunal_cov_w = max(0, (cat_cov - baseline_cov1_w)) / cat_cov
                prop_faunal_cov   = (cat_cov - baseline_cov1) / cat_cov
                prop_faunal_cov_w = (cat_cov - baseline_cov1_w) / cat_cov
                n_reads_faunal += site_results[category]['n_reads'] * prop_faunal_cov
                n_reads_faunal_w += site_results[category]['n_reads'] * prop_faunal_cov_w
                n_reads_total += site_results[category]['n_reads']
                n_reads_faunal_in_cat = site_results[category]['n_reads'] * prop_faunal_cov
            else:
                prop_faunal_cov = '-'
                prop_faunal_cov_w = '-'
                pass
                
            print('CATSITES',
                  '\t'.join(category),
                  site_category_counts[category],
                  site_results[category]['n_reads'] / site_category_counts[category] if site_category_counts[category] > 0 else '-',
                  site_results[category]['n_snps'] / site_category_counts[category] if site_category_counts[category] > 0 else '-',
                  prop_faunal_cov,
                  prop_faunal_cov_w,
                  '\t'.join(str(site_results[category][k]) for k in all_stats_keys),
                  # '\t'.join(str(site_results[category][k]) for k in all_stats_keys), THIS DIDN'T HAVE STR() BEFORE, BUT IT WORKED?? WHY???
                  sep='\t')
            pass

        print('cov-estimates', int(n_reads_faunal), n_reads_total, n_reads_faunal / n_reads_total)
        print('cov-estimates-w', int(n_reads_faunal_w), n_reads_total, n_reads_faunal_w / n_reads_total)

        pass

    print('\n')
    print('ALLSITES', '\t'.join(all_stats_keys), sep='\t')
    print('ALLSITES',
          '\t'.join(str(sum_stats[k]) for k in all_stats_keys),
          sep='\t')

    ######
    ## very basic stats (anc, der, a3, a4)

    ## config
    config_results = results['config']
    all_cats = sorted(config_results.keys())
    all_stats_keys = config_results[all_cats[0]].keys()
    sum_stats = {k: sum(d[k] for d in config_results.values() if k in d) for k in all_stats_keys}

    print('\n')
    print('ALLELES', '\t'.join(all_stats_keys), sep='\t')
    print('ALLELES', '\t'.join(str(sum_stats[k]) for k in all_stats_keys), sep='\t')

    if 'full_report_str' not in results:
        results['full_report_str'] = ''
        pass
    results['full_report_str'] += '\t' + 'alleles:' + '\t' + \
        str(sum_stats['n_anc']) + '\t' + str(sum_stats['n_der']) + \
        '\t' + str(sum_stats['n_der'] / (sum_stats['n_der'] + sum_stats['n_anc']))

    
    ######
    ## simulation "truth"
    if args.report_sim_truth is not None:

        sim_results = results['sim_source']
        all_spc = sorted(sim_results.keys())
        all_stats_keys = sim_results[all_spc[0]].keys()
        sum_stats = {k: sum(d[k] for d in sim_results.values() if k in d) for k in all_stats_keys}

        print('\n')
        print('SIMS_SPC', 'species', '\t'.join(all_stats_keys), sep='\t')
        for spc in all_spc:
            print('SIMS_SPC', spc,
                  '\t'.join(str(sim_results[spc][k]) for k in all_stats_keys),
                  sep='\t')
            pass

        print('\n')
        print('SIMS_SUM', '\t'.join(args.report_sim_truth), sep='\t')
        print('SIMS_SUM', '\t'.join(str(sim_results[spc]['N'] if spc in sim_results else 0) for spc in args.report_sim_truth), sep='\t')

        spc_sum = sum((sim_results[spc]['N'] if spc in sim_results else 0) for spc in args.report_sim_truth)
        print('SIMS_PROP', '\t'.join(args.report_sim_truth), sep='\t')
        print('SIMS_PROP', '\t'.join(str(sim_results[spc]['N']/spc_sum if spc in sim_results else 0) for spc in args.report_sim_truth), sep='\t')

        print('\n')
        print(sim_results)

        if 'human' in sim_results:
            h = sim_results['human']['N']
            h_a1a2 = sim_results['human']['anc_base'] + sim_results['human']['der_base']
        else:
            h = 0
            h_a1a2 = 0
            pass
        results['full_report_str'] += '\t' + 'truth:' + \
            '\t' + str(h) + \
            '\t' + str(h / sum(sim_results[spc]['N'] for spc in sim_results))
        results['full_report_str'] += '\t' + 'truth_a1a2:' + \
            '\t' + str(h_a1a2) + \
            '\t' + str(h_a1a2 / sum(sim_results[spc]['anc_base'] + sim_results[spc]['der_base'] for spc in sim_results)) + \
            '\t' + 'meth:' + str(args.lik_method)
        


    pass



def calc_lik_for_params(all_cats, config_results, f_anc, f_der, f_a4, subs, args):

    lik = 0

    #####
    ## A configuration is e.g. A,C,G,T for anc, der, a3, a4.
    ## So there are 16 possible configurations

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

            if args.lik_method == 0:
                ## anc - is this correct, to mult f_a4 with the subs rate? or at least, with this sub rate?
                temp_lik += (f_anc + f_a4) * subs[(anc,obs)]
                ## der
                temp_lik += (f_der + f_a4) * subs[(der,obs)]
                ## a3
                temp_lik += (1 - f_anc - f_der - 3*f_a4) * subs[(a3,obs)]
                ## a4
                temp_lik += (f_a4) * subs[(a4,obs)]

            elif args.lik_method == 1:
                ## THIRD.l1 does not work as well - this was an attempt at a different MLE function.
                
                ## anc - is this correct, to mult f_a4 with the subs rate? or at least, with this sub rate?
                temp_lik += (f_anc) * subs[(anc,obs)] + (f_a4 if anc == obs else 0)
                ## der
                temp_lik += (f_der) * subs[(der,obs)] + (f_a4 if der == obs else 0)
                ## a3
                temp_lik += (1 - f_anc - f_der - 3*f_a4) * subs[(a3,obs)]
                ## a4
                ## shouldn't have f_a4*subs here either?
                temp_lik += (f_a4) * subs[(a4,obs)]  + (f_a4 if a4 == obs else 0)


            elif args.lik_method == 2:
                ## THIRD.l2 is currently the best method!

                ## HERE f_a4 is the "multiple mutation" rate, from the ancestral allele.
                ## We assume that there isn't enough time for the derived allele to get multiple mutations
                
                ## anc
                temp_lik += (f_anc) * subs[(anc,obs)]
                ## der
                temp_lik += (f_der + f_a4) * subs[(der,obs)]
                ## a3
                temp_lik += (1 - f_anc - f_der - 2*f_a4) * subs[(a3,obs)]
                ## a4
                temp_lik += (f_a4) * subs[(a4,obs)]

            elif args.lik_method == 3:
                ## THIRD.l3 does not work as well - this was an attempt at a different MLE function.

                ## HERE f_a4 is the "multiple mutation" rate, from the ancestral allele.
                ## We assume that there isn't enough time for the derived allele to get multiple mutations

                ## given an observed base "obs", this base could have originated from any of the 4 bases.
                ## If it originates from a different base than 
                
                ## anc
                temp_lik += (f_anc) * subs[(anc,obs)]
                ## der
                temp_lik += (f_der + f_a4) * subs[(der,obs)]
                ## a3
                temp_lik += (1 - f_anc - f_der - 1*f_a4) * subs[(a3,obs)]
                ## a4
                temp_lik += (f_a4) * subs[(a4,obs)]
                pass

            
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

    if args.breaks_anc_der is None:
        breaks_anc_der = np.arange(0, 1.1, args.step)
    else:
        breaks_anc_der = args.breaks_anc_der
        pass
    
    if args.breaks_a4 is None:
        breaks_a4 = np.arange(0, .1, args.step)
    else:
        breaks_a4 = args.breaks_a4
        pass
        
    for f_anc in breaks_anc_der:
        for f_der in breaks_anc_der:
            for f_a4 in breaks_a4:

                if f_anc + f_der + 3*f_a4 > 1:
                    continue
                if f_anc + f_der == 0:
                    continue
                
                lik = calc_lik_for_params(all_cats, config_results, f_anc, f_der, f_a4, subs, args)
                # print('LIKELIHOOD', f_anc, f_der, f_a4, lik, f_anc/(f_anc+f_der), f_der/(f_anc+f_der), sep='\t')
                all_params += [(f_anc, f_der, f_a4)]
                all_liks += [lik]
                
                pass
            pass
        pass

    print('LIKELIHOOD', 'f_anc', 'f_der', 'f_a4', 'lik', 'rel_f_anc', 'rel_f_der', sep='\t')

    for lik, params in sorted(zip(all_liks, all_params)):
        f_anc, f_der, f_a4 = params
        print('LIKELIHOOD', f_anc, f_der, f_a4, lik, f_anc/(f_anc+f_der), f_der/(f_anc+f_der), sep='\t')
        pass

    max_lik_i = np.argmax(all_liks)
    f_anc, f_der, f_a4 = all_params[max_lik_i]
    max_lik = all_liks[max_lik_i]
    print('MAXIMUM_LIK', f_anc, f_der, f_a4, max_lik, f_anc/(f_anc+f_der), f_der/(f_anc+f_der), sep='\t', end='')
    if 'full_report_str' in results:
        print(results['full_report_str'])
    else:
        print()
        pass

    return(max_lik, max_lik_i, all_liks, all_params)


def compute_subs_probs(subs_matrix):
    subs_probs = dict()
    print(' Substitution matrix + probabilities')
    if args.debug: print('subs_matrix', subs_matrix)
    for b1 in list('ACGT'):
        s = sum(subs_matrix[(b1,b2)] for b2 in list('ACGT'))
        for b2 in list('ACGT'):
            subs_probs[(b1,b2)] = subs_matrix[(b1,b2)] / s 
            print('', b1, b2, subs_matrix[(b1,b2)], subs_probs[(b1,b2)], sep='\t')
            pass
        pass
    return subs_probs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Call alleles from a BAM file using various criteria")
    parser.add_argument("--bam", help="BAM file to sample from", required=True)
    parser.add_argument("--chrom", help="Chromosome to sample from")
    # parser.add_argument("--proportion", help="Required proportion of the majority allele", type=check_range, default=0.5)
    parser.add_argument("--seed", help="Set seed for random allele sampling [random]", type=int)
    # parser.add_argument("--mincov", help="Minimum coverage", type=int, default=1)
    parser.add_argument("--minbq", help="Minimum base quality (default 13)", type=int, default=13)
    parser.add_argument("--minmq", help="Minimum read mapping quality (default 25)", type=int, default=25)
    parser.add_argument("--minlen", help="Minimum read length (default 35)", type=int, default=35)
    parser.add_argument("--add-chr", help="Add 'chr' to the beginining of every site's chromosome, to match bam file.", action='store_true')
    parser.add_argument("--sample-name", help="Add sample ID as an additional column to output tables [NEEDS TO BE IMPLEMENTED!]")
    parser.add_argument("--output", help="Output file name [NOT SURE IF WE USE THIS? COULD BE OUTPUT BASE?]")
    parser.add_argument("--sites", help="Restrict output to this list of sites (see documentation for file format)")
    # parser.add_argument("--flush", help="Print to file every N bases", type=int, default=100000)
    parser.add_argument("--limit", help="Stop after processing N bases from bam", type=int, default=None)
    parser.add_argument("--lik-method", help="Method for computing likelihood (devel)", type=int, default=0)
    # parser.add_argument("--third_file", help="Defines which nucleotide is the third allele")
    parser.add_argument("--all-reads", help="Use all reads (default is to do random read sampling at each SNP)", action='store_true')
    parser.add_argument("--step", help="Step for grid likelihood search", type=float, default=0.01)
    parser.add_argument("--breaks-anc-der", help="Explicit breaks for grid likelihood search (anc and der proportions)", nargs='+', type=float, default=None)
    parser.add_argument("--breaks-a4", help="Explicit breaks for grid likelihood search (a4 allele proportion)", nargs='+', type=float, default=None)
    parser.add_argument("--report-reads", help="Report base state (ancestral, derived, a3, a4) for every read, along with read name. Mostly for debugging.", action='store_true')
    parser.add_argument("--report-sim-truth",
                        help='Try to parse the simulated species from a read name, and report the "truth". ' +
                        'Assumes format "species_rest_of_name". Give a list of species to report.',
                        nargs = '+', default = None)
    # parser.add_argument("--by-cats", help="Print statistics for each SNP category", action=argparse.BooleanOptionalAction)
    parser.add_argument("-use-cats", "--use-categories", help="Print statistics for each SNP category (categories taken from columns 6+ in SNP file; e.g. b_3 [mammalian diversity])", action='store_true')
    parser.add_argument('-cats', "--categories", help="Use these categories (subset of columns from header)", nargs = '+', default = None)
    parser.add_argument('-cov-baseline', "--cov-baseline", help="Use these category states to calculate a coverage baseline (mostly debugging)", nargs = '+', default = [str(x) for x in range(14,30)] + ['30+', 'only-primate-msa'])
    parser.add_argument('-debug', "--debug", help="Print some (confusing) debugging information.", action='store_true')
    parser.add_argument("--strategy", help="How to calculate faunal mismaping proportions?",
                        choices=["likelihood", "derived-alleles", "coverage", "none"],
                        required=False,
                        nargs='+',
                        default=['derived-alleles'])


    args = parser.parse_args()


    if args.output is None:
        new_file = True
        output = sys.stdout
    else:
        new_file = not os.path.isfile(args.output)
        output = open(args.output, "w" if new_file else "a")
        pass

    bam = pysam.AlignmentFile(args.bam)

    if not bam.has_index():
        print("BAM file index is missing", file=sys.stderr)
        sys.exit(1)
        pass


    # if args.third_file:
    #     third = read_third_file(args.third_file)
    #     pass


    ###########
    ## DO WE EVER WANT TO RUN THIS PROGRAM WITHOUT A SET OF SITES TO CONSIDER??
    ##
    
    if args.sites is not None:
        print('\nReading SNP file..', file=sys.stderr)
        sites, site_category_counts = read_sites(args.sites, args.add_chr)
        # print('all sites', sites)
    else:
        sites = None
        site_category_counts = None
        pass

    # check_dicts(third, sites)

    # subs_matrix = dict()
    results = dict()
    results['subs_matrix'] = dict()
    results['config'] = dict()
    results['sim_source'] = dict()

    results['sites'] = dict()
    results['sites'] = defaultdict(lambda : dict({'n_reads_a12' : 0,
                                                  'n_reads' : 0,
                                                  'n_snps' : 0,
                                                  'n_der' : 0,
                                                  'n_anc' : 0,
                                                  'ref_bam' : 0,
                                                  # 'ref_nonmatch' : 0,
                                                  'n_first' : 0,
                                                  'n_second' : 0,
                                                  'n_third' : 0,
                                                  'n_fourth' : 0,
                                                  'ref_bam_not_a3' : 0}))

    # call_fun, out_fun, 
    # call_bases(bam, args.mincov,
    #            args.minbq, args.minmq, args.minlen,
    #            args.chrom, sites, args.flush, args.limit,
    #            args.random_bases, third, subs = subs_matrix)
    print('\nProcessing bam file..', file=sys.stderr)
    call_bases(bam = bam,
               sites = sites,
               # third = third,
               results = results)

    print('\nComputing substitution matrix..', file=sys.stderr)
    subs_probs = compute_subs_probs(results['subs_matrix'])

    print('\nReporting basic statistics..', file=sys.stderr)
    report_stats(args, sites, site_category_counts, results)

    if 'likelihood' in args.strategy:
        print('\nCalculating maximum likelihood mismapping estimate..', file=sys.stderr)
        calc_lik(args, sites, results['config'], subs_probs)
        pass

    pass
