import math
import pysam
import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chisquare
from itertools import groupby
from collections import Counter

warnings.simplefilter("error", RuntimeWarning)

"""Open VCF"""

vcf = pysam.VariantFile('data/experiment_11.vcf')

"""Count samples"""

samples = vcf.header.samples
cases = {s for s in vcf.header.samples if s.startswith('case')}
controls = {s for s in vcf.header.samples if s.startswith('control')}

len_samples, len_cases, len_controls = len(samples), len(cases), len(controls)

print('Total={} Case={} Control={}'.format(len_samples, len_cases, len_controls))

"""Fetch all mutations"""

variants = list(vcf.fetch())

"""Count mutations and show total count per chromosome table"""

muts_per_chrom = Counter([v.chrom for v in variants])
len_variants = sum(muts_per_chrom.values())

freq_data = [[k, v] for k, v in muts_per_chrom.items()]
freq_data.append(['Total mutations', len_variants])

freq_labels = ['Chromosome', 'Mutations']
freq_table = plt.table(cellText=freq_data, colLabels=freq_labels, loc='center')

plt.axis('off')
plt.show()

"""Extract Minor Allele Frequency (MFA) and sampled genotypes per mutation"""

def get_genotypes():
    maf_list, ref_list, alt_list, het_list = [], [], [], []

    for v in variants:
        counts = Counter([sum(v.samples[s]['GT']) for s in samples])
        refs, alts, hets = counts[0], counts[2], counts[1]

        maf_list.append((2 * min(refs, alts) + hets) / (2 * len_samples))
        ref_list.append(refs / len_samples)
        alt_list.append(alts / len_samples)
        het_list.append(hets / len_samples)

    return maf_list, ref_list, alt_list, het_list

maf, ref, alt, het = get_genotypes()

"""Show histograms"""

def show_hist(data, label, index):
    plt.subplot(1, 4, index)
    plt.hist(data, bins=25, label=label)
    plt.legend()

plt.figure(figsize=(16, 4))
show_hist(maf, 'Minor Allele Frequency', 1)
show_hist(ref, 'Homozygous Reference', 2)
show_hist(alt, 'Homozygous Alternative', 3)
show_hist(het, 'Heterozygous', 4)
plt.show()

"""Run Chi-Square test for each mutation"""

def run_chi2(v):
    case_counts = Counter([sum(v.samples[s]['GT']) for s in cases])
    case_refs, case_alts = case_counts[0], case_counts[2]

    control_counts = Counter([sum(v.samples[s]['GT']) for s in controls])
    control_refs, control_alts = control_counts[0], control_counts[2]

    total = case_refs + case_alts + control_refs + control_alts

    try:
        observed = [case_refs, case_alts, control_refs, control_alts]
        expected = [(case_refs + case_alts) * (case_refs + control_refs) / total,
                    (case_refs + case_alts) * (case_alts + control_alts) / total,
                    (control_refs + control_alts) * (case_refs + control_refs) / total,
                    (control_refs + control_alts) * (case_alts + control_alts) / total]
        return chisquare(observed, expected, ddof=2).pvalue
    except Warning:
        return None

chi2_raw = [(v.chrom, v.pos, run_chi2(v)) for v in variants]
chi2 = [x for x in chi2_raw if x[2] is not None]
results = {k: list(v) for k, v in groupby(chi2, key=lambda x: x[0])}

"""Calculate chromosome lengths"""

def chrom_lengths():
    chrom_groups = groupby(variants, key=lambda v: v.chrom)
    return [max([v.pos for v in g]) for _, g in chrom_groups]

"""Show Manhattan plot"""

def show_manhattan(bonferroni = False):
    lengths = chrom_lengths()
    widths = [math.ceil(100 * x / sum(lengths)) for x in lengths]

    plt.figure(figsize=(16, 4))

    for i, chrom in enumerate(results):
        pos = [x[1] for x in results[chrom]]
        p = [-np.log10(x[2]) for x in results[chrom]]

        plt.subplot2grid((1, sum(widths)), (0, sum(widths[:i])), colspan=widths[i])
        plt.scatter(pos, p)
        plt.ylim(0, 13)
        plt.axhline(-np.log10(0.01 / len_variants if bonferroni else 0.01), c='r')

    plt.show()

show_manhattan()
show_manhattan(True)

"""Extract regions with most significant mutations"""

def get_regions():
    regions = []
    
    for i, chrom in enumerate(results):
        pos = [x[1] for x in results[chrom]]
        p = [-np.log10(x[2]) for x in results[chrom]]

        signif = [x for x in p if x > -math.log10(0.01 / len_variants)]

        pos_min, pos_max = float('inf'), float('-inf')
    
        for j, x in enumerate(p):
            if x in signif:
                pos_min, pos_max = min(pos_min, pos[j]), max(pos_max, pos[j])

        if pos_min != float('inf') and pos_max != float('-inf'):
            regions.append({'chr' : i + 1, 'min' : pos_min, 'max' : pos_max})

    return regions

"""Check if mutation is inside region"""

def isinside(v, r):
    return int(v.chrom) == r['chr'] and r['min'] <= v.pos < r['max']

signif = [v for r in get_regions() for v in variants if isinside(v, r)]

"""Check Hardy-Weinberg equilibrium for each significant mutation"""

def check_hwe(v):
    counts = Counter([sum(s['GT']) for s in v.samples.values() if None not in s['GT']])
    total = sum(counts.values())
    
    p = (2 * counts[0] + counts[1]) / (2 * total)
    q = (2 * counts[2] + counts[1]) / (2 * total)

    try:
        observed = [counts[0], counts[1], counts[2]]
        expected = [p * p * total, 2 * p * q * total, q * q * total]
        return chisquare(observed, expected, ddof=1).pvalue
    except Warning:
        return None

hwe_raw = [check_hwe(v) for v in signif]
hwe = [x > (0.01 / len_variants) for x in hwe_raw if x is not None]

"""Show HWE results for each significant mutation"""

hwe_data = [[v.chrom, v.pos, h] for (v, h) in zip(signif, hwe)]
hwe_labels = ['Mutation Chromosome', 'Mutation Position', 'HWE Applies']
hwe_table = plt.table(cellText=hwe_data, colLabels=hwe_labels, loc='center')

plt.axis('off')
plt.show()

"""Show all significant regions

Extract relevant genes in this regions using [UCSC Genome browser](https://genome.ucsc.edu/cgi-bin/hgGateway?redirect=manual&source=genome.ucsc.edu)
"""

for r in get_regions():
    print('chr{}:{}-{}'.format(r['chr'], r['min'], r['max']))
