#!/usr/bin/env python3
"""Call six-SNP MAPT haplotypes from a phased multi-sample VCF.

Uses VCF record IDs, not fixed positions. Input must include six rsID-labelled
biallelic or multiallelic SNP records and phased diploid GTs. Output has two
rows per sample (sample_id_A and sample_id_B), one per phased chromosome.

Classification references (allele order: rs1467967, rs242557, rs3785883,
rs2471738, rs8070723, rs7521):
  Allen et al., Alzheimer's Research & Therapy 2014, Table 2.
  Labbé et al., Parkinsonism & Related Disorders 2016, Table 3.
  Babić Leko et al., Brain and Behavior 2018, Table 2.

An empty study cell means that the exact six-base string is absent from that
study's cited table. The H1_or_H2 column is a separate rs8070723 allele call.
This script neither infers phase nor assigns copy-number/structural haplotypes.
"""

import argparse
import csv
import gzip
import subprocess
import sys

SNPS = ('rs1467967', 'rs242557', 'rs3785883', 'rs2471738', 'rs8070723', 'rs7521')

ALLEN = {
    'AGGCGG': 'A (H2a)', 'GGGCAA': 'B (H1b)', 'AAGTAG': 'C (H1c)',
    'AAGCAA': 'D (H1d)', 'AGGCAA': 'E (H1e)', 'GAACAA': 'G',
    'AGACAA': 'H', 'GAGCAA': 'I', 'AGGCAG': 'J', 'AGACAG': 'L',
    'GAGCAG': 'M', 'GGACAG': 'N', 'AAACAA': 'O', 'GGGTAG': 'P',
    'AGGTAG': 'R', 'AAGCAG': 'U', 'GGATAG': 'V', 'GGGCAG': 'W',
    'GAATAG': 'X', 'AAATAG': 'Y*', 'GAGTAG': 'Z*',
}

LABBE = {
    'AAGTAG': 'H1C', 'AGGCGG': 'H2', 'GGGCAA': 'H1B',
    'AGGCAA': 'H1E', 'AAGCAA': 'H1D', 'GAGCAA': 'H1I',
    'AGACAA': 'H1H', 'AGACAG': 'H1L', 'GAGCAG': 'H1M',
    'AAGCAG': 'H1U', 'AAACAA': 'H1O', 'AAATAG': 'H1y',
    'GGGTAG': 'H1P', 'GAATAG': 'H1x', 'GGACAA': 'H1F',
    'GGATAG': 'H1v', 'AGGTAG': 'H1R', 'GAACAA': 'H1G',
    'AAGTAA': 'H1Q', 'AGGCAG': 'H1J', 'GGGCAG': 'H1S',
    'GAGTAG': 'H1z', 'GGACAG': 'H1N', 'AAACAG': 'H1K',
}

LEKO = {
    'GGGCAA': 'H1B', 'AAGTAG': 'H1c', 'AGGCAA': 'H1E',
    'AAGCAA': 'H1D', 'AGACAG': 'H1l', 'AAGCAG': 'H1u',
    'AGACAA': 'H1h', 'GAGCAA': 'H1i', 'AGGCAG': 'H1J',
    'GAGCAG': 'H1m', 'AGATAG': 'H1T', 'AAATAG': 'H1k',
    'GGGCAG': 'H1jj', 'AAACAA': 'H1o', 'GGATAG': 'H1v',
    'AAGTAA': 'H1q', 'GAACAA': 'H1g', 'GAATAG': 'H1x',
    'GAACAG': 'H1y', 'GGACAA': 'H1F', 'GGGTAG': 'H1p',
    'GAGTAG': 'H1aa', 'AGGTAG': 'H1r', 'AGGCGG': 'H2A',
    'AGACGG': 'H2gg', 'AAGCGG': 'H2ff', 'AGGCGA': 'H2kk',
    'GGGCGA': 'H2w',
}


def vcf_lines(path, region):
    if region or path.endswith('.bcf'):
        command = ['bcftools', 'view', '-Ov']
        if region:
            command += ['-r', region]
        command.append(path)
        try:
            process = subprocess.Popen(command, stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE, text=True)
        except FileNotFoundError as exc:
            raise ValueError('bcftools is needed for --region or BCF input') from exc
        try:
            yield from process.stdout
        finally:
            process.stdout.close()
            error = process.stderr.read()
            returncode = process.wait()
            process.stderr.close()
            if returncode:
                raise ValueError(f'bcftools failed: {error.strip()}')
    else:
        opener = gzip.open if path.endswith('.gz') else open
        with opener(path, 'rt') as handle:
            yield from handle


def load_six_variants(path, region):
    samples = None
    records = {}
    for line in vcf_lines(path, region):
        if line.startswith('#CHROM'):
            samples = line.rstrip('\r\n').split('\t')[9:]
            if not samples or len(samples) != len(set(samples)):
                raise ValueError('VCF must contain unique sample IDs')
            continue
        if line.startswith('#'):
            continue
        if samples is None:
            raise ValueError('No VCF #CHROM header found')
        fields = line.rstrip('\r\n').split('\t')
        matches = set(fields[2].split(';')).intersection(SNPS)
        if not matches:
            continue
        if len(matches) != 1:
            raise ValueError(f'A record matches multiple target rsIDs: {fields[2]}')
        snp = matches.pop()
        if snp in records:
            raise ValueError(f'Multiple records have the ID {snp}')
        if len(fields) != len(samples) + 9:
            raise ValueError(f'Unexpected number of sample fields at {snp}')
        alts = [] if fields[4] == '.' else [a.upper() for a in fields[4].split(',')]
        alleles = [fields[3].upper(), *alts]
        if any(len(a) != 1 or a not in 'ACGT' for a in alleles):
            raise ValueError(f'{snp} is not a single-base SNP: {alleles}')
        formats = fields[8].split(':')
        if 'GT' not in formats:
            raise ValueError(f'{snp} has no GT field')
        gt_idx = formats.index('GT')
        ps_idx = formats.index('PS') if 'PS' in formats else None
        calls = []
        for sample, value in zip(samples, fields[9:]):
            tokens = value.split(':')
            gt = tokens[gt_idx] if gt_idx < len(tokens) else ''
            if '|' not in gt:
                raise ValueError(f'{sample} at {snp}: unphased or missing GT {gt!r}')
            indexes = gt.split('|')
            if len(indexes) != 2 or any(not i.isdigit() or int(i) >= len(alleles)
                                        for i in indexes):
                raise ValueError(f'{sample} at {snp}: invalid diploid GT {gt!r}')
            pair = tuple(alleles[int(i)] for i in indexes)
            ps = tokens[ps_idx] if ps_idx is not None and ps_idx < len(tokens) else None
            calls.append((pair, ps))
        records[snp] = calls
    if samples is None:
        raise ValueError('No VCF #CHROM header found')
    missing = [s for s in SNPS if s not in records]
    if missing:
        raise ValueError('Missing rsID records: ' + ', '.join(missing) +
                         '. Check VCF IDs and, if used, --region and genome build.')
    return samples, records


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='multi-sample phased VCF, VCF.gz, or BCF')
    parser.add_argument('-o', '--output', required=True, help='output TSV path')
    parser.add_argument('--region', help='optional indexed region, e.g. chr17:45900000-46050000'
                        ' for GRCh38; requires bcftools and an indexed VCF/BCF')
    args = parser.parse_args()
    samples, records = load_six_variants(args.vcf, args.region)
    # Validate all samples before opening the output, to avoid partial results.
    for i, sample in enumerate(samples):
        phase_sets = {records[s][i][1] for s in SNPS
                      if records[s][i][0][0] != records[s][i][0][1]
                      and records[s][i][1] not in (None, '.', '')}
        if len(phase_sets) > 1:
            raise ValueError(f'{sample}: heterozygous SNPs are in different PS blocks: '
                             + ', '.join(sorted(phase_sets)))
    with open(args.output, 'w', newline='') as handle:
        out = csv.writer(handle, delimiter='\t', lineterminator='\n')
        out.writerow(('sample_id', *SNPS, 'hap_as_bases', 'allen', 'labbe',
                      'leko', 'H1_or_H2'))
        for i, sample in enumerate(samples):
            for chromosome, suffix in enumerate(('A', 'B')):
                bases = [records[s][i][0][chromosome] for s in SNPS]
                seq = ''.join(bases)
                clade = {'A': 'H1', 'G': 'H2'}[bases[4]]
                out.writerow((f'{sample}_{suffix}', *bases, seq, ALLEN.get(seq, ''),
                              LABBE.get(seq, ''), LEKO.get(seq, ''), clade))
    print(f'Wrote {len(samples) * 2} haplotype rows for {len(samples)} samples '
          f'to {args.output}', file=sys.stderr)


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError) as exc:
        sys.exit(f'Error: {exc}')
