#!/usr/bin/env python3
"""Issue #175 regression tests using exclusively generated, invented DNA.

python3 test/test_nm_md.py --bwa ./bwa-mem2 [--baseline /path/to/original]
Optionally pass --samtools to compare the raw output with stock calmd.
"""
import argparse
from pathlib import Path
import random
import re
import subprocess
import tempfile
import unittest


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def expected_tags(sequence, reference, position, cigar):
    """Calculate SAM NM/MD independently from original FASTA and full read."""
    query, ref, nm, matches, md = 0, position - 1, 0, 0, []
    for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        length = int(length)
        if op in 'M=X':
            for a, b in zip(sequence[query:query+length], reference[ref:ref+length]):
                b = b.upper()
                if a in 'ACGT' and b in 'ACGT' and a == b:
                    matches += 1
                else:
                    nm += 1
                    md.extend((str(matches), b))
                    matches = 0
            query += length
            ref += length
        elif op == 'D':
            nm += length
            md.extend((str(matches), '^'+reference[ref:ref+length].upper()))
            matches = 0
            ref += length
        elif op == 'I':
            nm += length
            query += length
        elif op in 'SH':
            # sequence is the full input read, including hard-clipped bases.
            query += length
        elif op == 'N':
            ref += length
    md.append(str(matches))
    assert query == len(sequence)
    assert ref <= len(reference)
    return nm, ''.join(md)


def parse_sam(path):
    return [line.split('\t') for line in path.read_text().splitlines()
            if line and not line.startswith('@')]


def tags(row):
    return {field.split(':', 2)[0]: field.split(':', 2)[2] for field in row[11:]}


class NmMdTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='bwa175-invented-')
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def command(self, args, output):
        with output.open('w') as handle:
            result = subprocess.run([str(a) for a in args], stdout=handle,
                                    stderr=subprocess.PIPE, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)

    def align(self, references, reads, mates=None, options=()):
        reference = self.root/'reference.fa'
        reference.write_text(''.join('>'+n+'\n'+s+'\n' for n, s in references.items()))
        inputs = {}
        read_paths = []
        for mate, sequences in enumerate((reads,) if mates is None else (reads, mates)):
            path = self.root/f'reads{mate}.fq'
            path.write_text(''.join('@'+name+'\n'+sequence+'\n+\n'+'I'*len(sequence)+'\n'
                                    for name, sequence in sequences.items()))
            read_paths.append(path)
            for name, sequence in sequences.items():
                inputs[name, mate+1 if mates is not None else 0] = sequence.upper()
        self.command([ARGS.bwa, 'index', reference], self.root/'index.log')

        def run(aligner, prefix):
            sam = self.root/(prefix+'.sam')
            self.command([aligner, 'mem', '-t', '1', *options, reference, *read_paths], sam)
            return sam, parse_sam(sam)

        sam, rows = run(ARGS.bwa, 'fixed')
        self.assertTrue(rows)
        mapped = [row for row in rows if not int(row[1]) & 4]
        self.assertEqual({row[0] for row in mapped}, set(reads))
        for row in mapped:
            flag = int(row[1])
            mate = (1 if flag & 64 else 2) if flag & 1 else 0
            sequence = inputs[row[0], mate]
            oriented = reverse_complement(sequence) if flag & 16 else sequence
            nm, md = expected_tags(oriented, references[row[2]], int(row[3]), row[5])
            with self.subTest(read=row[0], contig=row[2], cigar=row[5], flag=flag):
                self.assertEqual(int(tags(row)['NM']), nm)
                self.assertEqual(tags(row)['MD'], md)
            # NM copies embedded in alternative and supplementary tags must
            # also be based on the original reference, not the packed bases.
            for tag in ('SA', 'XA'):
                for hit in tags(row).get(tag, '').strip(';').split(';'):
                    if not hit:
                        continue
                    fields = hit.split(',')
                    if tag == 'SA':
                        contig, pos, strand, cigar, _, hit_nm = fields
                    else:
                        contig, signed_pos, cigar, hit_nm = fields
                        strand, pos = signed_pos[0], signed_pos[1:]
                    hit_read = reverse_complement(sequence) if strand == '-' else sequence
                    expected_nm, _ = expected_tags(hit_read, references[contig], int(pos), cigar)
                    self.assertEqual(int(hit_nm), expected_nm, (tag, hit))

        if ARGS.baseline:
            _, baseline = run(ARGS.baseline, 'baseline')
            self.assertEqual(len(rows), len(baseline))
            for old, new in zip(baseline, rows):
                self.assertEqual(old[:11], new[:11], 'alignment fields changed')
                old_tags, new_tags = tags(old), tags(new)
                for tag in ('NM', 'MD'):
                    old_tags.pop(tag, None)
                    new_tags.pop(tag, None)
                for tag in ('SA', 'XA'):
                    for value in (old_tags, new_tags):
                        if tag in value:
                            value[tag] = [''.join(hit.rsplit(',', 1)[:1])
                                          for hit in value[tag].strip(';').split(';')]
                self.assertEqual(old_tags, new_tags, 'non-NM/MD alignment tags changed')
            if all(set(s.upper()) <= set('ACGT') for s in references.values()):
                self.assertEqual(rows, baseline, 'no-ambiguity control changed')

        if ARGS.samtools:
            # Use SAM text input so no BAM library is needed by these tests.
            calmd = self.root/'calmd.sam'
            self.command([ARGS.samtools, 'calmd', sam, reference], calmd)
            for old, new in zip(rows, parse_sam(calmd)):
                if old[9] != '*' and not int(old[1]) & 4:
                    self.assertEqual(tags(old)['NM'], tags(new)['NM'])
                    self.assertEqual(tags(old)['MD'], tags(new)['MD'])
        return mapped

    def test_tiny_reference_N_on_both_strands(self):
        rng = random.Random(175)
        reference = list(''.join(rng.choices('ACGT', k=240)))
        reference[119] = 'N'
        reference = ''.join(reference)
        reads = {}
        for base in 'ACGTN':
            sequence = reference[69:119] + base + reference[120:170]
            reads['synthetic_'+base] = sequence
            reads['synthetic_'+base+'_reverse'] = reverse_complement(sequence)
        rows = self.align({'synthetic': reference}, reads)
        self.assertEqual(len(rows), 10)
        for row in rows:
            self.assertEqual((row[3], row[5], tags(row)['NM'], tags(row)['MD']),
                             ('70', '101M', '1', '50N50'))

    def test_no_ambiguity_control(self):
        rng = random.Random(176)
        reference = ''.join(rng.choices('ACGT', k=1000))
        read = reference[200:400]
        self.align({'control': reference}, {'exact': read, 'reverse': reverse_complement(read),
                   'read_N': read[:75]+'N'+read[76:]})

    def test_iupac_runs_indels_and_clipping(self):
        rng = random.Random(177)
        references, reads = {}, {}
        for index, code in enumerate('NnRYSWKMBDHVXrys'):
            truth = ''.join(rng.choices('ACGT', k=800))
            reference = truth[:299] + code + truth[300:]
            references['ambiguity_'+str(index)] = reference
            read = truth[200:450]
            reads['forward_'+str(index)] = read
            reads['reverse_'+str(index)] = reverse_complement(read)
        for kind in ('insertion', 'deletion', 'clipping', 'runs'):
            truth = ''.join(rng.choices('ACGT', k=1000))
            reference = truth[:295]+'NNN'+truth[298:350]+'nRyn'+truth[354:]
            references[kind] = reference
            read = truth[200:500]
            if kind == 'insertion':
                read = read[:150]+'AACCGG'+read[150:]
            elif kind == 'deletion':
                read = read[:94]+read[105:]
            elif kind == 'clipping':
                read = 'A'*45+read+'T'*45
            reads[kind] = read
            reads[kind+'_reverse'] = reverse_complement(read)
        rows = self.align(references, reads)
        for name, op in (('insertion', 'I'), ('deletion', 'D'), ('clipping', 'S')):
            self.assertTrue(any(row[0] == name and op in row[5] for row in rows))
        self.assertTrue(any('^' in tags(row)['MD'] and 'N' in tags(row)['MD'].split('^')[1]
                            for row in rows if row[0] == 'deletion'))

    def test_paired_reads(self):
        rng = random.Random(178)
        truth = ''.join(rng.choices('ACGT', k=16000))
        reference = list(truth)
        reads, mates = {}, {}
        for index in range(24):
            start = 200 + index*600
            for offset in (70, 71, 310, 311):
                reference[start+offset] = 'N'
            reads['pair_'+str(index)] = truth[start:start+150]
            mates['pair_'+str(index)] = reverse_complement(truth[start+250:start+400])
        rows = self.align({'paired': ''.join(reference)}, reads, mates)
        self.assertEqual(len(rows), 48)
        self.assertTrue(all(int(row[1]) & 1 for row in rows))

    def test_secondary_supplementary_and_alternative_tags(self):
        rng = random.Random(179)
        duplicate = ''.join(rng.choices('ACGT', k=3000))
        unique = ''.join(rng.choices('ACGT', k=3000))
        refs = {'copyA': duplicate[:1070]+'N'+duplicate[1071:],
                'copyB': duplicate[:1070]+'n'+duplicate[1071:],
                'unique': unique[:2070]+'R'+unique[2071:]}
        reads = {'duplicate': duplicate[1000:1200],
                 'chimeric': duplicate[1000:1150]+unique[2000:2150]}
        rows = self.align(refs, reads)
        self.assertTrue(any('XA' in tags(row) for row in rows))
        self.assertTrue(any('SA' in tags(row) and int(row[1]) & 2048 for row in rows))
        rows = self.align(refs, reads, options=('-a',))
        self.assertTrue(any(int(row[1]) & 256 for row in rows))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bwa', type=Path, required=True)
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--samtools', type=Path)
    ARGS, remaining = parser.parse_known_args()
    for name in ('bwa', 'baseline', 'samtools'):
        if getattr(ARGS, name):
            setattr(ARGS, name, getattr(ARGS, name).resolve())
    unittest.main(argv=[__file__, *remaining])
