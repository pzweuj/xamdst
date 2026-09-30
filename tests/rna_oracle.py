#!/usr/bin/env python3
"""Independent per-base RNA oracle and transactional integration fixtures."""
import gzip
import json
import os
from pathlib import Path
import re
import resource
import signal
import shutil
import subprocess
import sys
import tempfile
from collections import Counter

BINARY = str(Path(sys.argv[1]).resolve())
DNA_FILES = {'coverage.report', 'coverage.report.json', 'cumu.plot', 'insert.plot',
             'chromosome.report', 'region.tsv.gz', 'depth.tsv.gz', 'uncover.bed'}
RNA_FILES = DNA_FILES | {'splice.tsv.gz', 'distribution.tsv'}
LENGTHS = {'chr1': 120, 'chr2': 80, 'chr3': 40}
# This oracle uses small explicit sets of genomic positions, independently of
# the production interval parser, merging, lookup and CIGAR implementation.
EXONS = {'chr1': set(range(15)) | set(range(20, 30)) | set(range(60, 70)),
         'chr2': set(range(8)) | set(range(18, 25)), 'chr3': set(range(2, 5))}
BODIES = {'chr1': set(range(40)) | set(range(60, 70)),
          'chr2': set(range(25)), 'chr3': set()}
PLUS = {'chr1': set(range(10)) | set(range(20, 30)), 'chr2': set(), 'chr3': set()}
MINUS = {'chr1': set(range(5, 15)), 'chr2': EXONS['chr2'], 'chr3': set()}
GTF = '''chr1\tsrc\tgene\t1\t40\t.\t+\t.\tgene_id "a";
chr1\tsrc\texon\t1\t10\t.\t+\t.\tgene_id "a"; transcript_id "a1";
chr1\tsrc\texon\t3\t8\t.\t+\t.\tgene_id "a"; transcript_id "a2";
chr1\tsrc\texon\t21\t30\t.\t+\t.\tgene_id "a";
chr1\tsrc\tgene\t6\t35\t.\t-\t.\tgene_id "b";
chr1\tsrc\texon\t6\t15\t.\t-\t.\tgene_id "b";
chr1\tsrc\texon\t61\t70\t.\t.\t.\tgene_id "unknown";
chr2\tsrc\texon\t1\t8\t.\t-\t.\tgene_id "c";
chr2\tsrc\texon\t19\t25\t.\t-\t.\tgene_id "c";
absent\tsrc\texon\t1\t10\t.\t+\t.\tgene_id "ignored";
chr3\tsrc\texon\t3\t5\t.\t.\t.\t.
'''

def record(name, chrom='chr1', pos=1, cigar='5M', flag=0, nh=1, hi=1, xs=None):
    tags = ([] if nh is None else [f'NH:i:{nh}']) + ([] if hi is None else [f'HI:i:{hi}'])
    if xs: tags.append(f'XS:A:{xs}')
    n = sum(int(n) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar) if op in 'MIS=X')
    return [name, flag, chrom, pos, cigar, n, tags]

def write_sam(path, records):
    records = sorted(records, key=lambda r: (list(LENGTHS).index(r[2]) if r[2] != '*' else 99, r[3]))
    lines = ['@HD\tVN:1.6\tSO:coordinate'] + [f'@SQ\tSN:{c}\tLN:{n}' for c, n in LENGTHS.items()]
    for name, flag, chrom, pos, cigar, n, tags in records:
        lines.append('\t'.join(map(str, [name, flag, chrom, pos, 60, cigar, '*', 0, 0,
                                        'A'*n if n else '*', 'I'*n if n else '*', *tags])))
    path.write_text('\n'.join(lines)+'\n')

def expected(records, copies=1):
    depth = {c: [[0, 0, 0] for _ in range(n)] for c, n in LENGTHS.items()}
    classes = [[0, 0] for _ in range(3)]
    counts = Counter()
    multimap = set()
    introns = Counter()
    quadrants = [0]*4
    xs_quadrants = [0]*4
    for input_id in range(copies):
        for name, flag, chrom, start, cigar, n, tags in records:
            if flag & (256 | 2048): continue
            counts['raw'] += 1
            if flag & 4: continue
            counts['mapped'] += 1
            vals = {t[:2]: t[5:] for t in tags}
            if 'NH' not in vals: counts['missing'] += 1
            if int(vals.get('NH', 1)) > 1:
                counts['multimap_records'] += 1
                multimap.add((input_id, name, flag & 192))
                continue
            counts['unique'] += 1
            matches, deletions = [], []
            pos = start-1
            for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
                length = int(length)
                if op in 'M=X': matches += range(pos, pos+length)
                if op == 'D': deletions += range(pos, pos+length)
                if op == 'N': introns[length] += 1
                if op in 'M=XDN': pos += length
            if 'N' in cigar: counts['spliced'] += 1
            for pos in matches:
                depth[chrom][pos][0] += 1
                depth[chrom][pos][1] += not (flag & (512 | 1024))
                depth[chrom][pos][2] += 1
            for pos in deletions: depth[chrom][pos][2] += 1
            bases = [0, 0, 0]
            for pos in matches:
                category = 0 if pos in EXONS[chrom] else 1 if pos in BODIES[chrom] else 2
                bases[category] += 1
            if matches: classes[max(range(3), key=lambda i: bases[i])][0] += 1
            for i in range(3): classes[i][1] += bases[i]
            plus = bool(set(matches) & PLUS[chrom])
            minus = bool(set(matches) & MINUS[chrom])
            if plus and minus: counts['ambiguous'] += 1
            elif plus or minus:
                same = bool(flag & 16) == minus
                quadrants[(2 if flag & 128 else 0) + same] += 1
            if 'XS' in vals:
                same = ('-' if flag & 16 else '+') == vals['XS']
                xs_quadrants[(2 if flag & 128 else 0) + same] += 1
    return depth, classes, counts, len(multimap), introns, quadrants, xs_quadrants

def run(root, sam, annotation=None, args=(), out='out', copies=1, preexec=None, success=True):
    output = root/out
    cmd = [BINARY, '--rna', '-p', str(root/'target.bed'), '-o', str(output), '--flank', '0']
    if annotation: cmd += ['--annotation', str(annotation)]
    result = subprocess.run(cmd+list(args)+[str(sam)]*copies, capture_output=True, text=True, preexec_fn=preexec)
    assert (result.returncode == 0) == success, (cmd, result.returncode, result.stderr)
    assert not list(output.glob('.xamdst-rna-dedup-*')), list(output.iterdir())
    assert not [p for p in output.iterdir() if '.tmp.' in p.name or '.bak.' in p.name]
    return output, result

def check(output, records, annotation=True, copies=1, summary=False):
    assert {p.name for p in output.iterdir()} == RNA_FILES - ({'depth.tsv.gz'} if summary else set())
    report = json.loads((output/'coverage.report.json').read_text())
    depth, classes, counts, multi, introns, q, xq = expected(records, copies)
    rna = report['rna']
    assert report['total']['raw_reads'] == counts['raw']
    assert report['total']['mapped_reads'] == counts['mapped']
    assert rna['unique_mapping_reads'] == counts['unique']
    assert rna['multimap_records'] == counts['multimap_records']
    assert rna['multimap_reads'] == multi
    assert abs(rna['unique_mapping_rate'] - counts['unique']/(counts['unique']+multi)) < 1e-6
    assert abs(rna['multimap_read_fraction'] - multi/(counts['unique']+multi)) < 1e-6
    assert rna['nh_missing_records'] == counts['missing']
    assert rna['spliced_reads'] == counts['spliced']
    assert rna['introns']['count'] == sum(introns.values())
    splice_rows=[line.split() for line in gzip.open(output/'splice.tsv.gz','rt') if not line.startswith('#')]
    assert {int(row[0]):int(row[1]) for row in splice_rows} == dict(introns)
    for row in splice_rows: assert abs(float(row[2])-int(row[1])/sum(introns.values())) < 1e-6
    if not summary:
        lines = gzip.open(output/'depth.tsv.gz', 'rt').readlines()
        observed = [(v[0], int(v[1])-1, list(map(int, v[2:]))) for line in lines
                    if not line.startswith('#') and (v := line.split())]
        assert len(observed) == sum(LENGTHS.values())
        for chrom, pos, values in observed: assert values == depth[chrom][pos], (chrom,pos,values,depth[chrom][pos])
    expected_sum = sum(row[0] for c in depth.values() for row in c)
    assert abs(report['target']['average_depth'] - expected_sum/sum(LENGTHS.values())) < 1e-6
    assert '[RNA] Enabled\ttrue\n' in (output/'coverage.report').read_text()
    strand = rna['strand']
    assert strand['xs_reads'] == sum(xq)
    assert [strand[k] for k in ('q_read1_antisense','q_read1_sense','q_read2_antisense','q_read2_sense')] == xq
    if annotation:
        for i, name in enumerate(('exonic', 'intronic', 'intergenic')):
            assert rna['distribution'][name]['reads'] == classes[i][0]
            assert rna['distribution'][name]['bases'] == classes[i][1]
        distribution=[line.split() for line in (output/'distribution.tsv').read_text().splitlines()
                      if not line.startswith('#')]
        for i,row in enumerate(distribution[:3]):
            assert [int(row[1]),int(row[3])] == classes[i]
        assert int(distribution[3][1]) == sum(c[0] for c in classes)
        assert int(distribution[3][3]) == sum(c[1] for c in classes)
        assert rna['distribution']['exon_bases'] == sum(map(len, EXONS.values()))
        assert rna['distribution']['intron_bases'] == sum(len(BODIES[c]-EXONS[c]) for c in LENGTHS)
        assert rna['distribution']['annotated_span'] == sum(LENGTHS.values())
        assert strand['evidence'] == 'annotation_exon'
        assert strand['effective_reads'] == sum(q)
        assert strand['ambiguous_reads'] == counts['ambiguous']
        assert [strand[k] for k in ('annotation_read1_antisense','annotation_read1_sense',
                                   'annotation_read2_antisense','annotation_read2_sense')] == q
    else:
        assert rna['distribution'] is None
        assert strand['evidence'] == 'none' and strand['effective_reads'] == 0
        assert strand['inference'] == 'insufficient_data'
    return report

def snapshot(path): return {p.name: p.read_bytes() for p in path.iterdir() if p.is_file()}

with tempfile.TemporaryDirectory(prefix='xamdst-rna-') as tmp:
    root = Path(tmp)
    (root/'target.bed').write_text(''.join(f'{c}\t0\t{n}\n' for c, n in LENGTHS.items()))
    annotation = root/'genes.gtf'; annotation.write_text(GTF)
    sam = root/'rna.sam'
    records = [record('star_default', nh=3, hi=2),
               record('star_default', pos=12, flag=256, nh=3, hi=1),
               record('star_default', chrom='chr2', flag=256, nh=3, hi=3),
               record('allbest', pos=2, nh=3, hi=2),
               record('allbest', pos=15, nh=3, hi=1),
               record('allbest', chrom='chr2', nh=3, hi=3),
               record('hi_missing', nh=2, hi=None),
               record('hi_missing', chrom='chr2', nh=2, hi=None),
               record('pair', flag=65, nh=2), record('pair', flag=129, nh=2),
               record('pair', chrom='chr2', flag=65, nh=2, hi=2),
               record('pair', chrom='chr2', flag=129, nh=2, hi=2),
               record('fallback', pos=16, nh=None, hi=None),
               record('splice', cigar='5M10N5M', xs='+'),
               record('opposite', pos=6, xs='-'),
               record('opposite_blocks', cigar='3M8N3M',xs='+'),
               record('annotation_differs_xs', flag=81, cigar='3M', xs='-'),
               record('mate2', flag=129, cigar='3M', xs='+'),
               record('exon_intron_tie', pos=26, cigar='10M'),
               record('intron_intergenic_tie', pos=36, cigar='10M'),
               record('unknown_strand', pos=61),
               record('unannotated_chrom', chrom='chr3'),
               record('cigar', chrom='chr2', pos=4, cigar='2S2=1I3X2D4M2N2M', flag=1024),
               record('unmapped', chrom='*', pos=0, flag=4, cigar='*', nh=None, hi=None),
               record('supplement', flag=2048)]
    write_sam(sam, records)
    baseline, result = run(root, sam, annotation)
    assert 'lack NH' in result.stderr
    check(baseline, records)
    for name, options in [('parallel', ['--compute-threads','4']),
                          ('spill', ['--rna-dedup-mem','128']),
                          ('parallel-spill', ['--compute-threads','4','--rna-dedup-mem','128'])]:
        output, _ = run(root,sam,annotation,args=options,out=name)
        check(output,records)
        assert snapshot(output) == snapshot(baseline), name
    output, _ = run(root,sam,annotation,args=['--summary-only'],out='summary')
    check(output,records,summary=True)
    # Successful summary replacement removes stale depth; failures keep it.
    before = snapshot(baseline)
    output, _ = run(root,sam,annotation,args=['--summary-only'],out='out')
    check(output,records,summary=True)
    output, _ = run(root,sam,annotation,out='out')
    assert snapshot(output) == before
    output, _ = run(root,sam,annotation,copies=2,out='multi-input',args=['--rna-dedup-mem','128'])
    check(output,records,copies=2)
    output, _ = run(root,sam,out='no-annotation')
    check(output,records,annotation=False)
    gz = root/'genes.gtf.gz'; gz.write_bytes(gzip.compress(GTF.encode()))
    output, _ = run(root,sam,gz,out='gzip')
    check(output,records)
    gff = root/'genes.gff2'
    # GFF2 transcript grouping, with explicit gene spans retained.
    gff.write_text(GTF.replace('gene_id ', 'gene '))
    output, _ = run(root,sam,gff,out='gff2')
    check(output,records)
    transcript_gff=root/'transcripts.gff2'
    transcript_gff.write_text(re.sub(r' transcript_id "[^"]+";', '', GTF).replace('gene_id ', 'Transcript '))
    output,_=run(root,sam,transcript_gff,out='gff2-transcript')
    check(output,records)
    if shutil.which('samtools'):
        reference=root/'reference.fa'
        reference.write_text(''.join(f'>{chrom}\n'+ 'A'*n+'\n' for chrom,n in LENGTHS.items()))
        subprocess.run(['samtools','faidx',str(reference)],check=True)
        for suffix,options in [('bam',['-b']),('cram',['-C','-T',str(reference)])]:
            converted=root/('rna.'+suffix)
            subprocess.run(['samtools','view',*options,'-o',str(converted),str(sam)],check=True)
            output,_=run(root,converted,annotation,out=suffix,
                         args=['--reference',str(reference),'--compute-threads','4'])
            check(output,records)
    for tag in ['NH:i:0','NH:i:-1','NH:Z:1','NH:f:1','HI:i:0','HI:i:-1','HI:Z:1','HI:i:2']:
        bad = record('bad'); bad[6] = [t for t in bad[6] if t[:2] != tag[:2]]+[tag]
        broken = root/'bad.sam'; write_sam(broken, records+[bad])
        for workers in ('1','4'):
            _, failure = run(root,broken,annotation,out='out',args=['--compute-threads',workers,'--rna-dedup-mem','128'],success=False)
            assert 'invalid NH/HI' in failure.stderr
            assert snapshot(baseline) == before
    for text in [GTF.replace('\texon\t1\t10', '\texon\t0\t10'),
                 GTF.replace('\t1\t40', '\t1\t999'),
                 GTF.replace('gene_id "a";', 'gene_id "a;'), 'garbage\n']:
        bad_gtf=root/'bad.gtf'; bad_gtf.write_text(text)
        run(root,sam,bad_gtf,out='out',success=False)
        assert snapshot(baseline) == before
    truncated=root/'truncated.gtf.gz'; truncated.write_bytes(gz.read_bytes()[:-12])
    run(root,sam,truncated,out='out',success=False)
    assert snapshot(baseline) == before
    def limit_file_size():
        signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
        resource.setrlimit(resource.RLIMIT_FSIZE,(256,256))
    spill_failure = root/'spill-failure.sam'
    write_sam(spill_failure, [record('long'+'q'*180,nh=2), record('other'+'q'*180,nh=2)])
    _, failure=run(root,spill_failure,annotation,out='out',args=['--rna-dedup-mem','128'],
                   preexec=limit_file_size,success=False)
    assert 'RNA dedup partition' in failure.stderr, failure.stderr
    assert snapshot(baseline) == before
    unique_failure=root/'unique-failure.sam'
    write_sam(unique_failure,[record('unique-only')])
    run(root,unique_failure,annotation,out='out',preexec=limit_file_size,success=False)
    assert snapshot(baseline) == before
    # Reusing an RNA directory in DNA mode removes the two RNA reports only
    # after successful commit. A protected input in that namespace is refused.
    mode_dir=root/'mode-switch'; shutil.copytree(baseline,mode_dir)
    dna_cmd=[BINARY,'-p',str(root/'target.bed'),'-o',str(mode_dir),'--flank','0']
    protected=mode_dir/'distribution.tsv'
    protected.write_text(sam.read_text())
    protected_before=snapshot(mode_dir)
    failure=subprocess.run(dna_cmd+[str(protected)],capture_output=True,text=True)
    assert failure.returncode != 0 and 'protected input' in failure.stderr
    assert snapshot(mode_dir) == protected_before
    subprocess.run(dna_cmd+[str(sam)],check=True)
    assert {p.name for p in mode_dir.iterdir()} == DNA_FILES
    dna=json.loads((mode_dir/'coverage.report.json').read_text())
    assert 'rna' not in dna and dna['total']['raw_reads'] == expected(records)[2]['raw']
    # Force several merge levels with repeated read names on different references.
    many_multi=[record(f'multi{i}',chrom=c,nh=3,hi=j+1) for i in range(513)
                for j,c in enumerate(LENGTHS)]
    stress=root/'multi-stress.sam'; write_sam(stress,many_multi)
    output,_=run(root,stress,annotation,out='multi-stress',args=['--rna-dedup-mem','4096','--compute-threads','4'])
    check(output,many_multi)
    # More than 1000 informative read ends are required for a strand call.
    for label, flags in [('fr-firststrand',(81,129)),('fr-secondstrand',(65,145)),
                         ('unstranded',(65,81,129,145))]:
        many=[record(f'strand{i}',flag=flags[i%len(flags)],cigar='3M',xs='-') for i in range(1004)]
        strand_sam=root/'strand.sam'; write_sam(strand_sam,many)
        output,_=run(root,strand_sam,annotation,out=label)
        assert check(output,many)['rna']['strand']['inference'] == label
        if label == 'fr-firststrand':
            output,_=run(root,strand_sam,out='xs-alone')
            check(output,many,annotation=False)
print('xamdst RNA oracle and rollback tests passed')
