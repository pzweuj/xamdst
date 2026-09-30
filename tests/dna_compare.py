#!/usr/bin/env python3
"""Compare DNA outputs and interleaved timing against an archived Git revision.

Not part of make test: at the default scale this generates about 500 MB of
temporary SAM data and runs serial/parallel, full/summary measurements.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import statistics
import subprocess
import tempfile
import time

ROOT = Path(__file__).resolve().parent.parent
FLAGS = '-O2 -g -std=c99 -Wall -Wextra -Wpedantic -Wformat=2 -Werror'

def invoke(args, **kw):
    return subprocess.run(list(map(str,args)),check=True,**kw)

def compare_outputs(left,right):
    assert {p.name for p in left.iterdir()} == {p.name for p in right.iterdir()}
    for p in left.iterdir():
        other=right/p.name
        if p.name == 'coverage.report.json':
            a=json.loads(p.read_text()); b=json.loads(other.read_text())
            for key in ('schema_version','version'): a.pop(key,None); b.pop(key,None)
            assert a == b, p.name
        elif p.name == 'coverage.report':
            a=[line for line in p.read_text().splitlines() if not line.startswith('## Version')]
            b=[line for line in other.read_text().splitlines() if not line.startswith('## Version')]
            assert a == b, p.name
        elif p.name in ('cumu.plot','insert.plot'):
            a=[line.split() for line in p.read_text().splitlines()]
            b=[line.split() for line in other.read_text().splitlines()]
            assert len(a) == len(b)
            for old,new in zip(a,b):
                assert old[:3] == new[:3], (p.name,old,new)
                assert int(new[3]) == int(old[3])+int(old[1]), (p.name,old,new)
                assert abs(float(new[4])-float(old[4])-float(old[2])) < 2e-9
        else:
            assert p.read_bytes() == other.read_bytes(), p.name

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('binary',type=Path)
    parser.add_argument('--baseline',default='af559ba')
    parser.add_argument('--rounds',type=int,default=7)
    parser.add_argument('--scale',type=int,default=1)
    parser.add_argument('--pin',action='store_true',help='pin each child to a fixed set of CPUs')
    parser.add_argument('--datasets',nargs='+',choices=['panel','broad'],default=['panel','broad'])
    parser.add_argument('--modes',nargs='+',choices=['full1','summary1','full4','summary4'],
                        default=['full1','summary1','full4','summary4'])
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.rounds < 1 or args.scale < 1: parser.error('rounds and scale must be positive')
    candidate=args.binary.resolve()
    env=dict(os.environ,CFLAGS=FLAGS)
    result={'baseline_revision':args.baseline,'candidate_version':invoke([candidate,'--version'],capture_output=True,text=True).stdout.strip(),
            'platform':platform.platform(),'compiler':invoke(['cc','--version'],capture_output=True,text=True).stdout.splitlines()[0],
            'flags':FLAGS,'rounds':args.rounds,'scale':args.scale,'cases':[]}
    result['candidate_sha256']=hashlib.sha256(candidate.read_bytes()).hexdigest()
    cpus=sorted(os.sched_getaffinity(0))[-4:] if args.pin else []
    result['affinity_cpus']=cpus
    with tempfile.TemporaryDirectory(prefix='xamdst-dna-compare-') as temp:
        work=Path(temp); source=work/'baseline'; source.mkdir()
        archive=work/'baseline.tar'
        with archive.open('wb') as f:
            invoke(['git','-c',f'safe.directory={ROOT}','archive',args.baseline],cwd=ROOT,stdout=f)
        invoke(['tar','xf',archive,'-C',source])
        invoke(['make','-j2'],cwd=source,env=env,stdout=subprocess.DEVNULL)
        baseline=source/'xamdst'
        result['baseline_version']=invoke([baseline,'--version'],capture_output=True,text=True).stdout.strip()
        for dataset,length,depth in [('panel',40000,500),('broad',2000000,10)]:
            if dataset not in args.datasets: continue
            depth *= args.scale
            data=work/dataset
            invoke(['python3',ROOT/'tests/generate_panel.py','--outdir',data,'--length',length,'--depth',depth],stdout=subprocess.DEVNULL)
            for workers,summary in [(1,False),(1,True),(4,False),(4,True)]:
                if ('summary' if summary else 'full')+str(workers) not in args.modes: continue
                options=['--compute-threads',str(workers),'-p',str(data/'panel.bed')]
                if summary: options+=['--summary-only']
                binaries={'baseline':baseline,'candidate':candidate}
                preexec=(lambda: os.sched_setaffinity(0,cpus[:workers])) if cpus else None
                times={name:[] for name in binaries}
                cpu_times={name:[] for name in binaries}
                # Warm both executables and the input page cache, then verify
                # every report before measuring. Both use identical options.
                for name,binary in binaries.items():
                    invoke([binary,*options,'-o',work/(name+'-out'),data/'panel.sam'],stdout=subprocess.DEVNULL,preexec_fn=preexec)
                compare_outputs(work/'baseline-out',work/'candidate-out')
                orders=(['baseline','candidate'], ['candidate','baseline'],
                        ['candidate','baseline'], ['baseline','candidate'],
                        ['candidate','baseline'], ['baseline','candidate'],
                        ['baseline','candidate'], ['candidate','baseline'],
                        ['candidate','baseline'], ['baseline','candidate'])
                for index in range(args.rounds):
                    order=orders[index%len(orders)]
                    for name in order:
                        usage=resource.getrusage(resource.RUSAGE_CHILDREN)
                        cpu_start=usage.ru_utime+usage.ru_stime
                        start=time.perf_counter()
                        invoke([binaries[name],*options,'-o',work/(name+'-out'),data/'panel.sam'],stdout=subprocess.DEVNULL,preexec_fn=preexec)
                        times[name].append(time.perf_counter()-start)
                        usage=resource.getrusage(resource.RUSAGE_CHILDREN)
                        cpu_times[name].append(usage.ru_utime+usage.ru_stime-cpu_start)
                medians={name:statistics.median(values) for name,values in times.items()}
                change=(medians['candidate']/medians['baseline']-1)*100
                cpu_medians={name:statistics.median(values) for name,values in cpu_times.items()}
                cpu_change=(cpu_medians['candidate']/cpu_medians['baseline']-1)*100
                item={'dataset':dataset,'reference_length':length,'read_length':100,
                      'records':len(range(0,length-100+1,20))*depth,'input_bytes':(data/'panel.sam').stat().st_size,
                      'compute_threads':workers,'summary_only':summary,'depth_and_counts_equivalent':True,
                      'expected_report_changes':['version/schema_version','plot cumulative > to >='],
                      'seconds':times,'median_seconds':medians,'change_percent':change,
                      'cpu_seconds':cpu_times,'median_cpu_seconds':cpu_medians,'cpu_change_percent':cpu_change}
                result['cases'].append(item)
                print(f'{dataset} threads={workers} summary={summary}: baseline={medians["baseline"]:.4f}s candidate={medians["candidate"]:.4f}s change={change:+.2f}% cpu_change={cpu_change:+.2f}%',flush=True)
                args.output.parent.mkdir(parents=True,exist_ok=True)
                args.output.write_text(json.dumps(result,indent=2)+'\n')
    return 0

if __name__ == '__main__': raise SystemExit(main())
