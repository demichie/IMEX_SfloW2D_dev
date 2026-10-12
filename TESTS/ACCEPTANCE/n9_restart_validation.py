"""Validate the v5 diagnostic restart against continuous and pre-change runs.

The prior dated N9 artifact remains historical evidence of the v4 defect.
This driver records fresh full/split runs, exact CSV history, unchanged native
fields and unchanged canonical checkpoint records, without new tolerances.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from zipfile import ZipFile

import numpy as np
import n9_long_cases as n9

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'TEST_BINARY_RESTART'))
from test_runout_restart import records


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--repo', type=Path, required=True)
    cli.add_argument('--build-root', type=Path, required=True)
    cli.add_argument('--prior-evidence', type=Path, required=True)
    cli.add_argument('--output', type=Path, required=True)
    cli.add_argument('--etna-data', type=Path, help='Reuse verified preprocessing rasters, without modifying them')
    args = cli.parse_args()
    repo, root = args.repo.resolve(), args.output.resolve()
    assert not root.exists()
    root.mkdir()
    prior = json.loads(args.prior_evidence.read_text())
    assert prior['summary']['physical_continuations_passed'] == 8
    source_hashes = {str(p.relative_to(repo)):n9.digest(p) for p in (repo/'src').rglob('*')
                     if p.suffix in ('.f90','.inc')}
    source_hashes_by_basename = {Path(name).name: digest for name, digest in source_hashes.items()}
    evidence = {'schema_version':1, 'started_utc':datetime.now(timezone.utc).isoformat(),
        'source_base_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
        'source_sha256':source_hashes, 'prior_evidence_sha256':n9.digest(args.prior_evidence),
        'runner_sha256':n9.digest(Path(__file__)), 'runs':[], 'thread_pairs':[],
        'scientific_baseline_approved':False, 'solver_tolerances_changed':False,
        'contract':'exact same-build full/restart CSVs and checkpoints; native fields and canonical records exact against pre-change physical baseline'}
    for name, value in prior['source_sha256'].items():
        basename = Path(name).name
        if basename not in ('IMEX_SfloW2D.f90','inpout_2d.f90','runtime_2d.f90') and basename in source_hashes_by_basename:
            assert source_hashes_by_basename[basename] == value, ('unexpected numerical source change',name)
    data = root/'ETNA-data'
    data.mkdir(); (data/'DEM').mkdir()
    for name in ('IMEX_SfloW2D.inp','create_input_ellipsoid.py','input_ellipsoid.py'):
        shutil.copy2(repo/'EXAMPLES/EXAMPLE_ETNA'/name,data/name)
    if args.etna_data:
        for name in ('new_dem.asc','init_volume.asc','ellipsoid.asc','Etna2014_crop.asc'):
            source = args.etna_data/'DEM'/name
            expected = prior['ETNA_fixture']['files']['DEM/'+name]
            assert n9.digest(source) == expected
            (data/'DEM'/name).symlink_to(source.resolve())
    else:
        with ZipFile(repo/'EXAMPLES/EXAMPLE_ETNA/DEM/Etna2014_crop.zip') as archive:
            with archive.open('Etna2014_crop.asc') as src, (data/'DEM/Etna2014_crop.asc').open('xb') as dst:
                shutil.copyfileobj(src,dst)
        subprocess.run([sys.executable,str(repo/'TESTS/ACCEPTANCE/n9_etna_preprocess.py')],
                       cwd=data,check=True,timeout=300)
        for name in ('new_dem.asc','init_volume.asc','ellipsoid.asc'):
            assert n9.digest(data/'DEM'/name) == prior['ETNA_fixture']['files']['DEM/'+name]

    def save():
        (root/'evidence.json').write_text(json.dumps(evidence,indent=2,allow_nan=False)+'\n')

    def run(exe, directory, threads, logname):
        start = time.monotonic()
        with (directory/logname).open('w') as log:
            proc = subprocess.run([str(exe)],cwd=directory,timeout=1800,
                env={**os.environ,'OMP_NUM_THREADS':str(threads),'OMP_DYNAMIC':'FALSE'},
                stdout=log,stderr=subprocess.STDOUT)
        output = (directory/logname).read_text()
        team = re.search(r'number of threads used\s+(\d+)',output)
        assert proc.returncode == 0 and team and int(team[1]) == threads, (directory,output[-2000:])
        assert not (directory/'runout.txt').exists()
        return {'returncode':proc.returncode,'actual_threads':int(team[1]),
                'seconds':time.monotonic()-start, 'log_sha256':n9.digest(directory/logname),
                'input_sha256':n9.digest(directory/'IMEX_SfloW2D.inp')}

    def job(case, profile, threads):
        exe = args.build_root/profile/'src/IMEX_SfloW2D'
        tag = f'{case}-{profile}-{threads}'
        full, split = root/(tag+'-full'), root/(tag+'-split')
        canonical = (repo/f'EXAMPLES/EXAMPLE_{case}/IMEX_SfloW2D.inp').read_text()
        for directory in (full,split):
            directory.mkdir(); (directory/'DEM').mkdir()
            names = ('new_dem.asc','init_volume.asc','ellipsoid.asc') if case=='ETNA' else ('Etna2014_crop.asc',)
            for name in names:
                (directory/'DEM'/name).symlink_to((data/'DEM'/name).resolve())
            (directory/'IMEX_SfloW2D.inp').write_text(canonical)
        full_run = run(exe,full,threads,'full.log')
        half = 100 if case=='ETNA' else 1800
        (split/'IMEX_SfloW2D.inp').write_text(n9.changed_parameter(canonical,'T_END',str(half)))
        first_run = run(exe,split,threads,'first.log')
        shutil.copy2(split/'restart.bin',split/'split_checkpoint.bin')
        resumed = n9.changed_parameter(canonical,'RESTART','T')
        if case=='ETNA':
            resumed = n9.changed_parameter(resumed,'RESTART_FILES','"restart.bin"')
        else:
            resumed += '\n&RESTART_PARAMETERS\n RESTART_FILES="restart.bin",\n/\n'
        (split/'IMEX_SfloW2D.inp').write_text(resumed)
        second_run = run(exe,split,threads,'resumed.log')
        old = next(r for r in prior['runs'] if r['case']==case and r['profile']==profile
                   and r['requested_threads']==threads)
        baseline = Path(old['reference_directory'])
        return {'case':case,'profile':profile,'requested_threads':threads,
            'executable_sha256':n9.digest(exe), 'full_directory':str(full),
            'split_directory':str(split),'baseline_directory':str(baseline),
            'split_time_s':half,'full_run':full_run,'first_run':first_run,'second_run':second_run}

    reader = n9.NetCDF()
    jobs = [(c,p,t) for p in ('historical_optimized','strict_debug')
            for c in ('ETNA','ETNA_LAVA') for t in (1,4)]
    save()
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = {pool.submit(job,*job_args):job_args for job_args in jobs}
        for future in as_completed(futures):
            try:
                result = future.result()
                full, split, baseline = (Path(result[k]) for k in
                                         ('full_directory','split_directory','baseline_directory'))
                name = 'exampleEtna' if result['case']=='ETNA' else 'exampleEtna_lava'
                left, lm = reader.read(full/(name+'.nc'))
                right, rm = reader.read(split/(name+'.nc'))
                original, om = reader.read(baseline/(name+'.nc'))
                result['native_fields_restart_exact'] = lm==rm and all(np.array_equal(left[k],right[k]) for k in left)
                result['native_fields_prechange_exact'] = lm==om and all(np.array_equal(left[k],original[k]) for k in left)
                result['field_differences_prechange'] = {k:float(np.max(np.abs(left[k]-original[k]))) for k in left}
                result['decoded_sha256'] = {k:hashlib.sha256(v.tobytes()).hexdigest() for k,v in left.items()}
                full_records, split_records, old_records = (records(d/'restart.bin') for d in (full,split,baseline))
                result['checkpoint_v5_restart_exact'] = full_records==split_records
                result['canonical_records_prechange_exact'] = full_records[2:-2]==old_records[2:]
                result['checkpoint_sha256'] = n9.digest(full/'restart.bin')
                result['diagnostic_csv'] = n9.diagnostic_comparison(full,split,name)
                result['CSV_prechange_full_exact'] = all(n9.digest(full/(name+'_'+s+'.csv'))==
                    n9.digest(baseline/(name+'_'+s+'.csv')) for s in ('runout','mass_center'))
                rasters = [{p.name:n9.digest(p) for p in d.glob('*.asc')} for d in (full,split,baseline)]
                result['ESRI_restart_and_prechange_exact'] = rasters[0]==rasters[1]==rasters[2]
                result['ESRI_files_compared'] = len(rasters[0])
                result['prior_mass_observations_retained'] = True
                result['budget_scope'] = 'ETNA full 200 s; lava saved snapshots through 3180 s, before boundary contact'
                gates = ('native_fields_restart_exact','native_fields_prechange_exact','checkpoint_v5_restart_exact',
                         'canonical_records_prechange_exact','CSV_prechange_full_exact','ESRI_restart_and_prechange_exact')
                result['status'] = 'pass' if all(result[k] for k in gates) and all(
                    item['status']=='pass' and item['full_sha256']==item['resumed_sha256']
                    for item in result['diagnostic_csv'].values()) else 'fail'
                print('FINISHED',result['case'],result['profile'],result['requested_threads'],result['status'],flush=True)
            except Exception as error:
                c,p,t = futures[future]
                result = {'case':c,'profile':p,'requested_threads':t,'status':'fail','error':repr(error)}
                print('FAILED',c,p,t,repr(error),flush=True)
            evidence['runs'].append(result); save()
    for case in ('ETNA','ETNA_LAVA'):
        for profile in ('historical_optimized','strict_debug'):
            pair = [r for r in evidence['runs'] if r['case']==case and r['profile']==profile]
            passed = len(pair)==2 and all(r['status']=='pass' for r in pair) and (
                pair[0]['decoded_sha256']==pair[1]['decoded_sha256'] and
                pair[0]['checkpoint_sha256']==pair[1]['checkpoint_sha256'])
            evidence['thread_pairs'].append({'case':case,'profile':profile,'status':'pass' if passed else 'fail'})
    evidence['summary'] = {'passed':sum(r['status']=='pass' for r in evidence['runs']),
                           'failed':sum(r['status']!='pass' for r in evidence['runs'])}
    save()
    assert evidence['summary']=={'passed':8,'failed':0} and all(
        p['status']=='pass' for p in evidence['thread_pairs'])


if __name__=='__main__':
    main()
