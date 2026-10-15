"""Exercise restartable runout and safe CSV histories using the actual solver.

All inputs and outputs are temporary copies. No user simulation is modified.
Use --output to retain the cases and a machine-readable verification manifest.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import tempfile


def records(path):
    result = []
    with path.open('rb') as stream:
        while marker := stream.read(4):
            assert len(marker) == 4
            size = struct.unpack('<i', marker)[0]
            assert size >= 0
            payload = stream.read(size)
            assert len(payload) == size and stream.read(4) == marker
            result.append(payload)
    assert result[0] == b'IMEX_SFLOW2D_RST'
    assert struct.unpack('<i', result[1])[0] in (4, 5)
    return result


def write_records(path, payloads):
    with path.open('wb') as stream:
        for payload in payloads:
            marker = struct.pack('<i', len(payload))
            stream.write(marker+payload+marker)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tree_hashes(directory):
    return {str(path.relative_to(directory)): sha(path)
            for path in directory.rglob('*') if path.is_file()}


def parameter(text, name, value):
    text, count = re.subn(r'(?mi)^(\s*'+name+r'\s*=)[^\n]*',
                         lambda m: m[1]+' '+value+',', text)
    assert count == 1, (name, count)
    return text


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('executable', type=Path)
    cli.add_argument('--netcdf-helper', type=Path, required=True)
    cli.add_argument('--threads', type=int, default=int(os.environ.get('OMP_NUM_THREADS', 1)))
    cli.add_argument('--output', type=Path)
    args = cli.parse_args()
    repo = Path(__file__).resolve().parents[2]
    root = args.output or Path(tempfile.mkdtemp(prefix='imex-runout-restart-'))
    if args.output:
        root.mkdir()
    exe = args.executable.resolve()
    fixture = repo/'TESTS/TEST_2D'
    text = (fixture/'IMEX_sfloW2D.inp').read_text()
    text = parameter(text, 'OUTPUT_RUNOUT_FLAG', 'T')
    text = text.replace(' OUTPUT_RUNOUT_FLAG= T,', ' OUTPUT_RUNOUT_FLAG= T,\n SERIAL_FLAG=F,')
    text = parameter(text, 'T_END', '0.015')
    text = parameter(text, 'DT_OUTPUT', '0.015')
    text = parameter(text, 'VERBOSE_LEVEL', '0')
    text += '\n&RUNOUT_PARAMETERS\n DT_RUNOUT=0.005,\n EPS_STOP=0,\n/\n'
    checks = []
    manifest = {'executable_sha256': sha(exe), 'threads': args.threads,
                'directory': str(root), 'checks': checks}

    def prepare(name, source=text):
        case = root/name
        shutil.copytree(fixture, case)
        (case/'IMEX_SfloW2D.inp').write_text(source)
        return case

    def run(case, logname, expected_error=None):
        with (case/logname).open('w') as log:
            proc = subprocess.run([str(exe)], cwd=case, timeout=120,
                env={**os.environ, 'OMP_NUM_THREADS':str(args.threads), 'OMP_DYNAMIC':'FALSE'},
                stdout=log, stderr=subprocess.STDOUT)
        output = (case/logname).read_text()
        if expected_error:
            assert proc.returncode != 0 and expected_error.lower() in output.lower(), (case, output[-3000:])
        else:
            assert proc.returncode == 0, (case, output[-3000:])
            assert re.search(r'number of threads used\s+'+str(args.threads)+r'\b', output)
            assert not (case/'runout.txt').exists()
        return output

    def require_output(case, filename):
        path = case/'OUTPUT'/'TEST_2D'/'diagnostics'/filename
        if path.is_file():
            return path
        details = [f'Missing expected output: {path}',
                   f'Files created: {sorted(str(p.relative_to(case)) for p in case.rglob("*") if p.is_file())}']
        inp = case/'IMEX_SfloW2D.inp'
        if inp.is_file():
            relevant_input = [line for line in inp.read_text().splitlines()
                              if any(key in line.upper() for key in
                                     ('RUN_NAME', 'OUTPUT_RUNOUT_FLAG', 'SERIAL_FLAG',
                                      'T_END', 'DT_OUTPUT', 'RUNOUT_PARAMETERS',
                                      'DT_RUNOUT', 'EPS_STOP'))]
            details.extend(['Effective input parameters:', *relevant_input])
        for logname in ('full.log', 'first.log', 'resumed.log'):
            logfile = case/logname
            if logfile.is_file():
                details.extend([f'Last lines of {logname}:',
                                *logfile.read_text(errors='replace').splitlines()[-60:]])
        raise AssertionError('\n'.join(details))

    def continuation(case, end=0.03):
        inp = case/'IMEX_SfloW2D.inp'
        value = parameter(inp.read_text(), 'T_END', str(end))
        inp.write_text(parameter(value, 'RESTART_FILES', '"restart.bin"'))

    def allow_output_reuse(case):
        inp = case/'IMEX_SfloW2D.inp'
        source, count = re.subn(r'(?mi)^&RUN_PARAMETERS\s*$',
                                '&RUN_PARAMETERS\n OUTPUT_OVERWRITE=T,', inp.read_text())
        assert count == 1
        inp.write_text(source)

    def checkpoint(case):
        return case/'OUTPUT'/'TEST_2D'/'restart'/'restart.bin'

    def netcdf_path(case):
        return case/'OUTPUT'/'TEST_2D'/'fields'/'TEST_2D.nc'

    def netcdf_helper(*arguments):
        proc = subprocess.run([str(args.netcdf_helper.resolve()), *map(str, arguments)],
                              capture_output=True, text=True, timeout=30)
        assert proc.returncode == 0, proc.stderr or proc.stdout

    def diagnostics(case):
        path = case/'OUTPUT'/'TEST_2D'/'diagnostics'
        path.mkdir(parents=True, exist_ok=True)
        return path

    def snapshot(case, filename):
        matches = list((case/'OUTPUT'/'TEST_2D'/'snapshots').glob('*/'+filename))
        assert len(matches) == 1, (case, filename, matches)
        return matches[0]

    def passed(name):
        checks.append({'name': name, 'status':'pass'})
        print('PASS:', name, flush=True)

    for mode in ('automatic', 'explicit'):
        source = text if mode == 'automatic' else text.replace(' EPS_STOP=0,',
                                         ' EPS_STOP=0,\n X0_RUNOUT=5,\n Y0_RUNOUT=0,')
        full = prepare(mode+'-full', parameter(source, 'T_END', '0.03'))
        split = prepare(mode+'-split', source)
        run(full, 'full.log'); run(split, 'first.log')
        shutil.copy2(checkpoint(split), split/'split_checkpoint.bin')
        continuation(split); run(split, 'resumed.log')
        assert records(checkpoint(full)) == records(checkpoint(split))
        for suffix in ('runout', 'mass_center'):
            filename = 'TEST_2D_'+suffix+'.csv'
            assert require_output(full, filename).read_bytes() == \
                require_output(split, filename).read_bytes()
        assert snapshot(full, 'TEST_2D_0002.q_2d').read_bytes() == \
            snapshot(split, 'TEST_2D_0002.q_2d').read_bytes()
        passed(mode+' origin: exact checkpoint, physical output and complete CSV history')

    source_case = root/'automatic-full'
    case = prepare('high-index-manifest', parameter(text, 'T_END', '0.045'))
    (case/'OUTPUT').mkdir()
    shutil.copytree(source_case/'OUTPUT'/'TEST_2D', case/'OUTPUT'/'TEST_2D')
    manifest_path = case/'OUTPUT'/'TEST_2D'/'metadata'/'output_manifest.csv'

    payloads = records(checkpoint(case))
    original_output_index = struct.unpack('<i', payloads[5][:4])[0]
    assert original_output_index < 10000
    checkpoint_time = struct.unpack('<4d', payloads[4])[0]
    synthetic_rows = ['snapshot_index,time_seconds'] + [
        f'{index},{checkpoint_time*index/10000:.16E}' for index in range(10000)
    ]
    manifest_path.write_text('\n'.join(synthetic_rows)+'\n')
    previous_manifest_rows = manifest_path.read_text().splitlines()
    payloads[5] = struct.pack('<i', 10000)+payloads[5][4:]
    write_records(checkpoint(case), payloads)
    continuation(case, end=0.045)
    run(case, 'high-index.log')

    manifest_rows = manifest_path.read_text().splitlines()
    assert manifest_rows[:len(previous_manifest_rows)] == previous_manifest_rows
    new_manifest_rows = manifest_rows[len(previous_manifest_rows):]
    assert new_manifest_rows, 'Restart did not append high-index manifest records'
    previous_time = float(previous_manifest_rows[-1].split(',')[1])
    appended_times = []
    for offset, row in enumerate(new_manifest_rows):
        fields = row.split(',')
        assert len(fields) == 2 and fields[0] == str(10000+offset), row
        appended_times.append(float(fields[1]))
    assert all(later > earlier for earlier, later in
               zip([previous_time]+appended_times, appended_times))
    manifest_indices = [int(row.split(',')[0]) for row in manifest_rows[1:]]
    assert len(manifest_indices) == len(set(manifest_indices))
    assert list((case/'OUTPUT'/'TEST_2D'/'snapshots').glob(
        '10000_t_*/TEST_2D_10000.q_2d'))
    passed('high-index restart: manifest and snapshot append at index 10000')

    snapshot_names = list((case/'OUTPUT'/'TEST_2D'/'snapshots').glob('10000_t_*'))
    assert len(snapshot_names) == 1
    time_label = snapshot_names[0].name.split('_t_', 1)[1]
    case = prepare('interrupted-snapshot', parameter(text, 'T_END', '0.045'))
    (case/'OUTPUT').mkdir()
    shutil.copytree(source_case/'OUTPUT'/'TEST_2D', case/'OUTPUT'/'TEST_2D')
    snapshot_dir = case/'OUTPUT'/'TEST_2D'/'snapshots'/f'0003_t_{time_label}'
    snapshot_dir.mkdir()
    partial_snapshot = snapshot_dir/'TEST_2D_0003.q_2d'
    partial_snapshot.write_bytes(b'partial snapshot from interrupted output\n')
    output_hashes_before = tree_hashes(case/'OUTPUT'/'TEST_2D')
    intact_manifest = (case/'OUTPUT'/'TEST_2D'/'metadata'/'output_manifest.csv').read_bytes()
    checkpoint_before = sha(checkpoint(case))
    continuation(case, end=0.045)
    run(case, 'snapshot-rejected.log', 'uncommitted snapshot directory')
    assert tree_hashes(case/'OUTPUT'/'TEST_2D') == output_hashes_before
    assert (case/'OUTPUT'/'TEST_2D'/'metadata'/'output_manifest.csv').read_bytes() == \
        intact_manifest
    assert sha(checkpoint(case)) == checkpoint_before

    shutil.rmtree(snapshot_dir)
    run(case, 'snapshot-recovered.log')
    assert list((case/'OUTPUT'/'TEST_2D'/'snapshots').glob(
        '0003_t_*/TEST_2D_0003.q_2d'))
    resumed_manifest_rows = (case/'OUTPUT'/'TEST_2D'/'metadata'/
                             'output_manifest.csv').read_text().splitlines()
    assert resumed_manifest_rows[-1].split(',')[0] == '3'
    passed('interrupted snapshot: preserve orphan, then resume after explicit cleanup')

    netcdf_source, netcdf_flag_count = re.subn(
        r'(?mi)^(\s*OUTPUT_CONS_FLAG\s*=[^\n]*)',
        lambda match: ' OUTPUT_NETCDF_FLAG= T,\n'+match[1], text)
    assert netcdf_flag_count == 1
    netcdf_split = prepare('netcdf-interrupted-append', netcdf_source)
    run(netcdf_split, 'netcdf-first.log')
    continuation(netcdf_split, end=0.03)
    run(netcdf_split, 'netcdf-second.log')

    netcdf_control = prepare('netcdf-control',
                             parameter(netcdf_source, 'T_END', '0.045'))
    (netcdf_control/'OUTPUT').mkdir()
    shutil.copytree(netcdf_split/'OUTPUT'/'TEST_2D',
                    netcdf_control/'OUTPUT'/'TEST_2D')
    continuation(netcdf_control, end=0.045)
    run(netcdf_control, 'netcdf-control.log')

    nc_file = netcdf_path(netcdf_split)
    netcdf_helper('partial', nc_file, '0.045')
    partial_netcdf_hash = sha(nc_file)
    manifest_before = (netcdf_split/'OUTPUT'/'TEST_2D'/'metadata'/
                       'output_manifest.csv').read_bytes()
    checkpoint_before = sha(checkpoint(netcdf_split))
    continuation(netcdf_split, end=0.045)
    run(netcdf_split, 'netcdf-recovered.log')
    assert sha(nc_file) != partial_netcdf_hash
    assert (netcdf_split/'OUTPUT'/'TEST_2D'/'metadata'/'output_manifest.csv').read_bytes() != \
        manifest_before
    assert sha(checkpoint(netcdf_split)) != checkpoint_before
    assert records(checkpoint(netcdf_control)) == records(checkpoint(netcdf_split))
    netcdf_helper('compare', netcdf_path(netcdf_control), nc_file)
    passed('interrupted NetCDF append: overwrite one uncommitted record and match full run')

    case = prepare('reject-multiple-netcdf-tails',
                   parameter(netcdf_source, 'T_END', '0.045'))
    (case/'OUTPUT').mkdir()
    shutil.copytree(netcdf_split/'OUTPUT'/'TEST_2D', case/'OUTPUT'/'TEST_2D')
    nc_file = netcdf_path(case)
    netcdf_helper('partial', nc_file, '0.06')
    netcdf_helper('partial', nc_file, '0.075')
    netcdf_hash_before = sha(nc_file)
    continuation(case, end=0.06)
    run(case, 'netcdf-rejected.log', 'more than one uncommitted time record')
    assert sha(nc_file) == netcdf_hash_before
    passed('reject multiple uncommitted NetCDF records without modifying the dataset')

    case = prepare('interrupted-manifest', parameter(text, 'T_END', '0.045'))
    (case/'OUTPUT').mkdir()
    shutil.copytree(root/'automatic-full'/'OUTPUT'/'TEST_2D', case/'OUTPUT'/'TEST_2D')
    manifest_path = case/'OUTPUT'/'TEST_2D'/'metadata'/'output_manifest.csv'
    intact_manifest = manifest_path.read_bytes()
    checkpoint_before = sha(checkpoint(case))
    manifest_path.write_bytes(intact_manifest+b'3,')
    interrupted_manifest = manifest_path.read_bytes()
    output_hashes_before = tree_hashes(case/'OUTPUT'/'TEST_2D')
    continuation(case, end=0.045)
    run(case, 'manifest-rejected.log', 'Output manifest has a malformed row')
    assert manifest_path.read_bytes() == interrupted_manifest
    assert sha(checkpoint(case)) == checkpoint_before
    assert tree_hashes(case/'OUTPUT'/'TEST_2D') == output_hashes_before

    manifest_path.write_bytes(intact_manifest)
    run(case, 'manifest-resumed.log')
    resumed_manifest_rows = manifest_path.read_text().splitlines()
    assert resumed_manifest_rows[:-1] == intact_manifest.decode().splitlines()
    assert resumed_manifest_rows[-1].split(',')[0] == '3'
    assert float(resumed_manifest_rows[-1].split(',')[1]) > \
        float(resumed_manifest_rows[-2].split(',')[1])
    passed('interrupted manifest: reject partial row unchanged, then resume after repair')

    source_case = root/'automatic-split'
    full = root/'automatic-full'
    for mode in ('both_missing', 'one_missing'):
        case = prepare(mode)
        shutil.copy2(source_case/'split_checkpoint.bin', case/'restart.bin')
        if mode == 'one_missing':
            diagnostics(case)
            original = require_output(full, 'TEST_2D_runout.csv').read_text().splitlines(keepends=True)
            prefix = [original[0]]+[line for line in original[1:] if float(line.split(',')[0]) <= 0.015]
            (diagnostics(case)/'TEST_2D_runout.csv').write_text(''.join(prefix))
        continuation(case)
        if mode == 'one_missing':
            allow_output_reuse(case)
        run(case, 'resumed.log')
        for suffix in ('runout', 'mass_center'):
            filename = 'TEST_2D_'+suffix+'.csv'
            original = require_output(full, filename).read_text().splitlines()
            expected = [original[0]]+[line for line in original[1:] if
                        (mode == 'one_missing' and suffix == 'runout') or float(line.split(',')[0]) > 0.015]
            assert require_output(case, filename).read_text().splitlines() == expected
        assert records(checkpoint(full)) == records(checkpoint(case))
        passed(mode+': headered segments with exact diagnostic suffix and final checkpoint')

    mutations = {
        'future': lambda s: s+'       4.000,       0.000,         0.000,         0.000,         0.000\n',
        'old': lambda s: ''.join(s.splitlines(keepends=True)[:-1]),
        'wrong_tail': lambda s: s.rsplit(',',1)[0]+',        42.000\n',
        'header': lambda s: s.replace('Time (s)', 'Wrong column'),
        'duplicate': lambda s: s+s.splitlines(keepends=True)[-1],
        'truncated': lambda s: s[:-12]+'\n',
        'nonfinite': lambda s: s.replace('       0.000,', '         NaN,', 1),
        'overflow': lambda s: s.replace('       0.000,', '       1e999,', 1),
        'empty': lambda s: '',
    }
    for name, mutate in mutations.items():
        case = prepare('reject-'+name)
        shutil.copy2(source_case/'split_checkpoint.bin', case/'restart.bin')
        diagnostics(case)
        for suffix in ('runout', 'mass_center'):
            filename = 'TEST_2D_'+suffix+'.csv'
            rows = require_output(full, filename).read_text().splitlines(keepends=True)
            prefix = [rows[0]]+[line for line in rows[1:] if float(line.split(',')[0]) <= 0.015]
            content = ''.join(prefix)
            if suffix == 'mass_center':
                content = mutate(content)
            (diagnostics(case)/filename).write_text(content)
        continuation(case)
        allow_output_reuse(case)
        before = {p.name: sha(p) for p in (case/'OUTPUT'/'TEST_2D'/'diagnostics').glob('*.csv')}
        checkpoint_before = sha(case/'restart.bin')
        run(case, 'rejected.log', 'Diagnostic CSV')
        assert before == {p.name:sha(p) for p in (case/'OUTPUT'/'TEST_2D'/'diagnostics').glob('*.csv')}
        assert checkpoint_before == sha(case/'restart.bin')
        passed('reject '+name+': both CSVs and checkpoint preserved')

    for name, key, value in (('cadence','DT_RUNOUT','0.25'),
                              ('threshold','EPS_STOP','1'), ('origin','X0_RUNOUT','6')):
        case = prepare('reject-'+name)
        shutil.copy2(source_case/'split_checkpoint.bin', case/'restart.bin')
        continuation(case)
        inp = case/'IMEX_SfloW2D.inp'
        content = inp.read_text()
        content = parameter(content, key, value) if name != 'origin' else content.replace(
            ' EPS_STOP=0,', ' EPS_STOP=0,\n X0_RUNOUT=6,\n Y0_RUNOUT=0,')
        inp.write_text(content)
        run(case, 'rejected.log', 'Restart runout cadence')
        assert not list((case/'OUTPUT'/'TEST_2D'/'diagnostics').glob('*.csv'))
        passed('reject incompatible '+name+' before opening CSVs')

    case = prepare('reject-v4')
    payloads = records(source_case/'split_checkpoint.bin')
    payloads[1] = struct.pack('<i',4)
    with (case/'restart.bin').open('wb') as stream:
        for payload in payloads[:-2]:
            marker = struct.pack('<i',len(payload))
            stream.write(marker+payload+marker)
    continuation(case); run(case, 'rejected.log', 'Version 5 is required')
    passed('reject version 4 with an explicit missing-history explanation')

    source = text.replace(' EPS_STOP=0,', ' EPS_STOP=1.0E99,')
    full = prepare('positive-stop-full', parameter(source, 'T_END','0.03'))
    split = prepare('positive-stop-split', source)
    run(full, 'full.log'); run(split,'first.log')
    assert struct.unpack('<4d',records(checkpoint(full))[4])[0] == 0.015
    continuation(split); run(split,'resumed.log')
    assert records(checkpoint(full)) == records(checkpoint(split))
    for suffix in ('runout','mass_center'):
        filename = 'TEST_2D_'+suffix+'.csv'
        assert require_output(full, filename).read_bytes() == \
            require_output(split, filename).read_bytes()
    passed('positive EPS_STOP: restored stop deadline prevents unintended additional evolution')
    manifest['status'] = 'pass'
    (root/'evidence.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Evidence:',root/'evidence.json',flush=True)


if __name__ == '__main__':
    main()
