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


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parameter(text, name, value):
    text, count = re.subn(r'(?mi)^(\s*'+name+r'\s*=)[^\n]*',
                         lambda m: m[1]+' '+value+',', text)
    assert count == 1, (name, count)
    return text


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('executable', type=Path)
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
        path = case/filename
        if path.is_file():
            return path
        details = [f'Missing expected output: {path}',
                   f'Files created: {sorted(p.name for p in case.iterdir())}']
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

    def passed(name):
        checks.append({'name': name, 'status':'pass'})
        print('PASS:', name, flush=True)

    for mode in ('automatic', 'explicit'):
        source = text if mode == 'automatic' else text.replace(' EPS_STOP=0,',
                                         ' EPS_STOP=0,\n X0_RUNOUT=5,\n Y0_RUNOUT=0,')
        full = prepare(mode+'-full', parameter(source, 'T_END', '0.03'))
        split = prepare(mode+'-split', source)
        run(full, 'full.log'); run(split, 'first.log')
        shutil.copy2(split/'restart.bin', split/'split_checkpoint.bin')
        continuation(split); run(split, 'resumed.log')
        assert records(full/'restart.bin') == records(split/'restart.bin')
        for suffix in ('runout', 'mass_center'):
            filename = 'TEST_2D_'+suffix+'.csv'
            assert require_output(full, filename).read_bytes() == \
                require_output(split, filename).read_bytes()
        assert (full/'TEST_2D_0002.q_2d').read_bytes() == (split/'TEST_2D_0002.q_2d').read_bytes()
        passed(mode+' origin: exact checkpoint, physical output and complete CSV history')

    source_case = root/'automatic-split'
    full = root/'automatic-full'
    for mode in ('both_missing', 'one_missing'):
        case = prepare(mode)
        shutil.copy2(source_case/'split_checkpoint.bin', case/'restart.bin')
        if mode == 'one_missing':
            original = (full/'TEST_2D_runout.csv').read_text().splitlines(keepends=True)
            prefix = [original[0]]+[line for line in original[1:] if float(line.split(',')[0]) <= 0.015]
            (case/'TEST_2D_runout.csv').write_text(''.join(prefix))
        continuation(case); run(case, 'resumed.log')
        for suffix in ('runout', 'mass_center'):
            filename = 'TEST_2D_'+suffix+'.csv'
            original = (full/filename).read_text().splitlines()
            expected = [original[0]]+[line for line in original[1:] if
                        (mode == 'one_missing' and suffix == 'runout') or float(line.split(',')[0]) > 0.015]
            assert (case/filename).read_text().splitlines() == expected
        assert records(full/'restart.bin') == records(case/'restart.bin')
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
        for suffix in ('runout', 'mass_center'):
            filename = 'TEST_2D_'+suffix+'.csv'
            rows = (full/filename).read_text().splitlines(keepends=True)
            prefix = [rows[0]]+[line for line in rows[1:] if float(line.split(',')[0]) <= 0.015]
            content = ''.join(prefix)
            if suffix == 'mass_center':
                content = mutate(content)
            (case/filename).write_text(content)
        continuation(case)
        before = {p.name: sha(p) for p in case.glob('*.csv')}
        checkpoint_before = sha(case/'restart.bin')
        run(case, 'rejected.log', 'Diagnostic CSV')
        assert before == {p.name:sha(p) for p in case.glob('*.csv')}
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
        assert not list(case.glob('*.csv'))
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
    assert struct.unpack('<4d',records(full/'restart.bin')[4])[0] == 0.015
    continuation(split); run(split,'resumed.log')
    assert records(full/'restart.bin') == records(split/'restart.bin')
    for suffix in ('runout','mass_center'):
        filename = 'TEST_2D_'+suffix+'.csv'
        assert (full/filename).read_bytes() == (split/filename).read_bytes()
    passed('positive EPS_STOP: restored stop deadline prevents unintended additional evolution')
    manifest['status'] = 'pass'
    (root/'evidence.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Evidence:',root/'evidence.json',flush=True)


if __name__ == '__main__':
    main()
