"""Non-mutating N9 long-case and continuation characterization.

All runs live under a fresh temporary root. The production source, canonical
example inputs and numerical tolerances are untouched. Native NetCDF is read
through its C API, sequentially, without a Python netCDF4 dependency.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import ctypes as C
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import sys
import tempfile
import time
from zipfile import ZipFile
import numpy as np


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(chunk)
    return result.hexdigest()


def changed_parameter(text, name, value):
    result, count = re.subn(r'(?mi)^(\s*' + name + r'\s*=)[^\n]*',
                            lambda m: m[1] + ' ' + value + ',', text)
    assert count == 1, (name, count)
    return result


class NetCDF:
    """Read numeric arrays and descriptive attributes using installed libnetcdf."""
    def __init__(self):
        libdir = Path(subprocess.check_output(['nc-config', '--libdir'], text=True).strip())
        library = next(p for p in (libdir/'libnetcdf.dylib', libdir/'libnetcdf.so') if p.exists())
        self.library = str(library)
        self.nc = C.CDLL(str(library))
        pointer = C.POINTER(C.c_int)
        for name, args in {
            'nc_open': [C.c_char_p, C.c_int, pointer], 'nc_close': [C.c_int],
            'nc_inq': [C.c_int, pointer, pointer, pointer, pointer],
            'nc_inq_var': [C.c_int, C.c_int, C.c_char_p, pointer, pointer, pointer, pointer],
            'nc_inq_dim': [C.c_int, C.c_int, C.c_char_p, C.POINTER(C.c_size_t)],
            'nc_get_var_double': [C.c_int, C.c_int, C.POINTER(C.c_double)],
            'nc_inq_attlen': [C.c_int, C.c_int, C.c_char_p, C.POINTER(C.c_size_t)],
            'nc_get_att_text': [C.c_int, C.c_int, C.c_char_p, C.c_char_p]
        }.items():
            getattr(self.nc, name).argtypes = args
            getattr(self.nc, name).restype = C.c_int
        self.nc.nc_strerror.argtypes = [C.c_int]
        self.nc.nc_strerror.restype = C.c_char_p

    def check(self, status):
        if status:
            raise AssertionError(self.nc.nc_strerror(status).decode())

    def read(self, path):
        handle = C.c_int()
        self.check(self.nc.nc_open(os.fsencode(path), 0, C.byref(handle)))
        fields, metadata = {}, {}
        try:
            ndim, nvar, natt, unlim = (C.c_int() for _ in range(4))
            self.check(self.nc.nc_inq(handle, C.byref(ndim), C.byref(nvar),
                                      C.byref(natt), C.byref(unlim)))
            for var in range(nvar.value):
                name = C.create_string_buffer(257)
                xtype, rank, attributes = (C.c_int() for _ in range(3))
                ids = (C.c_int * 1024)()
                self.check(self.nc.nc_inq_var(handle, var, name, C.byref(xtype),
                                              C.byref(rank), ids, C.byref(attributes)))
                shape, dimensions = [], []
                for index in ids[:rank.value]:
                    label = C.create_string_buffer(257)
                    length = C.c_size_t()
                    self.check(self.nc.nc_inq_dim(handle, index, label, C.byref(length)))
                    shape.append(length.value)
                    dimensions.append(label.value.decode())
                data = np.empty(shape, dtype=np.float64)
                self.check(self.nc.nc_get_var_double(handle, var,
                    data.ctypes.data_as(C.POINTER(C.c_double))))
                description = {}
                for key in ('units', 'long_name'):
                    length = C.c_size_t()
                    self.check(self.nc.nc_inq_attlen(handle, var, key.encode(), C.byref(length)))
                    value = C.create_string_buffer(length.value + 1)
                    self.check(self.nc.nc_get_att_text(handle, var, key.encode(), value))
                    description[key] = value.value.decode()
                key = name.value.decode()
                metadata[key] = dict(shape=shape, dimensions=dimensions,
                                     nc_type=xtype.value, **description)
                fields[key] = data
        finally:
            self.check(self.nc.nc_close(handle))
        return fields, metadata


def checkpoint(path):
    """Decode GNU Fortran little-endian sequential records, explicitly scoped."""
    records = []
    with path.open('rb') as stream:
        while marker := stream.read(4):
            assert len(marker) == 4
            length = struct.unpack('<i', marker)[0]
            assert length >= 0, 'split Fortran records unsupported by this comparator'
            record = stream.read(length)
            assert len(record) == length and stream.read(4) == marker
            records.append(record)
    assert records[0] == b'IMEX_SFLOW2D_RST'
    assert struct.unpack('<i', records[1])[0] == 4
    dims = struct.unpack('<9i', records[2])
    nx, ny, nv, ns, ng, nz, npore, nth, npress = dims
    state = np.frombuffer(records[6], '<f8').reshape((nv, nx, ny), order='F')
    bed = np.frombuffer(records[7], '<f8').reshape((nx+1, ny+1), order='F')
    assert np.isfinite(state).all() and np.isfinite(bed).all()
    assert state[0].min() >= 0
    return {'magic': records[0].decode().rstrip(), 'version': 4,
            'layout': list(dims), 'model_flags': list(struct.unpack('<2i', records[3])),
            'time_state': list(struct.unpack('<4d', records[4])),
            'records': len(records), 'sha256': digest(path),
            'conservative_sha256': hashlib.sha256(records[6]).hexdigest(),
            'nodal_bed_sha256': hashlib.sha256(records[7]).hexdigest()}, state


def diagnostic_csv(path):
    """Read either the normal CSV header or the headerless resumed output."""
    with path.open() as stream:
        first = stream.readline()
    header = 'Time' in first
    values = np.loadtxt(path, delimiter=',', skiprows=int(header), ndmin=2)
    assert values.shape[1] == 5 and np.isfinite(values).all()
    assert np.all(np.diff(values[:,0]) > 0), 'duplicate or unordered diagnostic times'
    return values, header


def diagnostic_comparison(full, resumed, run_name):
    """Retain lost history and wrong common-time diagnostics as failed checks."""
    checks = {}
    for suffix in ('runout', 'mass_center'):
        filename = run_name+'_'+suffix+'.csv'
        left, left_header = diagnostic_csv(full/filename)
        right, right_header = diagnostic_csv(resumed/filename)
        common, il, ir = np.intersect1d(left[:,0],right[:,0],return_indices=True)
        errors = np.max(np.abs(left[il]-right[ir]),axis=0) if len(common) else np.full(5,np.inf)
        same_history = left.shape == right.shape and np.array_equal(left[:,0],right[:,0])
        same_values = bool(len(common) and np.array_equal(left[il],right[ir]))
        checks[suffix] = {'full_rows':len(left),'resumed_rows':len(right),
            'full_has_header':left_header,'resumed_has_header':right_header,
            'history_complete':same_history,'common_time_values_equal':same_values,
            'common_times':len(common), 'common_time_Linf_by_column':errors.tolist() if len(common) else None,
            'full_sha256':digest(full/filename),'resumed_sha256':digest(resumed/filename),
            'status':'pass' if same_history and same_values and left_header==right_header else 'fail'}
    return checks


def observations(fields, metadata, restart, case):
    """Check native physical schema and observe budgets without new tolerances."""
    times, h = fields['time'], fields['h']
    end, interval = (200, 5) if case == 'ETNA' else (3600, 60)
    assert np.array_equal(times, np.arange(0, end+1, interval))
    assert np.array_equal(fields['x'], 499813+(np.arange(250)+.5)*10)
    assert np.array_equal(fields['y'], 4176700+(np.arange(150)+.5)*10)
    for name, values in fields.items():
        assert np.isfinite(values).all() and not np.any(values == 9.969209968386869e36), name
        assert metadata[name]['nc_type'] == 6, (name, 'not NC_DOUBLE')
        expected = [name] if name in ('x', 'y', 'time') else ['time', 'y', 'x']
        assert metadata[name]['dimensions'] == expected, name
        assert metadata[name]['long_name'], name
    for name, unit in [('h', 'meters'), ('b', 'meters'), ('eta', 'meters'),
                       ('T', 'kelvin'), ('u', 'meters/seconds'), ('v', 'meters/seconds'),
                       ('w', 'meters/seconds'), ('rhomix', 'kg/m3')]:
        assert metadata[name]['units'] == unit, (name, metadata[name])
    assert h.min() >= 0
    assert np.array_equal(fields['b'], np.broadcast_to(fields['b'][0], h.shape))
    eta_error = float(np.max(np.abs(fields['eta']-(fields['b']+h))))
    # Existing field-roundoff multiplier, not a changed solver/test limit.
    assert eta_error <= 8192*np.finfo(float).eps*max(1., float(np.max(np.abs(fields['eta']))))
    volume = 100*np.sum(h, axis=(1,2))
    # The native output intentionally zeroes density/composition below dry
    # depth, not conservative mass. Pure liquid has constant density, so h
    # alone retains its entire inventory, including positive dry residues.
    masked_mass = 100*np.sum(h*fields['rhomix'], axis=(1,2))
    mass = 2600*volume if case == 'ETNA_LAVA' else masked_mass
    edge = np.maximum(np.max(h[:,(0,-1),:], axis=(1,2)),
                      np.max(h[:,:,(0,-1)], axis=(1,2)))
    contacts = np.flatnonzero(edge > 0)
    stop = int(contacts[0]) if contacts.size else len(times)
    assert stop > 0
    budget = {'scope': 'saved snapshots before first strictly positive boundary depth',
              'last_precontact_output_s': float(times[stop-1]),
              'first_boundary_positive_output_s': float(times[stop]) if contacts.size else None,
              'exact_first_contact_time_known': False,
              'boundary_thickness_max_m': float(edge.max()),
              'time_s': times[:stop].tolist(), 'volume_m3': volume[:stop].tolist(),
              'mass_kg': mass[:stop].tolist(),
              'mass_tolerance_changed': False,
              'new_scientific_budget_threshold_approved': False}
    if case == 'ETNA_LAVA':
        injected = 9.24*times[:stop]
        vres = volume[:stop]-volume[0]-injected
        mres = mass[:stop]-mass[0]-2600*injected
        budget.update(volume_rate_m3_s=9.24, density_kg_m3=2600,
                      expected_volume_m3=injected.tolist(),
                      volume_residual_m3=vres.tolist(), mass_residual_kg=mres.tolist(),
                      max_abs_volume_residual_m3=float(np.max(np.abs(vres))),
                      max_abs_mass_residual_kg=float(np.max(np.abs(mres))),
                      max_scaled_volume_residual=float(np.max(np.abs(vres)/np.maximum(1,injected))),
                      max_scaled_mass_residual=float(np.max(np.abs(mres)/np.maximum(1,2600*injected))),
                      user_authorized_scope='2026-10-11: mass check only before flow reaches boundary cells')
        assert times[stop-1] == 3180 and times[stop] == 3240
    else:
        residual = mass[:stop]-mass[0]
        solid = 100*np.sum(h*fields['solid_frac_01']*2000, axis=(1,2))
        budget.update(mass_residual_kg=residual.tolist(),
                      max_abs_mass_residual_kg=float(np.max(np.abs(residual))),
                      max_scaled_mass_residual=float(np.max(np.abs(residual))/mass[0]),
                      max_abs_solid_mass_residual_kg=float(np.max(np.abs(solid[:stop]-solid[0]))),
                      volume_conservation_required=False,
                      volume_note='Gas density depends on T; mixture volume is not the conservative mass.',
                      intermediate_mass_note='Native rhomix and solid fractions exclude dry-depth residuals; these intermediate sums are observations, not the conservative budget gate.')
    info, q = checkpoint(restart)
    expected_nv, expected_ns = (5, 1) if case == 'ETNA' else (4, 0)
    assert info['layout'][:7] == [250, 150, expected_nv, expected_ns, 0, 0, 0]
    assert info['model_flags'] == [0,0] and info['time_state'][0] == end
    final_mass = 100*np.sum(q[0])
    wet = h[-1] > 1.e-10
    dry = ~wet
    assert np.all(fields['rhomix'][-1][dry] == 0)
    if case == 'ETNA':
        assert np.all(fields['solid_frac_01'][-1][dry] == 0)
        # Initial ETNA data are either exactly empty or safely above the dry
        # cutoff, so its initial native sum contains the whole initial mass.
        assert np.min(h[0][h[0] > 0]) > 1.e-10
        initial_solid = 100*np.sum(h[0]*fields['solid_frac_01'][0]*2000)
        budget.update(conservative_endpoint_initial_mass_kg=float(mass[0]),
                      conservative_endpoint_final_mass_kg=float(final_mass),
                      conservative_endpoint_mass_residual_kg=float(final_mass-mass[0]),
                      conservative_endpoint_solid_mass_residual_kg=float(100*np.sum(q[4])-initial_solid),
                      final_retained_dry_mass_kg=float(100*np.sum(q[0].T[dry])),
                      conservative_endpoint_scope='t=0 native fields (no positive dry cells) versus full raw conservative checkpoint at t=200; no boundary contact in all saved outputs')
    decoded_mass_error = float(np.max(np.abs(q[0].T-h[-1]*fields['rhomix'][-1])[wet]))
    thermal_capacity = q[0]*1150 if case == 'ETNA_LAVA' else q[4]*1617+(q[0]-q[4])*998
    thermal_error = float(np.max(np.abs(q[3].T-thermal_capacity.T*fields['T'][-1])))
    scale_mass = max(1., float(np.max(np.abs(q[0]))))
    scale_thermal = max(1., float(np.max(np.abs(q[3]))))
    assert decoded_mass_error <= 8192*np.finfo(float).eps*scale_mass
    assert thermal_error <= 8192*np.finfo(float).eps*scale_thermal
    if case == 'ETNA':
        component_error = float(np.max(np.abs(q[4].T-h[-1]*fields['solid_frac_01'][-1]*2000)[wet]))
        assert component_error <= 8192*np.finfo(float).eps*max(1,float(q[4].max()))
    else:
        component_error = None
    return {'snapshots': len(times), 'variables': len(fields), 'metadata': metadata,
            'decoded_sha256': {k: hashlib.sha256(v.tobytes()).hexdigest() for k,v in fields.items()},
            'eta_identity_Linf_m': eta_error, 'checkpoint': info, 'budget': budget,
            'final_volume_m3': float(volume[-1]), 'final_mass_kg': float(final_mass),
            'checkpoint_physical_mass_Linf_kg_m2': decoded_mass_error,
            'checkpoint_physical_mass_scope': 'wet cells only; native dry density and fractions are zero by existing schema',
            'checkpoint_thermal_energy_Linf_J_m2': thermal_error,
            'checkpoint_solid_mass_Linf_kg_m2': component_error}


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--repo', type=Path, required=True)
    cli.add_argument('--build-audit', type=Path, required=True)
    cli.add_argument('--lava-evidence', type=Path, required=True)
    cli.add_argument('--preprocess-adapter', type=Path,
                     default=Path(__file__).with_name('n9_etna_preprocess.py'))
    cli.add_argument('--reuse-first-halves', type=Path,
                     help='recover completed original-input full/first-half runs; verify their inputs/time/schema before continuation')
    cli.add_argument('--characterize-completed', type=Path,
                     help='only inspect completed full/split runs in this root; never rerun or overwrite them')
    cli.add_argument('--output', type=Path, required=True)
    cli.add_argument('--reference-zip', type=Path,
                     help='verify the frozen dry-safe clean cores; not the earlier September 24 package')
    args = cli.parse_args()
    if args.reuse_first_halves and args.characterize_completed:
        cli.error('choose one recovery mode only')
    repo, root = args.repo.resolve(), args.output.resolve()
    assert not root.exists() and repo not in root.parents
    root.mkdir()
    audit = json.loads((args.build_audit/'manifest.json').read_text())
    lava = json.loads(args.lava_evidence.read_text())
    for name, value in lava['source_sha256'].items():
        assert digest(repo/name) == value, ('changed production source', name)
        assert audit['source_and_fixture_sha256'][name] == value, ('build provenance mismatch', name)
    manifest = {'schema_version': 1, 'started_utc': datetime.now(timezone.utc).isoformat(),
                'source_revision': subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
                'source_sha256': lava['source_sha256'], 'production_sources_changed': False,
                'tools': audit['tools'], 'platform': audit['platform'],
                'build_audit_sha256': digest(args.build_audit/'manifest.json'),
                'prior_lava_evidence_sha256': digest(args.lava_evidence),
                'runner_sha256': digest(Path(__file__)),
                'runs': [], 'continuations': [], 'thread_pairs': [],
                'N9_scientific_approval': False,
                'comparison_contract': {'same_profile_thread_fields': 'bitwise',
                    'same_profile_continuation_fields': 'bitwise',
                    'checkpoint_records': 'bitwise', 'cross_profile': 'observations only',
                    'budget': 'raw absolute/scaled residuals; no new scientific threshold',
                    'lava_budget_end_s': 3180, 'roundoff_identity_multiplier': 8192}}
    if args.reference_zip:
        expected = json.loads((repo/'TESTS/TEST_PCCU_GATE_H/field_contract_v3.json').read_text())['reference_core_sha256']
        with ZipFile(args.reference_zip) as archive:
            members = {Path(name).name:hashlib.sha256(archive.read(name)).hexdigest()
                       for name in archive.namelist() if Path(name).name in expected}
        assert members == expected
        assert all(digest(repo/'TESTS/TEST_SPATIAL_OPERATOR/reference'/name)==value
                   for name,value in expected.items())
        manifest['clean_reference_provenance'] = {'archive':args.reference_zip.name,
            'archive_sha256':digest(args.reference_zip),'core_sha256':members,
            'matches_frozen_repository_cores':True,
            'external_full_suite_rerun':False}
    def save():
        (root/'evidence.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
    save()
    recovery = args.reuse_first_halves or args.characterize_completed
    case_root = recovery.resolve() if recovery else root
    data = case_root/'ETNA-data'
    example = repo/'EXAMPLES/EXAMPLE_ETNA'
    archive = example/'DEM/Etna2014_crop.zip'
    if not recovery:
        data.mkdir(); (data/'DEM').mkdir()
        for name in ('IMEX_SfloW2D.inp', 'create_input_ellipsoid.py','input_ellipsoid.py'):
            shutil.copy2(example/name, data/name)
        with ZipFile(archive) as z:
            with z.open('Etna2014_crop.asc') as src, (data/'DEM/Etna2014_crop.asc').open('xb') as dst:
                shutil.copyfileobj(src,dst)
        with (data/'preprocessing.log').open('w') as log:
            subprocess.run([sys.executable,str(args.preprocess_adapter.resolve())],cwd=data,
                           stdout=log,stderr=subprocess.STDOUT,check=True,timeout=300)
    else:
        # Recovery reuses output payloads, not an unverified scientific baseline.
        previous = json.loads((case_root/'evidence.json').read_text())
        assert previous['source_sha256'] == manifest['source_sha256']
        assert previous['build_audit_sha256'] == manifest['build_audit_sha256']
        for name in ('IMEX_SfloW2D.inp','create_input_ellipsoid.py','input_ellipsoid.py'):
            assert digest(data/name) == digest(example/name)
        for name, value in previous['ETNA_fixture']['files'].items():
            assert digest(data/name) == value
        manifest['recovery'] = {'first_half_root': str(case_root),
                                'previous_manifest_sha256': digest(case_root/'evidence.json'),
            'prior_status': 'incomplete harness; no scientific pass claimed',
            'mode': 'read-only completed characterization' if args.characterize_completed else 'continue first halves'}
    manifest['ETNA_fixture'] = {'input_sha256': digest(data/'IMEX_SfloW2D.inp'),
        'archive_sha256': digest(archive),
        'preprocessing': 'tracked geometry script unchanged; NumPy numeric-table parser replaces pandas reader only',
        'adapter_sha256': digest(args.preprocess_adapter),
        'files': {str(p.relative_to(data)): digest(p) for p in data.rglob('*') if p.is_file()}}
    old_data = Path('/Users/demichie/Codes/IMEX_SfloW2D_dev/EXAMPLES/EXAMPLE_ETNA/DEM')
    manifest['ETNA_fixture']['prior_user_raster_comparison'] = {
        name: digest(data/'DEM'/name) == digest(old_data/name)
        for name in ('new_dem.asc','init_volume.asc','ellipsoid.asc') if (old_data/name).exists()}
    save(); print('PREPROCESSING',manifest['ETNA_fixture']['prior_user_raster_comparison'],flush=True)
    spec = importlib.util.spec_from_file_location('lava_prepare',repo/'EXAMPLES/EXAMPLE_ETNA_LAVA/prepare.py')
    prep = importlib.util.module_from_spec(spec); spec.loader.exec_module(prep)
    profiles = ('historical_optimized','strict_debug')
    teams = (1,4)
    native = NetCDF()

    def prepare(case, directory):
        if case == 'ETNA':
            directory.mkdir()
            shutil.copy2(data/'IMEX_SfloW2D.inp',directory/'IMEX_SfloW2D.inp')
            (directory/'DEM').symlink_to(data/'DEM',target_is_directory=True)
        else:
            prep.prepare(directory)

    def run(exe, directory, threads, logname):
        start = time.monotonic()
        with (directory/logname).open('w') as log:
            proc = subprocess.run([str(exe)],cwd=directory,
                env={**os.environ,'OMP_NUM_THREADS':str(threads),'OMP_DYNAMIC':'FALSE'},
                stdout=log,stderr=subprocess.STDOUT,timeout=1800)
        text = (directory/logname).read_text()
        team = re.search(r'Parallel run: number of threads used\s+(\d+)',text)
        assert proc.returncode == 0 and team and int(team[1]) == threads, (directory,logname)
        return {'returncode': proc.returncode,'seconds': time.monotonic()-start,
                'actual_threads': int(team[1]),'log_sha256': digest(directory/logname),
                'input_sha256': digest(directory/'IMEX_SfloW2D.inp')}

    def recovered_run(directory, threads, logname, input_text, expected_time):
        log = directory/logname
        # An unfinished solver may still own these files; wait for its normal
        # final log marker, never read a checkpoint while it is being written.
        deadline = time.monotonic()+1800
        while not log.exists() or 'Total elapsed real time is' not in log.read_text():
            if time.monotonic() >= deadline:
                raise TimeoutError(f'Incomplete recovered run: {log}')
            time.sleep(1)
        text = log.read_text()
        team = re.search(r'Parallel run: number of threads used\s+(\d+)',text)
        assert team and int(team[1]) == threads
        assert (directory/'IMEX_SfloW2D.inp').read_text() == input_text
        info, _ = checkpoint(directory/'restart.bin')
        assert info['time_state'][0] == expected_time
        return {'recovered_completed_run': True,'actual_threads':threads,
                'seconds':None,'log_sha256':digest(log),
                'input_sha256':digest(directory/'IMEX_SfloW2D.inp')}

    def job(case, profile, threads):
        exe = args.build_audit/profile/'src/IMEX_SfloW2D'
        assert digest(exe) == audit['profiles'][profile]['executable_sha256']
        tag = f'{case}-{profile}-{threads}'
        reference = case_root/(tag+'-full')
        if case == 'ETNA':
            if recovery:
                full_run = recovered_run(reference,threads,'run.log',
                                         (data/'IMEX_SfloW2D.inp').read_text(),200)
            else:
                prepare(case,reference)
                full_run = run(exe,reference,threads,'run.log')
        else:
            prior = next(r for r in lava['results'] if r['profile']==profile and r['actual_threads']==threads)
            reference = Path(prior['directory'])
            assert digest(reference/'restart.bin') == prior['restart_sha256']
            assert digest(reference/'exampleEtna_lava.nc') == prior['NetCDF_sha256']
            full_run = {'reused_verified_full_run': True, 'actual_threads':threads,
                        'seconds':prior['seconds'],'input_sha256':prior['input_sha256']}
        resumed = case_root/(tag+'-split')
        if not recovery:
            prepare(case,resumed)
        input_path = resumed/'IMEX_SfloW2D.inp'
        canonical = (repo/('EXAMPLES/EXAMPLE_'+case)/'IMEX_SfloW2D.inp').read_text()
        split_time = 100 if case == 'ETNA' else 1800
        first_input = changed_parameter(canonical,'T_END',str(split_time))
        second_input = changed_parameter(canonical,'RESTART','T')
        if case == 'ETNA':
            second_input = changed_parameter(second_input,'RESTART_FILES','"restart.bin"')
        else:
            second_input += '\n&RESTART_PARAMETERS\n RESTART_FILES="restart.bin",\n/\n'
        if args.characterize_completed:
            second_run = recovered_run(resumed,threads,'resumed_run.log',second_input,
                                       200 if case=='ETNA' else 3600)
            first_log = (resumed/'first_run.log').read_text()
            assert 'Total elapsed real time is' in first_log
            team = re.search(r'Parallel run: number of threads used\s+(\d+)',first_log)
            assert team and int(team[1]) == threads
            first_run = {'recovered_completed_run':True,'actual_threads':threads,
                         'log_sha256':digest(resumed/'first_run.log'),
                         'expected_input_sha256':hashlib.sha256(first_input.encode()).hexdigest(),
                         'input_identity_scope':'first input generated by recorded harness; current file is the independently verified continuation input'}
            first_checkpoint,_ = checkpoint(resumed/'split_checkpoint.bin')
            assert first_checkpoint['time_state'][0] == split_time
            return {'case':case,'profile':profile,'requested_threads':threads,
                'flags':audit['profiles'][profile]['flags'],'executable_sha256':digest(exe),
                'reference_directory':str(reference),'resumed_directory':str(resumed),
                'reference_run':full_run,'first_run':first_run,'second_run':second_run,
                'split_time_s':split_time,'split_checkpoint':first_checkpoint,
                'input_override_policy':'temporary first T_END; second RESTART and RESTART_FILES only'}
        if args.reuse_first_halves:
            first_run = recovered_run(resumed,threads,'first_run.log',first_input,split_time)
        else:
            input_path.write_text(first_input)
            first_run = run(exe,resumed,threads,'first_run.log')
        first_checkpoint,_ = checkpoint(resumed/'restart.bin')
        assert first_checkpoint['time_state'][0] == split_time
        # Preserve the intermediate evidence before production overwrites it.
        shutil.copy2(resumed/'restart.bin',resumed/'split_checkpoint.bin')
        input_path.write_text(second_input)
        second_run = run(exe,resumed,threads,'resumed_run.log')
        print('RUNS FINISHED',tag,flush=True)
        return {'case':case,'profile':profile,'requested_threads':threads,
                'flags':audit['profiles'][profile]['flags'],'executable_sha256':digest(exe),
                'reference_directory':str(reference),'resumed_directory':str(resumed),
                'reference_run':full_run,'first_run':first_run,'second_run':second_run,
                'split_time_s':split_time,'split_checkpoint':first_checkpoint,
                'input_override_policy':'temporary first T_END; second RESTART and RESTART_FILES only'}

    jobs = [(case,p,t) for case in ('ETNA','ETNA_LAVA') for p in profiles for t in teams]
    # Solver processes may run concurrently. HDF5/NetCDF reads stay serial here.
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = {pool.submit(job,*j):j for j in jobs}
        for future in as_completed(futures):
            try:
                result = future.result()
            except Exception as error:
                case, profile, threads = futures[future]
                manifest['runs'].append({'case':case,'profile':profile,
                    'requested_threads':threads,'status':'fail','error':repr(error)})
                save(); print('FAILED RUN',futures[future],repr(error),flush=True)
                continue
            case = result['case']; name = 'exampleEtna' if case == 'ETNA' else 'exampleEtna_lava'
            directory = Path(result['reference_directory'])
            resumed = Path(result['resumed_directory'])
            try:
                left, lm = native.read(directory/(name+'.nc'))
                right, rm = native.read(resumed/(name+'.nc'))
                result['reference'] = observations(left,lm,directory/'restart.bin',case)
                result['resumed'] = observations(right,rm,resumed/'restart.bin',case)
            except Exception as error:
                result.update(status='fail',error=repr(error))
                manifest['runs'].append(result);save()
                print('FAILED OBSERVATION',case,result['profile'],repr(error),flush=True)
                continue
            result['netcdf_schema_bitwise_equal'] = lm == rm
            result['decoded_fields_bitwise_equal'] = result['reference']['decoded_sha256'] == result['resumed']['decoded_sha256']
            result['checkpoint_bitwise_equal'] = digest(directory/'restart.bin') == digest(resumed/'restart.bin')
            rasters = {p.name:digest(p) for p in directory.glob('*.asc')}
            resumed_rasters = {p.name:digest(p) for p in resumed.glob('*.asc')}
            result['ESRI_rasters_bitwise_equal'] = rasters == resumed_rasters
            result['ESRI_raster_files_compared'] = len(rasters)
            result['physical_continuation_status'] = 'pass' if all(result[k] for k in (
                'netcdf_schema_bitwise_equal','decoded_fields_bitwise_equal','checkpoint_bitwise_equal',
                'ESRI_rasters_bitwise_equal')) else 'fail'
            result['diagnostic_csv'] = diagnostic_comparison(directory,resumed,name)
            result['status'] = 'pass' if result['physical_continuation_status']=='pass' and all(
                v['status']=='pass' for v in result['diagnostic_csv'].values()) else 'fail'
            result['field_differences'] = {k:float(np.max(np.abs(left[k]-right[k]))) for k in left}
            result['full_NetCDF_sha256'] = digest(directory/(name+'.nc'))
            result['resumed_NetCDF_sha256'] = digest(resumed/(name+'.nc'))
            result['reference_output_sha256'] = {p.name:digest(p) for p in directory.iterdir()
                if p.is_file() and (p.suffix in ('.asc','.csv','.nc') or p.name=='restart.bin')}
            manifest['runs'].append(result); save()
            print('CONTINUATION',case,result['profile'],result['requested_threads'],result['status'],flush=True)
    for case in ('ETNA','ETNA_LAVA'):
        for profile in profiles:
            one = next(r for r in manifest['runs'] if (r['case'],r['profile'],r['requested_threads'])==(case,profile,1))
            four = next(r for r in manifest['runs'] if (r['case'],r['profile'],r['requested_threads'])==(case,profile,4))
            match = (one.get('physical_continuation_status')==four.get('physical_continuation_status')=='pass') and all(one[mode]['decoded_sha256']==four[mode]['decoded_sha256'] and
                        one[mode]['checkpoint']['sha256']==four[mode]['checkpoint']['sha256']
                        for mode in ('reference','resumed'))
            manifest['thread_pairs'].append({'case':case,'profile':profile,'status':'pass' if match else 'fail'})
    manifest['summary'] = {'continuations':len(manifest['runs']),
        'failed':sum(r['status']!='pass' for r in manifest['runs']),
        'physical_continuations_passed':sum(r.get('physical_continuation_status')=='pass' for r in manifest['runs']),
        'diagnostic_continuations_failed':sum(any(v['status']!='pass' for v in r.get('diagnostic_csv',{}).values()) for r in manifest['runs']),
        'thread_pairs':len(manifest['thread_pairs']),
        'N9_approved':False}
    manifest['finished_utc'] = datetime.now(timezone.utc).isoformat()
    save(); print('SUMMARY',manifest['summary'],flush=True)
    if manifest['summary']['failed'] or any(r['status']!='pass' for r in manifest['thread_pairs']):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
