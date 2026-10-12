"""Check N9 observation helpers without running or mutating the solver."""
import hashlib
import json
from pathlib import Path
import struct
import subprocess
import tempfile
import unittest

import numpy as np
import n9_long_cases as n9


REPO = Path(__file__).resolve().parents[2]


def archived_digest(revision, path):
    contents = subprocess.check_output(
        ["git", "-C", str(REPO), "show", f"{revision}:{path}"])
    return hashlib.sha256(contents).hexdigest()


class N9Helpers(unittest.TestCase):
    def test_dated_evidence_retains_physical_success_and_diagnostic_failure(self):
        """Keep the original v4 failure report as immutable historical evidence."""
        evidence = json.loads((REPO/'docs/acceptance/results_n9_long_cases_2026_10_11.json').read_text())
        self.assertFalse(evidence['scientific_baseline_approved'])
        self.assertFalse(evidence['production_sources_changed'])
        self.assertEqual(evidence['summary']['physical_continuations_passed'],8)
        self.assertEqual(evidence['summary']['diagnostic_continuations_failed'],8)
        self.assertEqual(evidence['summary']['failed'],8)
        self.assertEqual(len(evidence['thread_pairs']),4)
        self.assertTrue(all(v['status']=='pass' for v in evidence['thread_pairs']))
        self.assertRegex(evidence['runner_sha256'], r'^[0-9a-f]{64}$')
        self.assertTrue((REPO/evidence['report']).is_file())
        for name,value in evidence['retained_acceptance_contract_sha256'].items():
            self.assertEqual(n9.digest(REPO/name),value)
        self.assertTrue(evidence['clean_reference_provenance']['matches_frozen_repository_cores'])
        self.assertTrue(all(evidence['ETNA_fixture']['prior_user_raster_comparison'].values()))
        self.assertEqual({(r['case'],r['profile'],r['requested_threads']) for r in evidence['runs']},
                         {(c,p,t) for c in ('ETNA','ETNA_LAVA')
                          for p in ('historical_optimized','strict_debug') for t in (1,4)})
        for record in evidence['runs']:
            self.assertEqual(record['physical_continuation_status'],'pass')
            self.assertEqual(record['status'],'fail')
            self.assertTrue(record['ESRI_rasters_bitwise_equal'])
            self.assertTrue(record['checkpoint_bitwise_equal'])
            self.assertTrue(record['decoded_fields_bitwise_equal'])
            self.assertTrue(record['netcdf_schema_bitwise_equal'])
            self.assertTrue(all(value==0 for value in record['field_differences'].values()))
            for csv in record['diagnostic_csv'].values():
                self.assertEqual(csv['status'],'fail')
                self.assertFalse(csv['history_complete'])
            if record['case']=='ETNA_LAVA':
                budget=record['reference']['budget']
                self.assertEqual(budget['last_precontact_output_s'],3180)
                self.assertEqual(budget['first_boundary_positive_output_s'],3240)
                self.assertFalse(budget['exact_first_contact_time_known'])
                self.assertFalse(budget['mass_tolerance_changed'])
                self.assertFalse(budget['new_scientific_budget_threshold_approved'])

    def test_v5_restart_matrix_supersedes_historical_diagnostic_failures(self):
        evidence_path = REPO/'docs/acceptance/results_n9_restart_v5_2026_10_11.json'
        evidence = json.loads(evidence_path.read_text())
        self.assertEqual(evidence['summary'], {'passed': 8, 'failed': 0})
        self.assertTrue(evidence['scientific_baseline_approved'])
        sign_off = evidence['scientific_sign_off']
        self.assertEqual(sign_off['status'], 'approved')
        self.assertEqual(sign_off['approved_revision'], '3f7329c85d2e30e33cf68348721329f7726b01c1')
        self.assertIn('3180 s', sign_off['budget_policy']['ETNA_LAVA'])
        self.assertIn('not required', sign_off['raw_output_retention'])
        self.assertFalse(evidence['solver_tolerances_changed'])
        self.assertEqual(evidence['source_base_revision'], '0bab7d29daebe84c49d574ab2e6d6d0e293af5ac')
        for name, digest in evidence['source_sha256'].items():
            self.assertEqual(archived_digest(sign_off['approved_revision'], name), digest, name)
        self.assertRegex(evidence['runner_sha256'], r'^[0-9a-f]{64}$')
        self.assertEqual(len(evidence['runs']), 8)
        self.assertEqual(len(evidence['thread_pairs']), 4)
        self.assertTrue(all(pair['status'] == 'pass' for pair in evidence['thread_pairs']))
        for record in evidence['runs']:
            self.assertEqual(record['status'], 'pass')
            self.assertTrue(record['native_fields_restart_exact'])
            self.assertTrue(record['native_fields_prechange_exact'])
            self.assertTrue(record['checkpoint_v5_restart_exact'])
            self.assertTrue(record['canonical_records_prechange_exact'])
            self.assertTrue(record['CSV_prechange_full_exact'])
            self.assertTrue(record['ESRI_restart_and_prechange_exact'])
            for diagnostic in record['diagnostic_csv'].values():
                self.assertEqual(diagnostic['status'], 'pass')
                self.assertTrue(diagnostic['history_complete'])
                self.assertTrue(diagnostic['common_time_values_equal'])
        stochastic = evidence['stochastic_restart_matrix']
        self.assertTrue(stochastic['all_passed'])
        self.assertTrue(stochastic['actual_threads_checked_in_reference_first_and_resumed_runs'])
        self.assertRegex(stochastic['test_runner_sha256'], r'^[0-9a-f]{64}$')
        self.assertEqual({(run['profile'], run['actual_threads']) for run in stochastic['runs']},
                         {(profile, threads) for profile in ('historical_optimized', 'strict_debug')
                          for threads in (1, 4)})
        self.assertTrue(all(run['status'] == 'pass' and run['requested_threads'] == run['actual_threads']
                            for run in stochastic['runs']))
        for profile in ('historical_optimized', 'strict_debug'):
            pair = [run for run in stochastic['runs'] if run['profile'] == profile]
            self.assertEqual(len(pair), 2)
            self.assertEqual(pair[0]['checkpoint_sha256'], pair[1]['checkpoint_sha256'])
            self.assertEqual(pair[0]['snapshot_sha256'], pair[1]['snapshot_sha256'])
            self.assertTrue(all(len(run['harness_log_sha256']) == 64 for run in pair))
            long_run = next(run for run in evidence['runs'] if run['profile'] == profile)
            self.assertEqual(pair[0]['executable_sha256'], long_run['executable_sha256'])
        plan = json.loads((REPO/'TESTS/ACCEPTANCE/acceptance_plan.json').read_text())
        criteria = {criterion['id']: criterion for criterion in plan['criteria']}
        self.assertIn('results_n9_restart_v5_2026_10_11.json', criteria['N9-06']['gap'])
        self.assertEqual(criteria['N9-06']['coverage'], 'approved')
        self.assertIn('approved', criteria['N9-01']['coverage'])
        self.assertIn('runout origin/initial offsets', criteria['N9-07']['gap'])

    def test_parameter_override_is_unique_and_does_not_match_restart_files(self):
        text = '&RUN_PARAMETERS\n RESTART=F,\n T_END=3600,\n/\n&RESTART_PARAMETERS\n RESTART_FILES="a.asc",\n/\n'
        changed = n9.changed_parameter(text, 'RESTART', 'T')
        self.assertIn(' RESTART= T,', changed)
        self.assertIn(' RESTART_FILES="a.asc",', changed)
        with self.assertRaises(AssertionError):
            n9.changed_parameter(text+' T_END=100,\n', 'T_END', '200')
        with self.assertRaises(AssertionError):
            n9.changed_parameter(text, 'ABSENT_PARAMETER', '0')

    def records(self, path, *, version=4, bad_tail=False):
        q = np.zeros((4,2,3), dtype='<f8')
        q[0] = 2
        records = [b'IMEX_SFLOW2D_RST', struct.pack('<i',version),
                   struct.pack('<9i',2,3,4,0,0,0,0,6,5), struct.pack('<2i',0,0),
                   struct.pack('<4d',3600,.1,.1,.2), b'output schedule',
                   q.tobytes(order='F'), np.zeros((3,4),dtype='<f8').tobytes(order='F')]
        with path.open('wb') as stream:
            for i, record in enumerate(records):
                marker = struct.pack('<i',len(record))
                stream.write(marker+record+(b'bad!' if bad_tail and i==7 else marker))
        return q

    def test_checkpoint_header_layout_and_fortran_order(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)/'checkpoint.bin'
            expected = self.records(path)
            info, actual = n9.checkpoint(path)
            self.assertEqual(info['magic'], 'IMEX_SFLOW2D_RST')
            self.assertEqual(info['version'], 4)
            self.assertEqual(info['layout'][:3], [2,3,4])
            self.assertTrue(np.array_equal(actual,expected))
            self.assertEqual(info['conservative_sha256'],
                             hashlib.sha256(expected.tobytes(order='F')).hexdigest())

    def test_bad_checkpoint_version_and_record_markers_fail(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)/'checkpoint.bin'
            for options in ({'version':3},{'bad_tail':True}):
                self.records(path,**options)
                with self.assertRaises(AssertionError):
                    n9.checkpoint(path)

    def test_headerless_csv_is_read_but_lost_history_and_values_fail(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); full=root/'full';split=root/'split'
            full.mkdir();split.mkdir()
            for suffix in ('runout','mass_center'):
                (full/f'case_{suffix}.csv').write_text('Time (s), a,b,c,d\n0,0,1,2,3\n1,2,3,4,5\n')
                (split/f'case_{suffix}.csv').write_text('1,4000000,3,4,5\n')
            result = n9.diagnostic_comparison(full,split,'case')
            for item in result.values():
                self.assertEqual(item['status'],'fail')
                self.assertFalse(item['history_complete'])
                self.assertFalse(item['resumed_has_header'])
                self.assertEqual(item['common_time_Linf_by_column'][1],3999998)

    def test_exact_csv_and_nonfinite_or_duplicate_controls(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);full=root/'full';split=root/'split'
            full.mkdir();split.mkdir()
            for directory in (full,split):
                for suffix in ('runout','mass_center'):
                    (directory/f'case_{suffix}.csv').write_text('Time (s), a,b,c,d\n0,0,1,2,3\n1,2,3,4,5\n')
            result=n9.diagnostic_comparison(full,split,'case')
            self.assertTrue(all(v['status']=='pass' for v in result.values()))
            path=split/'bad.csv'
            for content in ('0,nan,1,2,3\n','0,0,1,2,3\n0,0,1,2,3\n'):
                path.write_text(content)
                with self.assertRaises(AssertionError):
                    n9.diagnostic_csv(path)


if __name__ == '__main__':
    unittest.main()
