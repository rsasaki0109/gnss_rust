"""Verify retained outputs, hashes, status counts and reference comparisons."""
import argparse
import collections
import csv
from decimal import Decimal
import gzip
import hashlib
import io
import json
from pathlib import Path
import runpy
import subprocess


def digest(data):
    return hashlib.sha256(data).hexdigest()


def rows(data):
    return list(csv.DictReader(io.StringIO(data.decode())))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-commit')
    a = p.parse_args()
    root = Path(__file__).resolve().parent
    repo = root.parents[2]
    manifest = json.loads((root/'predeclared-inputs.json').read_text())
    evidence = json.loads((root/'evidence-sha256.json').read_text())
    for f, h in evidence.items():
        assert digest((root/f).read_bytes()) == h, f
    def source(f):
        if a.source_commit:
            return subprocess.check_output(['git','show',a.source_commit+':'+f],cwd=repo)
        return (repo/f).read_bytes()
    for f, h in manifest['source_sha256'].items():
        assert digest(source(f)) == h, f
    native = json.loads((root/'native-generator.json').read_text())
    assert native['repeat_runs'][0]['fixture_sha256'] == native['repeat_runs'][1]['fixture_sha256']
    for f, h in native['repeat_runs'][0]['fixture_sha256'].items():
        assert digest(source('tests/fixtures/'+f)) == h
    for f, h in native['source_sha256'].items():
        assert digest(source(f)) == h, f
    assert native['actual_formMeasurementEquations_calls_per_repeat'] == 64
    assert native['actual_updateFilter_calls_per_repeat'] == 0
    runner = {}
    exec(compile(source('scripts/run_ppp_reference_update.py'),'retained-runner','exec'),runner)
    normalize = runner['normalized']
    scoring = {}
    exec(compile(source('scripts/run_ppc_rtk_float.py'),'retained-scorer','exec'),scoring)
    comparison = {}
    exec(compile(source('scripts/run_ppc_ppp_initial_prior.py'),'retained-comparison','exec'),comparison)
    scorer, compare = scoring['score'], comparison['compare']
    old = repo/manifest['frozen_control_archive']
    frozen = repo/manifest['frozen_native_archive']
    assert digest((frozen/'evidence-sha256.json').read_bytes()) == manifest['frozen_native_evidence_sha256']
    old_manifest = json.loads((old/'predeclared-inputs.json').read_text())
    report = json.loads((root/'report.json').read_text())
    assert len(report['profiles']) == 8
    total = 0
    for name, cmd in manifest['commands'].items():
        execution = json.loads((root/(name+'.execution.json')).read_text())
        data = gzip.decompress((root/(name+'.stdout.csv.gz')).read_bytes())
        stderr = gzip.decompress((root/(name+'.stderr.log.gz')).read_bytes())
        assert execution['command'] == cmd and execution['returncode'] == 0
        assert digest(data) == execution['stdout_sha256']
        assert digest(stderr) == execution['stderr_sha256']
        rr = rows(data)
        metric = report['profiles'][name]
        assert len(rr) == metric['epochs']
        assert dict(collections.Counter(r['status'] for r in rr)) == metric['statuses']
        assert sum(r['status']=='PPP_FLOAT' and not r['phase_rms_m'] for r in rr) == metric['PPP_code_only']
        assert sum(r['status']=='PPP_FLOAT' and bool(r['phase_rms_m']) for r in rr) == metric['PPP_with_phase']
        assert sum(r['status']=='PPP_FLOAT' and float(r['zenith_m'])<0 for r in rr) == metric['negative_PPP_ZTD']
        total += len(rr)
        city, run, pair, profile = name.split('-')
        case = city+'-'+run
        with Path(old_manifest['cases'][case]['reference']).open(newline='') as f:
            ref = [{k.strip():v for k,v in r.items()} for r in csv.DictReader(f)]
        refs = {(int(r['GPS Week']),Decimal(r['GPS TOW (s)'])):r for r in ref}
        nr = normalize(rr)
        view = [dict(r,status='FLOAT') if r['status']=='PPP_FLOAT' else r for r in nr]
        assert scorer(view,refs,city) == metric['reference_errors']
        control_path = old/(case+'-'+pair+'-default.stdout.csv.gz')
        control_data = gzip.decompress(control_path.read_bytes())
        if profile == 'control':
            assert data == control_data and metric['all_old_72_columns_byte_identical']
        else:
            assert compare(nr,normalize(rows(control_data)),refs,city,scorer) == metric['frozen_Rust_comparison']
            if pair == '1C2W':
                for seed, expected in metric['frozen_native_comparison'].items():
                    native_rows = normalize(rows(gzip.decompress((frozen/(case+'-'+seed+'.stdout.csv.gz')).read_bytes())))
                    assert compare(nr,native_rows,refs,city,scorer) == expected
    assert total == 69036
    assert json.loads((root/'checks.json').read_text())['passed_tests_including_doctests'] == 456
    print('Verified 69,036 retained epochs, 32 native component cases, old CSV control and frozen comparisons; no practical acceptance')


if __name__ == '__main__':
    main()
