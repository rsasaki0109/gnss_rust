"""Run the source-ported single-update path on PPC; no parameter search."""
import argparse
import collections
import csv
from decimal import Decimal
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import runpy
import subprocess
import time


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def rows(path):
    return list(csv.DictReader(io.StringIO(gzip.decompress(path.read_bytes()).decode())))


def normalized(rr):
    result = []
    seen = set()
    for r in rr:
        tow = Decimal(r['tow']).quantize(Decimal('.000001'))
        assert abs(tow-Decimal(r['tow'])) <= Decimal('.000001')
        key = (int(r['week']), tow)
        assert key not in seen
        seen.add(key)
        result.append(dict(r, week=str(key[0]), tow=str(tow)))
    return result


def save(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir', type=Path, required=True)
    a = p.parse_args()
    repo = Path(__file__).resolve().parents[1]
    out = a.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    assert not (out/'predeclared-inputs.json').exists()
    old = repo/'docs/comparison/2026-10-21-ppp-ztd-prior'
    prior = json.loads((old/'predeclared-inputs.json').read_text())
    native = repo/'docs/comparison/2026-10-21-native-ppp-common-inputs'
    compare = runpy.run_path(str(repo/'scripts/run_ppc_ppp_initial_prior.py'))['compare']
    scoring = runpy.run_path(str(repo/'scripts/run_ppc_rtk_float.py'))
    scorer, targets = scoring['score'], scoring['reference_targets']
    profiles = [(c, '1C2W', True) for c in prior['cases']]
    profiles += [('nagoya-run3', '1C2L', True), ('nagoya-run3', '1C2L', False)]
    sources = list((repo/'src').rglob('*.rs')) + [repo/f for f in [
        'examples/ppc_ppp.rs', 'examples/support/ppp_seed.rs',
        'examples/support/ppp_ztd_prior.rs', 'scripts/run_ppp_reference_update.py',
        'scripts/run_ppc_ppp_initial_prior.py', 'scripts/run_ppc_rtk_float.py',
        'Cargo.toml', 'Cargo.lock', 'rust-toolchain.toml']]
    commands = {}
    for case, pair, enabled in profiles:
        name = case+'-'+pair+('-reference' if enabled else '-control')
        cmd = json.loads((old/(case+'-'+pair+'-default.execution.json')).read_text())['command']
        commands[name] = cmd+(['-native-update'] if enabled else [])
    manifest = dict(
        base_commit=subprocess.check_output(['git','rev-parse','HEAD'], cwd=repo, text=True).strip(),
        source_sha256={str(f.relative_to(repo)): sha(f) for f in sources},
        binary_sha256=sha(commands[next(iter(commands))][0]),
        input_sha256=prior['input_sha256'], commands=commands,
        frozen_control_archive=str(old.relative_to(repo)),
        frozen_native_archive=str(native.relative_to(repo)),
        frozen_native_evidence_sha256=sha(native/'evidence-sha256.json'),
        initial_mean_ambiguity_and_measurement_update_changed=True,
        source_profile='native precise IF initialization, phase lock 1, one update, native prefit and 4-sigma screen',
        lifecycle='Rust accepted-only numeric time and transactional phase arcs; not yet native',
        native_feature_match=False, whole_native_PPP_executed=False,
        runtime_memory_acceptance=False, practical_parity_accepted=False,
        reference_used_by_solver=False,
    )
    assert all(sha(f)==h for f,h in manifest['input_sha256'].items())
    save(out/'predeclared-inputs.json', manifest)
    report = dict(profiles={}, practical_parity_accepted=False, native_feature_match=False,
                  no_new_native_whole_processor_runs=True)
    for (case, pair, enabled), (name, cmd) in zip(profiles, commands.items()):
        assert all(sha(repo/f)==h for f,h in manifest['source_sha256'].items())
        assert sha(cmd[0])==manifest['binary_sha256']
        started = time.monotonic()
        r = subprocess.run(cmd, capture_output=True, cwd=repo)
        for label, data in [('stdout.csv',r.stdout), ('stderr.log',r.stderr)]:
            (out/(name+'.'+label+'.gz')).write_bytes(gzip.compress(data, mtime=0))
        save(out/(name+'.execution.json'), dict(command=cmd, returncode=r.returncode,
             elapsed_seconds=time.monotonic()-started,
             stdout_sha256=hashlib.sha256(r.stdout).hexdigest(),
             stderr_sha256=hashlib.sha256(r.stderr).hexdigest()))
        assert r.returncode==0, r.stderr[-1000:].decode()
        rr = rows(out/(name+'.stdout.csv.gz'))
        assert all(None not in row for row in rr)
        nr = normalized(rr)
        assert len(nr)==len(rows(old/(case+'-'+pair+'-default.stdout.csv.gz')))
        with Path(prior['cases'][case]['reference']).open(newline='') as f:
            ref = [{k.strip():v for k,v in row.items()} for row in csv.DictReader(f)]
        references = {(int(row['GPS Week']),Decimal(row['GPS TOW (s)'])):row for row in ref}
        assert len(references)==len(ref)
        metric = dict(epochs=len(rr), statuses=dict(collections.Counter(row['status'] for row in rr)),
                      errors=dict(collections.Counter(row['error'] for row in rr if row['error'])),
                      PPP_code_only=sum(row['status']=='PPP_FLOAT' and not row['phase_rms_m'] for row in rr),
                      PPP_with_phase=sum(row['status']=='PPP_FLOAT' and bool(row['phase_rms_m']) for row in rr),
                      negative_PPP_ZTD=sum(row['status']=='PPP_FLOAT' and float(row['zenith_m'])<0 for row in rr))
        view = [dict(row,status='FLOAT') if row['status']=='PPP_FLOAT' else row for row in nr]
        metric['reference_errors'] = scorer(view,references,prior['cases'][case]['city'])
        if enabled:
            control = normalized(rows(old/(case+'-'+pair+'-default.stdout.csv.gz')))
            metric['frozen_Rust_comparison'] = compare(nr,control,references,prior['cases'][case]['city'],scorer)
            if pair=='1C2W':
                metric['frozen_native_comparison'] = {}
                for seed in ['native_builtin','native_header_zero']:
                    frozen = normalized(rows(native/(case+'-'+seed+'.stdout.csv.gz')))
                    metric['frozen_native_comparison'][seed] = compare(nr,frozen,references,prior['cases'][case]['city'],scorer)
            worst = {}
            for row in nr:
                if row['status']!='PPP_FLOAT':
                    continue
                refrow = references[(int(row['week']),Decimal(row['tow']))]
                hypotheses,_,_ = targets(refrow,prior['cases'][case]['city'])
                for hypothesis, target in hypotheses.items():
                    error = math.sqrt(sum((float(row[k])-v)**2 for k,v in zip(['x_m','y_m','z_m'],target)))
                    if hypothesis not in worst or worst[hypothesis]['error_3d_m'] < error:
                        worst[hypothesis] = dict(error_3d_m=error, fields=row)
            metric['worst_PPP_by_reference_hypothesis'] = worst
        else:
            assert r.stdout==gzip.decompress((old/(case+'-'+pair+'-default.stdout.csv.gz')).read_bytes())
            metric['all_old_72_columns_byte_identical'] = True
        report['profiles'][name] = metric
        save(out/'report.json', report)
        print(name,metric['statuses'],metric['reference_errors']['hypotheses'].get('reference_is_antenna'),flush=True)
    assert all(sha(f)==h for f,h in manifest['input_sha256'].items())


if __name__=='__main__':
    main()
