"""Acquire IGS final GPS products and CAS DSBs for the three PPC observation dates."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date, datetime, timedelta
import gzip
import hashlib
import json
from pathlib import Path
import subprocess


def sha(data):
    return hashlib.sha256(data).hexdigest()


def acquire(root, day, kind):
    stamp = day.strftime("%Y%j")
    names = {
        "orbit": f"IGS0OPSFIN_{stamp}0000_01D_15M_ORB.SP3.gz",
        "clock": f"IGS0OPSFIN_{stamp}0000_01D_30S_CLK.CLK.gz",
        "dsb": f"CAS0OPSRAP_{stamp}0000_01D_01D_DCB.BIA.gz",
    }
    filename = names[kind]
    week = (day - date(1980, 1, 6)).days // 7
    url = f"https://igs.bkg.bund.de/root_ftp/IGS/products/{week}/{filename}"
    compressed = root / filename
    headers = root / (filename + ".headers")
    result = subprocess.run([
        "curl", "--location", "--fail", "--silent", "--show-error",
        "--connect-timeout", "10", "--max-time", "60",
        "--dump-header", str(headers), "--output", str(compressed), url,
    ], capture_output=True)
    item = dict(date=str(day), kind=kind, filename=filename, url=url,
                returncode=result.returncode, acquired=False,
                stderr=result.stderr.decode())
    if result.returncode:
        return item
    data = compressed.read_bytes()
    try:
        with gzip.open(compressed, "rb") as stream:
            decoded = stream.read(128 * 1024 * 1024 + 1)
        if len(decoded) > 128 * 1024 * 1024:
            raise ValueError("decompressed product exceeds 128 MiB")
        lines = decoded.decode("ascii").splitlines()
        if kind == "orbit":
            assert lines[0].startswith("#")
            first_epoch = next(line for line in lines if line.startswith("*"))
            assert tuple(map(int, first_epoch.split()[1:4])) == (day.year, day.month, day.day)
            records = sum(line.startswith("PG") for line in lines)
            assert records > 0 and any(line.strip() == "EOF" for line in lines)
        elif kind == "clock":
            assert "RINEX VERSION / TYPE" in lines[0]
            first_epoch = next(line for line in lines if line.startswith("AS G"))
            assert tuple(map(int, first_epoch.split()[2:5])) == (day.year, day.month, day.day)
            records = sum(line.startswith("AS G") for line in lines)
            assert records > 0 and any("END OF HEADER" in line for line in lines)
        else:
            assert lines[0].startswith("%=BIA")
            assert any(line.startswith("-BIAS/SOLUTION") for line in lines)
            records = sum(line.startswith(" DSB") for line in lines)
            assert records > 0
        expanded = compressed.with_suffix("")
        expanded.write_bytes(decoded)
        item.update(acquired=True, compressed_bytes=len(data), compressed_sha256=sha(data),
                    expanded_filename=expanded.name, expanded_bytes=len(decoded),
                    expanded_sha256=sha(decoded), relevant_records=records,
                    first_line=lines[0])
    except (AssertionError, ValueError, OSError, EOFError, StopIteration) as error:
        item["validation_error"] = type(error).__name__ + ": " + str(error)
    return item


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    centres = [date(2024, 7, 20), date(2024, 7, 23), date(2024, 8, 3)]
    days = sorted({day + timedelta(days=offset) for day in centres for offset in [-1, 0, 1]})
    source = Path(__file__).resolve()
    manifest = dict(
        source_sha256=sha(source.read_bytes()), central_observation_dates=list(map(str, centres)),
        padding_days_each_side=1, source="BKG public IGS product mirror",
        orbit_family="IGS final GPS, 15-minute SP3",
        clock_family="IGS final GPS, 30-second CLK",
        bias_family="CAS rapid DSB; not CODE OSB",
        code_osb_acquired=False, receiver_clock_signal_reference_configured=False,
        applied_to_ppp=False, ppc_performance_accepted=False,
        started_utc=datetime.now().isoformat(timespec="seconds") + "Z", products=[],
    )
    path = root / "acquisition.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures = [pool.submit(acquire, root, day, kind)
                   for day in days for kind in ["orbit", "clock", "dsb"]]
        for future in as_completed(futures):
            item = future.result()
            manifest["products"].append(item)
            path.write_text(json.dumps(manifest, indent=2) + "\n")
            print(json.dumps({k: item[k] for k in ["filename", "returncode", "acquired"]}), flush=True)
    manifest["products"].sort(key=lambda item: (item["date"], item["kind"]))
    manifest["completed_utc"] = datetime.now().isoformat(timespec="seconds") + "Z"
    manifest["acquired_count"] = sum(p["acquired"] for p in manifest["products"])
    assert sha(source.read_bytes()) == manifest["source_sha256"]
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Acquired {manifest['acquired_count']}/27 structurally validated products.")
    raise SystemExit(0 if manifest["acquired_count"] == 27 else 1)


if __name__ == "__main__":
    main()
