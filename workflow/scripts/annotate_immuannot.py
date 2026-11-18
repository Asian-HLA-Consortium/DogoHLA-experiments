#!/usr/bin/env python3
"""Wrapper tool to run Immuannot on phased assembly FASTA to generate HLA annotations GTF.

Invokes frozen Immuannot (commit 31362c370123f6706a04ffdda6a3bf0c108328e0) using
pinned reference bundle (Data-ipd3590-kir2130).
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

# Add workflow to sys.path to access dogohla.verify
SCRIPT_DIR = Path(__file__).resolve().parent
WORKFLOW_DIR = SCRIPT_DIR.parent
if str(WORKFLOW_DIR) not in sys.path:
    sys.path.insert(0, str(WORKFLOW_DIR))

try:
    from dogohla.verify import verify
except ImportError:
    verify = None


def find_immuannot_sh() -> Path:
    """Locate frozen immuannot.sh script."""
    candidate = WORKFLOW_DIR / "dogohla" / "frozen" / "Immuannot" / "scripts.pub.v3" / "immuannot.sh"
    if candidate.is_file():
        return candidate
    env_path = os.environ.get("IMMUANNOT_SH")
    if env_path and Path(env_path).is_file():
        return Path(env_path)
    raise FileNotFoundError(f"Could not locate immuannot.sh at {candidate} or via IMMUANNOT_SH")


def main():
    parser = argparse.ArgumentParser(description="Run Immuannot annotation on assembly FASTA")
    parser.add_argument("--assembly", required=True, type=Path, help="Phased assembly or span FASTA")
    parser.add_argument("--immuannot-ref", required=True, type=Path, help="Reference data bundle (e.g. Data-ipd3590-kir2130)")
    parser.add_argument("--sample", required=True, help="Sample/donor ID")
    parser.add_argument("--haplotype", required=True, help="Haplotype ID (1 or 2)")
    parser.add_argument("--threads", type=int, default=4, help="Threads (default: 4)")
    parser.add_argument("--out-gtf", type=Path, default=None, help="Output GTF.gz path")
    args = parser.parse_args()

    # 1. Verify frozen code integrity
    if verify is not None:
        v_res = verify()
        print(f"[annotate_immuannot] Frozen manifest verified ({v_res['files_checked']} files).")

    immuannot_sh = find_immuannot_sh()
    print(f"[annotate_immuannot] Using immuannot.sh at: {immuannot_sh}")

    assembly = args.assembly.resolve()
    refdir = args.immuannot_ref.resolve()
    if not assembly.is_file():
        raise FileNotFoundError(f"Assembly FASTA not found: {assembly}")
    if not refdir.is_dir():
        raise FileNotFoundError(f"Immuannot reference dir not found: {refdir}")

    # Determine out prefix and final gtf destination
    cwd = Path.cwd()
    out_pref = cwd / f"{args.sample}.{args.haplotype}"
    expected_gtf = Path(f"{out_pref}.gtf.gz")

    final_gtf = args.out_gtf.resolve() if args.out_gtf else expected_gtf

    cmd = [
        "bash",
        str(immuannot_sh),
        "-c", str(assembly),
        "-r", str(refdir),
        "-o", str(out_pref),
        "-t", str(args.threads)
    ]

    print(f"[annotate_immuannot] Running command: {' '.join(cmd)}")
    sys.stdout.flush()

    res = subprocess.run(cmd, check=False)
    if res.returncode != 0:
        print(f"[annotate_immuannot] Immuannot exited with error code {res.returncode}", file=sys.stderr)
        sys.exit(res.returncode)

    if not expected_gtf.is_file() or expected_gtf.stat().st_size == 0:
        print(f"[annotate_immuannot] Expected output GTF missing or empty: {expected_gtf}", file=sys.stderr)
        sys.exit(1)

    if final_gtf != expected_gtf:
        final_gtf.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(expected_gtf, final_gtf)
        print(f"[annotate_immuannot] Copied {expected_gtf} to {final_gtf}")

    print(f"[annotate_immuannot] Success: {final_gtf} ({final_gtf.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
