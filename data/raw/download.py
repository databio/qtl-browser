#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = ["pyyaml"]
# ///
"""Download raw data sources listed in sources.yaml into <dest>/<source dir>/.

    ./download.py                 # fetch everything not yet present and verified
    ./download.py --only dbsnp    # one source (repeatable)
    ./download.py --list          # show status without downloading
    ./download.py --config other.yaml

Downloads resume (curl -C -) only when the server proves it honours Range,
md5 is verified when the config gives one, and a file that another process is
currently writing is skipped rather than clobbered.
"""
import argparse
import hashlib
import os
import subprocess
import sys
import urllib.request
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent


def md5_of(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def being_written(path: Path) -> bool:
    """True if some other process has this path open for writing (macOS/Linux)."""
    try:
        out = subprocess.run(["lsof", "-F", "a", "--", str(path)], capture_output=True, text=True).stdout
    except FileNotFoundError:
        return False
    return any(line.startswith("a") and ("w" in line or "u" in line) for line in out.splitlines())


def resolve(src: dict, entry: dict) -> tuple[str, str]:
    """Return (url, basename) for a file entry."""
    if "url" in entry:
        return entry["url"], entry.get("file") or entry["url"].rstrip("/").rsplit("/", 1)[-1]
    if "url_base" not in src:
        sys.exit(f"{src['name']}: file entry {entry} has no url and source has no url_base")
    return src["url_base"].format(file=entry["file"]), entry["file"]


def expected_md5(entry: dict) -> str | None:
    if "md5" in entry:
        return entry["md5"].lower()
    if "md5_url" in entry:
        with urllib.request.urlopen(entry["md5_url"], timeout=60) as r:
            return r.read().decode().split()[0].lower()
    return None


def status(dest: Path, entry: dict, md5: str | None) -> str:
    if not dest.exists():
        return "missing"
    if being_written(dest):
        return "in-progress"
    size = entry.get("size")
    if size is not None:
        have = dest.stat().st_size
        if have > size:
            return "oversize"
        if have < size:
            return "partial"
    if md5 is not None:
        return "verified" if md5_of(dest) == md5 else "md5-mismatch"
    return "present"


def honours_range(url: str) -> bool:
    """True if the server answers a one byte Range request with 206.

    We ask for one real byte instead of reading `Accept-Ranges` from a HEAD.
    The header is only advisory, and these URLs redirect, so the host that
    answers the HEAD is not always the host that serves the body.
    """
    out = subprocess.run(
        ["curl", "-sL", "-o", os.devnull, "-w", "%{http_code}", "--max-time", "120", "-r", "0-0", url],
        capture_output=True,
        text=True,
    ).stdout.strip()
    return out.endswith("206")


def curl(url: str, dest: Path, resume: bool) -> None:
    cmd = ["curl", "-L", "--fail", "--retry", "5", "--retry-delay", "10"]
    if resume:
        cmd += ["-C", "-"]
    cmd += ["-o", str(dest), url]
    subprocess.run(cmd, check=True)


def fetch(url: str, dest: Path, size: int | None = None) -> None:
    """Download url to dest, resuming only when resuming is safe.

    Do not drop this guard as redundant. On 2026-09-22 a 9.4 GB Zenodo file
    died at 4.2 GB with curl exit 18. The next run resumed it with `-C -`,
    Zenodo ignored the Range header and answered 200 with the whole file from
    byte 0, and curl appended those bytes onto the 4.2 GB already on disk. The
    result was a 13.4 GB file that still looked merely "partial", so the run
    after that would have appended again. So: never append onto a file that is
    already at or past its expected size, only resume when the server proves
    it honours Range, and treat any overshoot as corruption.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    have = dest.stat().st_size if dest.exists() else 0
    if have and size is not None and have >= size:
        print(f"  discarding    {dest.name} ({have} bytes, expected {size})")
        dest.unlink()
        have = 0
    if have and not honours_range(url):
        print(f"  no range      {dest.name}, server will not resume, starting over")
        dest.unlink()
        have = 0
    curl(url, dest, resume=bool(have))
    if size is not None and dest.exists() and dest.stat().st_size > size:
        # The server said 206 to the probe but served the whole body anyway.
        print(f"  overshot      {dest.name}, discarding and fetching from zero")
        dest.unlink()
        curl(url, dest, resume=False)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=HERE / "sources.yaml")
    ap.add_argument("--only", action="append", metavar="SOURCE", help="restrict to this source name (repeatable)")
    ap.add_argument("--list", action="store_true", help="report status only, download nothing")
    args = ap.parse_args()

    cfg = yaml.safe_load(args.config.read_text())
    root = args.config.resolve().parent
    # `dest` in the config is relative to the repo root, which is two levels above data/raw/sources.yaml
    dest_root = (root.parents[1] / cfg["dest"]).resolve() if not Path(cfg["dest"]).is_absolute() else Path(cfg["dest"])

    failures = 0
    for src in cfg["sources"]:
        if args.only and src["name"] not in args.only:
            continue
        print(f"== {src['name']}: {src.get('version', '')}")
        for entry in src["files"]:
            url, name = resolve(src, entry)
            dest = dest_root / src["dir"] / name
            md5 = expected_md5(entry)
            st = status(dest, entry, md5)
            if args.list or st in ("verified", "present", "in-progress"):
                print(f"  {st:13s} {name}")
                continue
            if st == "md5-mismatch":
                print(f"  md5 mismatch, re-downloading {name}")
                dest.unlink()
            if st == "oversize":
                print(f"  bigger than the {entry['size']} bytes expected, re-downloading {name}")
                dest.unlink()
            print(f"  fetching      {name}  <- {url}")
            try:
                fetch(url, dest, entry.get("size"))
            except subprocess.CalledProcessError as e:
                print(f"  FAILED        {name} (curl exit {e.returncode})")
                failures += 1
                continue
            st = status(dest, entry, md5)
            print(f"  {st:13s} {name}")
            if st not in ("verified", "present"):
                print(f"  FAILED        {name} ({st} after downloading)")
                failures += 1
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
