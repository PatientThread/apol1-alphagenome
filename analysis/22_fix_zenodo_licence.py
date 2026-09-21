#!/usr/bin/env python3
"""
22_fix_zenodo_licence.py

Changes the licence on the published Zenodo record from CC-BY-4.0 to
CC-BY-NC-4.0, and sets the description that states the licence split.

WHY. Zenodo guessed CC-BY-4.0 from the repository. CC-BY permits commercial
reuse with attribution. The AlphaGenome predictions in this archive, and
everything derived from them, are non-commercial under Google DeepMind's
AlphaGenome Output Terms of Use, so CC-BY asserts a right the depositor does
not hold and contradicts LEGALLY_BINDING_TERMS_OF_USE.txt inside the archive.

Editing metadata does NOT mint a new DOI. 10.5281/zenodo.22875737 is unchanged.

TOKEN. Create one at https://zenodo.org/account/settings/applications/tokens/new/
with scopes `deposit:write` and `deposit:actions`. Store it the way the other
keys in this portfolio are stored, outside any synced directory:

    mkdir -p ~/.zenodo && printf '%s' 'YOUR_TOKEN' > ~/.zenodo/token
    chmod 600 ~/.zenodo/token

Run with --dry-run first to see what would change.

Author: Christopher Lawrence
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path

RECORD_ID = "22875737"
API = "https://zenodo.org/api"
KEY_PATH = Path.home() / ".zenodo" / "token"
NEW_LICENCE = "cc-by-nc-4.0"

DESCRIPTION = (
    "<p>Analysis code and frozen results supporting &quot;Reference choice "
    "affects positional enrichment of AlphaGenome variant rankings at "
    "APOL1-MYH9&quot;.</p>"
    "<p><strong>Licensing is split.</strong> The analysis code is released "
    "under the MIT licence (see LICENSE). The AlphaGenome predictions in "
    "results/scored*, and every file derived from them, are AlphaGenome Output "
    "governed by Google DeepMind's AlphaGenome Output Terms of Use "
    "(<a href=\"http://deepmind.google.com/science/alphagenome/output-terms\">"
    "output terms</a>), which are non-commercial. The record-level licence is "
    "set to CC-BY-NC-4.0 so that it does not assert commercial rights over "
    "those outputs. LEGALLY_BINDING_TERMS_OF_USE.txt in the archive lists "
    "exactly which files fall under which terms; any file not listed should be "
    "treated as AlphaGenome Output. Model weights are not included and were "
    "never accessed. ENCODE peak files and gnomAD frequencies carry their own "
    "source terms.</p>"
)


def token() -> str:
    env = os.environ.get("ZENODO_TOKEN", "").strip()
    if env:
        return env
    if KEY_PATH.exists():
        mode = KEY_PATH.stat().st_mode & 0o777
        if mode & 0o077:
            raise SystemExit(f"{KEY_PATH} is mode {mode:o}, readable beyond "
                             f"the owner. Run: chmod 600 {KEY_PATH}")
        t = KEY_PATH.read_text().strip()
        if t:
            return t
    raise SystemExit(
        "No Zenodo token found.\n"
        f"  Expected in $ZENODO_TOKEN or {KEY_PATH} (mode 600).\n"
        "  Create one at "
        "https://zenodo.org/account/settings/applications/tokens/new/\n"
        "  with scopes deposit:write and deposit:actions.")


def call(method: str, url: str, tok: str, payload=None):
    data = json.dumps(payload).encode() if payload is not None else None
    req = urllib.request.Request(url, data=data, method=method, headers={
        "Authorization": f"Bearer {tok}",
        "Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=60) as fh:
            body = fh.read().decode()
            return json.loads(body) if body else {}
    except urllib.error.HTTPError as e:
        raise SystemExit(f"{method} {url} -> HTTP {e.code}\n{e.read().decode()[:600]}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    tok = token()
    rec = call("GET", f"{API}/records/{RECORD_ID}", tok)
    cur = rec.get("metadata", {}).get("license") or rec.get("metadata", {}).get("rights")
    print(f"record   : {RECORD_ID}  doi {rec.get('doi')}")
    print(f"title    : {rec.get('metadata', {}).get('title')}")
    print(f"licence  : {cur}  ->  {NEW_LICENCE}")
    if a.dry_run:
        print("\n--dry-run: nothing sent.")
        return

    # Open a draft of the published record, patch it, publish again.
    draft = call("POST", f"{API}/records/{RECORD_ID}/draft", tok)
    md = draft.get("metadata", {})
    md["rights"] = [{"id": NEW_LICENCE}]
    md["description"] = DESCRIPTION
    call("PUT", f"{API}/records/{RECORD_ID}/draft", tok, {**draft, "metadata": md})
    out = call("POST", f"{API}/records/{RECORD_ID}/draft/actions/publish", tok)
    print(f"\npublished. doi is still {out.get('doi')}")
    print("verify: https://zenodo.org/records/" + RECORD_ID)


if __name__ == "__main__":
    main()
