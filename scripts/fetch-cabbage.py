#!/usr/bin/env python3
"""Fetch a complete CABBAGE snapshot from the EMBL-EBI AMR portal API.

Pages POST /amr/api/amr-records for every view and writes one JSONL file per
view into databases/cabbage/raw/. The portal's /amr-records/download CSV
silently drops rows (1.34M of 1.71M for the phenotype view at release
2026-07), so paging the JSON API is the only faithful source.

Usage: python3 scripts/fetch-cabbage.py [--per-page 1000] [--views 1,2,3]
Polite by default: 4 concurrent requests, 0.2s spacing. ~3,200 requests for
the full snapshot (~10 minutes).
"""
import argparse
import json
import os
import sys
import time
import urllib.request
import urllib.error
from concurrent.futures import ThreadPoolExecutor

API = "https://www.ebi.ac.uk/amr/api/amr-records"
RELEASE_URL = "https://www.ebi.ac.uk/amr/api/release"
RAW_DIR = os.path.join(os.path.dirname(__file__), "..", "databases", "cabbage", "raw")

VIEW_NAMES = {1: "phenotypes", 2: "genotypes", 3: "combined"}


def post_json(url, payload, timeout=120):
    req = urllib.request.Request(
        url, data=json.dumps(payload).encode(), headers={"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read())


def fetch_page(view_id, page, per_page, retries=5):
    payload = {"view_id": view_id, "page": page, "per_page": per_page}
    for attempt in range(retries):
        try:
            return post_json(API, payload)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-page", type=int, default=1000)
    ap.add_argument("--views", default="1,2,3")
    ap.add_argument("--concurrency", type=int, default=4)
    args = ap.parse_args()

    release = json.loads(urllib.request.urlopen(RELEASE_URL, timeout=30).read())["label"]
    os.makedirs(RAW_DIR, exist_ok=True)
    print(f"release: {release}")

    for vid in [int(v) for v in args.views.split(",")]:
        name = VIEW_NAMES[vid]
        out_path = os.path.join(RAW_DIR, f"{name}.jsonl")
        if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
            print(f"{name}: {out_path} already exists, skipping (delete to refetch)")
            continue

        first = fetch_page(vid, 1, args.per_page)
        total = first["meta"]["total_hits"]
        columns = [c["id"] for c in first["meta"]["columns"]]
        pages = (total + args.per_page - 1) // args.per_page
        print(f"{name}: {total:,} rows, {pages:,} pages")

        # Meta for this view travels beside the rows.
        with open(out_path + ".meta.json", "w") as mf:
            json.dump({"view_id": vid, "name": name, "release": release,
                       "total_hits": total, "columns": columns}, mf)

        done = 0
        with open(out_path, "w") as out, ThreadPoolExecutor(args.concurrency) as pool:
            def grab(page):
                return page, fetch_page(vid, page, args.per_page)

            results = {}
            for page, res in pool.map(grab, range(1, pages + 1)):
                results[page] = res["data"]
                while done + 1 in results:  # write strictly in order
                    done += 1
                    for row in results.pop(done):
                        out.write(json.dumps(row, separators=(",", ":")) + "\n")
                if done % 50 == 0:
                    print(f"  {name}: {done}/{pages} pages", flush=True)
                time.sleep(0.2)  # per-response spacing

        rows_written = sum(1 for _ in open(out_path))
        print(f"{name}: wrote {rows_written:,} rows to {out_path}")
        if rows_written != total:
            print(f"WARNING: {name} row count mismatch (api={total}, wrote={rows_written})",
                  file=sys.stderr)


if __name__ == "__main__":
    main()
