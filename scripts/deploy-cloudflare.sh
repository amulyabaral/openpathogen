#!/usr/bin/env bash
# Deploy openpathogen to Cloudflare Pages (free tier, no card needed).
# Files over the Pages 25 MiB per-file cap ship as chunks that the browser
# reassembles — see scripts/build-dist.sh and js/assets.js.
# Prerequisite (once): npx wrangler login
set -euo pipefail
cd "$(dirname "$0")/.."

./scripts/build-dist.sh

echo "→ Cloudflare Pages"
npx wrangler pages project create openpathogen --production-branch main 2>/dev/null || true
# no positional dir: wrangler.toml's pages_build_output_dir supplies it
npx wrangler pages deploy --project-name openpathogen --branch main --commit-dirty=true

echo
echo "Site: https://openpathogen.pages.dev"
