# Deploying openpathogen.org (Cloudflare Pages, free tier — no card)

Everything — site, WASM modules, database indexes, example dataset — is served
from **one Cloudflare Pages project** with unlimited bandwidth.

Pages caps individual files at 25 MiB, so `scripts/build-dist.sh` ships the
three files over that cap (the 74 MB VFDB index and the two ~25 MB subsampled
example FASTQs) as `<name>.partNN` chunks; `js/assets.js` tries the plain URL
first (local dev against the repo root works unchanged) and otherwise fetches
the parts and reassembles them in the browser.

## One-time setup

1. **Cloudflare account** (free, no payment method) —
   https://dash.cloudflare.com/sign-up
2. **Login from this repo** — `npx wrangler login` (opens the browser).
3. **Add your domain to Cloudflare** — dashboard → Add a domain →
   `openpathogen.org` → Free plan. Cloudflare shows two nameservers.
4. **Namecheap** — Domain List → Manage → Nameservers → *Custom DNS* →
   paste the two Cloudflare nameservers. Registration stays at Namecheap;
   only DNS moves. Propagation usually takes minutes–hours.

## Deploy

```bash
./scripts/deploy-cloudflare.sh
```

Builds `dist/` (chunking oversized files) and deploys. Live at
https://openpathogen.pages.dev immediately.

## One-time domain attach (after nameservers point at Cloudflare)

Dashboard → Workers & Pages → openpathogen → Custom domains → add
`openpathogen.org` and `www.openpathogen.org`. DNS records and SSL
certificates are provisioned automatically since the zone is on Cloudflare.

## Costs

$0. Free-tier Pages includes unlimited bandwidth and unlimited requests;
nothing here needs a payment method.
