# Deploy HELIX on Vercel

HELIX is a Flask application and is prepared for Vercel's Python/Flask runtime.

## Git-based deployment

1. Push this directory to a GitHub repository named `helix-dna-rna-translator`.
2. In Vercel, choose **Add New → Project** and import that repository.
3. Vercel should detect **Flask** automatically.
4. Keep the repository root as the project root and deploy.
5. After deployment, test:
   - `/`
   - `POST /api/analyze`
   - `/api/fetch-ncbi?accession=NM_000546.6`

No secret is required for the core application. NCBI fetching depends on outbound network availability.
