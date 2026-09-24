from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import os
from urllib.parse import quote
from urllib.request import urlopen

from flask import Flask, jsonify, render_template, request

from helix.core import GENETIC_CODES, analyze_six_frames

app = Flask(__name__)

@app.get('/')
def index():
    return render_template('index.html', genetic_codes=[(k, v[0]) for k,v in GENETIC_CODES.items()])

@app.post('/api/analyze')
def analyze():
    payload = request.get_json(silent=True) or {}
    try:
        result = analyze_six_frames(
            payload.get('sequence',''),
            genetic_code=str(payload.get('genetic_code','1')),
            min_orf_length_nt=int(payload.get('min_orf_length_nt', payload.get('min_orf', 60))),
        )
        return jsonify(result)
    except (ValueError, TypeError) as exc:
        return jsonify({'error': str(exc)}), 400

@app.get('/api/fetch-ncbi')
def fetch_ncbi():
    accession = request.args.get('accession','').strip()
    if not accession:
        return jsonify({'error':'Accession is required.'}), 400
    url = (
        'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
        f'?db=nuccore&id={quote(accession)}&rettype=fasta&retmode=text'
    )
    try:
        with urlopen(url, timeout=15) as response:
            fasta = response.read().decode('utf-8')
    except Exception as exc:
        return jsonify({'error':f'NCBI request failed: {exc}'}), 502
    if not fasta.startswith('>'):
        return jsonify({'error':'NCBI did not return FASTA for that accession.'}), 404
    return jsonify({'fasta':fasta, 'title':fasta.splitlines()[0][1:80] if fasta.splitlines() else accession})

if __name__ == '__main__':
    app.run(host='127.0.0.1', port=int(os.getenv('PORT','3030')), debug=False)
