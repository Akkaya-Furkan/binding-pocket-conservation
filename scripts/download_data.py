"""Download public source snapshots; write provenance after both fetches succeed."""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]
QUERY = '(gene_exact:CA2) AND (reviewed:true) AND (taxonomy_id:7742)'
SOURCES = {
    'ca2_reviewed.fasta': 'https://rest.uniprot.org/uniprotkb/stream?' + urlencode({'query': QUERY, 'format': 'fasta'}),
    '3HS4.pdb': 'https://files.rcsb.org/download/3HS4.pdb',
}


def main():
    # Fetch everything before mutating the snapshot directory.
    payloads = {name: urlopen(url, timeout=120).read() for name, url in SOURCES.items()}
    if not payloads['ca2_reviewed.fasta'].startswith(b'>') or b'ATOM  ' not in payloads['3HS4.pdb']:
        raise ValueError('Unexpected source format; snapshot not changed')
    dest = ROOT / 'data/raw'
    dest.mkdir(parents=True, exist_ok=True)
    provenance = []
    for name, data in payloads.items():
        (dest / name).write_bytes(data)
        provenance.append(dict(file=name, url=SOURCES[name], bytes=len(data),
                               sha256=hashlib.sha256(data).hexdigest(),
                               retrieved_utc=datetime.now(timezone.utc).isoformat()))
    (dest / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')


if __name__ == '__main__':
    main()
