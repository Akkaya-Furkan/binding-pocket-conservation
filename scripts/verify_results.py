"""Verify regenerated scientific outputs against the committed snapshot."""
import csv
import io
import json
import math
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def committed(name):
    return subprocess.check_output(['git', 'show', f'HEAD:results/{name}'], cwd=ROOT, text=True)


def compare(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for key in a:
            compare(a[key], b[key])
    elif isinstance(a, (int, float)):
        assert math.isclose(a, b, rel_tol=1e-8, abs_tol=1e-8), (a, b)
    else:
        assert a == b, (a, b)


def main():
    expected = json.loads(committed('summary.json'))
    actual = json.loads((ROOT / 'results/summary.json').read_text())
    # Patch versions can differ; scientific results must agree.
    expected.pop('versions')
    actual.pop('versions')
    compare(expected, actual)
    names = ['sequences.csv', 'residues.csv', 'pocket.csv', 'sensitivity.csv', 'permutation_null.csv']
    for name in names:
        left = list(csv.DictReader(io.StringIO(committed(name))))
        right = list(csv.DictReader((ROOT / 'results' / name).open()))
        assert len(left) == len(right), name
        for a, b in zip(left, right):
            assert a.keys() == b.keys()
            for key in a:
                try:
                    x, y = float(a[key]), float(b[key])
                except ValueError:
                    assert a[key] == b[key], (name, key)
                else:
                    compare(x, y)
    for name in ['alignment.fasta', 'input.fasta']:
        assert committed(name) == (ROOT / 'results' / name).read_text(), name
    for name in ['comparison.png', 'conservation_profile.png']:
        assert (ROOT / 'results' / name).read_bytes().startswith(b'\x89PNG\r\n\x1a\n')
    print('Verified: 1 summary, 5 CSV tables, 2 FASTA files, 2 PNG figures (10/10).')


if __name__ == '__main__':
    main()
