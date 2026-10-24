"""Run the entire offline analysis: python -m pocket_conservation."""
import csv
import hashlib
import json
import platform
import re
import subprocess
from pathlib import Path

import Bio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from .core import AA, conservation, structure_site, map_columns, permutation_test

ROOT = Path(__file__).resolve().parents[2]


def write_csv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    raw, out = ROOT / 'data/raw', ROOT / 'results'
    out.mkdir(exist_ok=True)
    provenance = json.loads((raw / 'provenance.json').read_text())
    for item in provenance:
        if hashlib.sha256((raw / item['file']).read_bytes()).hexdigest() != item['sha256']:
            raise ValueError('Input checksum mismatch: ' + item['file'])
    records = sorted(SeqIO.parse(raw / 'ca2_reviewed.fasta', 'fasta'), key=lambda r: r.id)
    unique, metadata = {}, []
    for r in records:
        sequence = str(r.seq)
        accession = r.id.split('|')[1]
        species = re.search(r'OS=(.*?) OX=', r.description).group(1)
        keep = 240 <= len(sequence) <= 280 and set(sequence) <= set(AA)
        representative = unique.setdefault(sequence, accession) if keep else ''
        metadata.append(dict(accession=accession, species=species, length=len(sequence),
                             retained=keep, representative=representative))
    write_csv(out / 'sequences.csv', metadata)
    model = PDBParser(QUIET=True).get_structure('3HS4', raw / '3HS4.pdb')[0]
    residues, ligand, zn_distance, distances = structure_site(model)
    structure_seq = ''.join(seq1(r.resname) for r in residues)
    with (out / 'input.fasta').open('w') as handle:
        for sequence, accession in unique.items():
            handle.write(f'>{accession}\n{sequence}\n')
        handle.write(f'>structure_3HS4_A\n{structure_seq}\n')
    command = ['mafft', '--auto', '--thread', '1', str(out / 'input.fasta')]
    with (out / 'alignment.fasta').open('w') as handle, (out / 'mafft.log').open('w') as log:
        subprocess.run(command, stdout=handle, stderr=log, check=True)
    alignment = {r.id: str(r.seq) for r in SeqIO.parse(out / 'alignment.fasta', 'fasta')}
    assert len(alignment) == len(unique) + 1
    family = [alignment[a] for a in unique.values()]
    rows = []
    for (r, col), distance in zip(map_columns(alignment['structure_3HS4_A'], residues), distances):
        coverage, entropy, score = conservation([s[col] for s in family])
        rows.append(dict(chain='A', pdb_residue=r.id[1], insertion_code=r.id[2].strip(),
                         amino_acid=seq1(r.resname), alignment_column=col + 1,
                         ligand_distance_A=distance, pocket=distance <= 5,
                         coverage=coverage, entropy_bits=entropy, conservation=score,
                         eligible=coverage >= 0.8))
    write_csv(out / 'residues.csv', rows)
    write_csv(out / 'pocket.csv', [r for r in rows if r['pocket']])
    eligible = [r for r in rows if r['eligible']]
    scores = np.array([r['conservation'] for r in eligible])
    mask = np.array([r['pocket'] for r in eligible])
    delta, p, null = permutation_test(scores, mask)
    write_csv(out / 'permutation_null.csv', [dict(iteration=i+1, mean_difference=v) for i, v in enumerate(null)])
    sensitivity = []
    for cutoff in (4, 5, 6):
        m = np.array([r['ligand_distance_A'] <= cutoff for r in eligible])
        sensitivity.append(dict(cutoff_A=cutoff, pocket_n=int(m.sum()), nonpocket_n=int((~m).sum()),
                                pocket_mean=float(scores[m].mean()), nonpocket_mean=float(scores[~m].mean()),
                                difference=float(scores[m].mean()-scores[~m].mean())))
    write_csv(out / 'sensitivity.csv', sensitivity)
    summary = dict(source_sequences=len(records), unique_sequences=len(unique),
                   source_species=len({r['species'] for r in metadata if r['retained']}),
                   resolved_residues=len(rows), mapped_residues=len(rows),
                   eligible_residues=len(eligible), pocket_n=int(mask.sum()), nonpocket_n=int((~mask).sum()),
                   pocket_mean=float(scores[mask].mean()), nonpocket_mean=float(scores[~mask].mean()),
                   mean_difference=delta, permutation_p_two_sided=p, permutations=len(null), seed=29,
                   cutoff_A=5, coverage_threshold=0.8,
                   ligand=dict(name=ligand.resname, chain='A', residue=ligand.id[1],
                               zinc_distance_A=zn_distance),
                   versions=dict(python=platform.python_version(), biopython=Bio.__version__,
                                 numpy=np.__version__, matplotlib=matplotlib.__version__,
                                 mafft=subprocess.run(['mafft', '--version'], capture_output=True, text=True).stderr.strip()),
                   mafft_command=' '.join(command[:-1] + ['results/input.fasta']))
    assert len(unique) >= 6 and summary['source_species'] >= 5
    assert summary['pocket_n'] >= 10 and summary['nonpocket_n'] >= 100
    assert len({(r['pdb_residue'], r['insertion_code']) for r in rows}) == len(residues)
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
    axes[0].boxplot([scores[mask], scores[~mask]], tick_labels=['Pocket', 'Non-pocket'], showmeans=True)
    axes[0].set(ylabel='Conservation (1 − H / log₂20)', ylim=(-0.02, 1.02))
    axes[1].hist(null, bins=45, color='#547b9e')
    axes[1].axvline(delta, color='#b33b36', label=f'Observed Δ = {delta:.4f}')
    axes[1].set(xlabel='Permuted mean difference', ylabel='Draws', title=f'Two-sided p = {p:.4g}')
    axes[1].legend()
    fig.savefig(out / 'comparison.png', dpi=180)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(11, 3.6), layout='constrained')
    ax.plot([r['pdb_residue'] for r in eligible], scores, color='#637689', lw=1)
    ax.scatter([r['pdb_residue'] for r in eligible if r['pocket']], scores[mask], color='#b33b36', label='5 Å pocket', zorder=3)
    ax.set(xlabel='3HS4 chain A residue number', ylabel='Conservation', ylim=(-0.02, 1.05))
    ax.legend()
    fig.savefig(out / 'conservation_profile.png', dpi=180)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
