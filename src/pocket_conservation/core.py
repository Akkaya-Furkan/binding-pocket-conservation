"""Explicit geometry, alignment mapping and conservation statistics."""
import numpy as np
from Bio.SeqUtils import seq1

AA = 'ACDEFGHIKLMNPQRSTVWY'


def conservation(column):
    valid = [a for a in column if a in AA]
    coverage = len(valid) / len(column) if len(column) else 0.0
    if not valid:
        return coverage, float('nan'), float('nan')
    counts = np.array([valid.count(a) for a in AA], dtype=float)
    p = counts[counts > 0] / len(valid)
    entropy = float(-np.sum(p * np.log2(p)))
    return coverage, entropy, 1 - entropy / np.log2(20)


def heavy_atoms(residue):
    return np.array([a.coord for a in residue.get_atoms() if a.element not in ('H', 'D')])


def minimum_distance(a, b):
    if len(a) == 0 or len(b) == 0:
        raise ValueError('Cannot measure an empty atom set')
    return float(np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(axis=2)).min())


def structure_site(model, chain_id='A'):
    chain = model[chain_id]
    residues = [r for r in chain if r.id[0] == ' ' and seq1(r.resname) in AA]
    ligands = [r for r in chain if r.resname == 'AZM']
    zinc = [r for r in chain if r.resname == 'ZN']
    if not residues or not ligands or not zinc:
        raise ValueError('Missing protein, AZM or zinc')
    metal = np.concatenate([heavy_atoms(r) for r in zinc])
    ligand = min(ligands, key=lambda r: minimum_distance(heavy_atoms(r), metal))
    zinc_distance = minimum_distance(heavy_atoms(ligand), metal)
    if zinc_distance > 3:
        raise ValueError('No zinc-bound AZM within 3 Å')
    distances = [minimum_distance(heavy_atoms(r), heavy_atoms(ligand)) for r in residues]
    return residues, ligand, zinc_distance, distances


def map_columns(aligned_structure, residues):
    expected = ''.join(seq1(r.resname) for r in residues)
    if aligned_structure.replace('-', '') != expected:
        raise ValueError('Structure sequence and alignment do not match')
    columns = [i for i, a in enumerate(aligned_structure) if a != '-']
    return list(zip(residues, columns))


def permutation_test(scores, pocket, draws=10000, seed=29):
    scores = np.asarray(scores, dtype=float)
    pocket = np.asarray(pocket, dtype=bool)
    if len(scores) != len(pocket) or not np.isfinite(scores).all():
        raise ValueError('Invalid scores or mask')
    n = int(pocket.sum())
    if n == 0 or n == len(scores) or draws < 1:
        raise ValueError('Need two nonempty groups and positive draws')
    observed = float(scores[pocket].mean() - scores[~pocket].mean())
    rng = np.random.default_rng(seed)
    null = np.empty(draws)
    for i in range(draws):
        chosen = rng.choice(len(scores), n, replace=False)
        total = scores[chosen].sum()
        null[i] = total / n - (scores.sum() - total) / (len(scores) - n)
    p = (1 + np.count_nonzero(np.abs(null) >= abs(observed) - 1e-12)) / (draws + 1)
    return observed, float(p), null
