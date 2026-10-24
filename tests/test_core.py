import json
import hashlib
from pathlib import Path
import numpy as np
import pytest
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from pocket_conservation.core import conservation, minimum_distance, map_columns, permutation_test, structure_site

ROOT = Path(__file__).resolve().parents[1]


def test_invariant_column():
    assert conservation('AAAA') == (1, 0, 1)


def test_maximum_entropy():
    coverage, entropy, score = conservation('ACDEFGHIKLMNPQRSTVWY')
    assert coverage == 1
    assert entropy == pytest.approx(np.log2(20))
    assert score == pytest.approx(0)


def test_gaps_and_unknown_excluded():
    assert conservation('AA-X') == (0.5, 0, 1)
    assert np.isnan(conservation('---')[1])


def test_two_state_entropy():
    assert conservation('AACC')[1] == pytest.approx(1)


def test_heavy_atom_distance_boundary():
    assert minimum_distance(np.array([[0., 0, 0], [20, 0, 0]]), np.array([[3., 4, 0]])) == 5
    with pytest.raises(ValueError):
        minimum_distance([], [[0, 0, 0]])


@pytest.fixture
def site():
    model = PDBParser(QUIET=True).get_structure('CA2', ROOT / 'data/raw/3HS4.pdb')[0]
    return model, structure_site(model)


def test_active_site_selects_zinc_bound_ligand(site):
    model, (residues, ligand, distance, distances) = site
    assert len([r for r in model['A'] if r.resname == 'AZM']) == 3
    assert distance < 3
    assert sum(d <= 5 for d in distances) >= 10
    assert len(distances) == len(residues)


def test_mapping_does_not_assume_pdb_number_equals_sequence_position(site):
    _, (residues, _, _, _) = site
    sequence = ''.join(seq1(r.resname) for r in residues)
    aligned = '--' + sequence[:10] + '---' + sequence[10:]
    mapping = map_columns(aligned, residues)
    assert mapping[0][1] == 2
    assert mapping[10][1] == 15
    assert len(mapping) == len(residues)
    assert mapping[-1][0].id == residues[-1].id
    with pytest.raises(ValueError):
        map_columns('A', residues)


def test_permutation_constant_null():
    delta, p, null = permutation_test([1, 1, 1, 1], [True, True, False, False], draws=100)
    assert delta == 0 and p == 1
    np.testing.assert_array_equal(null, np.zeros(100))


def test_permutation_seed_and_plus_one():
    args = ([1, 1, 0, 0], [True, True, False, False])
    delta, p, null = permutation_test(*args, draws=1000)
    assert delta == 1
    assert p == (1 + sum(abs(null) >= 1)) / 1001
    np.testing.assert_array_equal(null, permutation_test(*args, draws=1000)[2])
    # Exact two-sided probability for this 2-of-4 case is 2/6.
    assert abs(p - 1/3) < 0.06


@pytest.mark.parametrize('mask', [[False, False], [True, True]])
def test_empty_groups_rejected(mask):
    with pytest.raises(ValueError):
        permutation_test([0.2, 0.8], mask)


def test_snapshots_match_provenance():
    manifest = json.loads((ROOT / 'data/raw/provenance.json').read_text())
    assert len(manifest) == 2
    for item in manifest:
        data = (ROOT / 'data/raw' / item['file']).read_bytes()
        assert len(data) == item['bytes']
        assert hashlib.sha256(data).hexdigest() == item['sha256']


def test_mapping_preserves_insertion_codes():
    from Bio.PDB.Residue import Residue
    residues = [Residue((' ', 10, ' '), 'ALA', ''), Residue((' ', 10, 'A'), 'CYS', '')]
    mapping = map_columns('A-C', residues)
    assert [(r.id, col) for r, col in mapping] == [((' ', 10, ' '), 0), ((' ', 10, 'A'), 2)]
