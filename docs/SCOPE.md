# One-day scope — 2026-10-24

## Question and system
Are ligand-contacting residues more conserved than other resolved residues in vertebrate carbonic anhydrase II (CA2)? Use reviewed UniProt CA2 sequences (taxonomy 7742) and experimental human CA2 structure 3HS4, chain A. This is an ortholog-focused case study, not a survey of all carbonic anhydrase paralogs.

## Frozen analysis choices
- Download and archive the full UniProt query response and RCSB PDB structure with source URLs, UTC timestamps and SHA-256 hashes. Sort by accession; retain canonical amino-acid sequences of length 240–280; collapse identical sequences without discarding their source metadata.
- Select the AZM molecule closest to the structural zinc (minimum heavy-atom distance). Require a zinc contact within 3 Å; exclude surface AZM molecules. A pocket residue has any protein heavy atom within **5 Å** of that ligand. First model, chain A, selected alternate conformers only; waters, zinc and other heterogens are not protein residues.
- Include the observed structure sequence as a mapping-only MAFFT input. Map its non-gap positions to PDB residue identifiers, including insertion codes; do not count it in conservation estimates. Run conda-installed MAFFT with `--auto --thread 1`.
- For each mapped column, compute coverage, Shannon entropy in bits and conservation `1 - H/log2(20)`, excluding gaps/unknowns. Compare columns with coverage >=0.8. Equal weight per unique sequence; no pseudocounts.
- Primary statistic: mean pocket conservation minus mean non-pocket conservation among resolved chain-A residues. Two-sided label-permutation test, 10,000 draws, seed 29, fixed pocket size, plus-one p-value. Report descriptive 4/5/6 Å sensitivity, with 5 Å as the sole primary analysis.

## Numerical acceptance criteria
1. At least 6 unique real CA2 sequences from at least 5 source species and one ligand-bound experimental structure, with provenance for both downloads.
2. At least 10 eligible pocket residues and 100 eligible non-pocket residues at 5 Å; map 100% of resolved canonical residues exactly once.
3. One aligned FASTA, per-residue CSV, pocket CSV, summary JSON, permutation null table, and two PNG figures; 10,000 permutations and a finite p-value in (0,1]. Statistical significance is not an acceptance criterion.
4. At least 8 passing pytest cases, including analytical entropy, gap handling, geometry, mapping and permutation behavior. CI runs tests and the complete offline analysis on committed inputs.
5. English README reports actual results, methods, reproduction commands and explicit limitations. Commit, push to origin/main and read remote commit and GitHub Actions conclusion.

## Out of scope
De novo fpocket detection and cross-structure pocket similarity: not done; a known-ligand geometric pocket directly answers the narrower question. Phylogenetic correction, solvent-accessibility matched nulls, experimental validation and proteome-scale generalization: not done. No thesis ligand or homelab data.
