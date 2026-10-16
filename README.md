# One generation of the Standard Model as a quantum error-correcting code

Scripts accompanying the paper by Shiroman Prakash and K. L. Ramya, https://arxiv.org/abs/2610.21879.


## Python

Requires Python 3.8 or later, `numpy` and `sympy`. The full run takes less than an hour.

```
cd python
pip install -r requirements.txt
./run_all.sh
```

Each script stops with an `AssertionError` if a check fails. `qutrit_triorth_search.py` rewrites
`search_output.txt`, its log.

| Script | Topic |
|---|---|
| `letter_checks.py` | The 16 and its anomalies, field labels and charges, the printed generator matrix, the code for odd p and its gate, the choice of removed fermion, the qubit case, the single-U(1) example |
| `doubling_checks.py` | Doubling and Dirac pairing |
| `interactions_15.py` | The minimum-weight logicals as four-fermion relations |
| `removals_16.py` | Removing several fermions |
| `spinor_checks.py` | The spinor of so(14) and the 27 of E6 |
| `distillation.py` | Distillation figures |
| `automorphisms_15.py`, `automorphisms_p.py` | Automorphism groups (p = 3; p = 5, 7, 11) |
| `qutrit_triorth_search.py` | Exhaustive search for n ≤ 14, with self-tests |
| `uniqueness_15.py` | Search at n = 15 and uniqueness |

`sm16.py` contains shared definitions.

## Magma

The `magma` folder contains scripts written in Magma, which requires a licence. It is optional: the Python scripts above are enough to reproduce every result. See [`magma/README.md`](magma/README.md).

## Citation

```bibtex
@misc{prakash2026generationstandardmodelquantum,
      title={One generation of the Standard Model as a quantum error-correcting code}, 
      author={Shiroman Prakash and K. L. Ramya},
      year={2026},
      eprint={2610.21879},
      archivePrefix={arXiv},
      primaryClass={hep-th},
      url={https://arxiv.org/abs/2610.21879}, 
}
```

## License

[MIT](LICENSE).
