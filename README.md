This repository is currently under construction.

# Overview
This repository contains the input files required to reproduce the results presented in:
[![arXiv](https://img.shields.io/badge/arXiv-2510.05020-b31b1b.svg)](https://doi.org/10.48550/arXiv.2511.22504)
[![DOI](https://img.shields.io/badge/DOI-10.1088%2F2515--7655%2Fae6c5a-blue)](https://doi.org/10.1088/2515-7639/aeae74)

**Screening novel cathode materials from the Energy-GNoME database using foundation machine learning force fields and DFT**· Nada Alghamdi, Paolo de Angelis, Pietro Asinari, Eliodoro Chiavazzo

<img width="672" height="551" alt="screening" src="https://github.com/user-attachments/assets/5d43fc14-0bd8-4b4e-ae55-179fe0348833" />

## 📝 Contents

### `0_validation_of_computational_method/`
Validates the computational workflow against known materials.

#### `phonons_database_MACE_test/`
Contains files for validating MACE predicted dynamical stability against phonon database [REF](https://mdr.nims.go.jp/collections/8g84ms862?locale=en) and [REF](https://doi.org/10.1038/s41524-025-01650-1).

#### `vol_profile_existing_cathodes/`
Voltage profile validation using existing experimentally characterized cathode materials to benchmark accuracy.

### `1_screening_energy-GNoME/`
Screening for novel cathodes from the [Energy-GNoME](https://paolodeangelis.github.io/Energy-GNoME/) database with AI-confidence > 90%.

#### `0_phonons/`
Phonon calculations to assess dynamical stability of candidate materials from Energy-GNoME.

#### `1_voltage_profile/`
Voltage profile calculations to evaluate electrochemical performance of candidates.

### `2_dft_refinement/`
DFT refinement inputs for the most promising cathode candidates identified through screening.

## Citation

If you use this work in your research, please cite the follwoing:

```bibtex
@article{Alghamdi_2026,
doi = {10.1088/2515-7639/aeae74},
url = {https://doi.org/10.1088/2515-7639/aeae74},
year = {2026},
month = {oct},
publisher = {IOP Publishing},
volume = {9},
number = {4},
pages = {045007},
author = {Alghamdi, Nada and de Angelis, Paolo and Asinari, Pietro and Chiavazzo, Eliodoro},
title = {Screening novel cathode materials from the Energy-GNoME database using foundation machine learning force fields and DFT},
journal = {Journal of Physics: Materials},
}
```

If you use the Energy-GNoME data, please follow citaiton indeicate by [Energy-GNoME](https://paolodeangelis.github.io/Energy-GNoME/)
