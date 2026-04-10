This repository is currently under construction.

# Overview
This repository contains the input files required to reproduce the results presented in:
[![arXiv](https://img.shields.io/badge/arXiv-2510.05020-b31b1b.svg)](https://doi.org/10.48550/arXiv.2511.22504)
[![DOI](https://img.shields.io/badge/DOI-10.48550/arXiv.2510.05020-blue)](https://doi.org/10.48550/arXiv.2511.22504)

**Screening novel cathode materials from the Energy-GNoME database using MACE machine learning force field and DFT**· Nada Alghamdi, Paolo de Angelis, Pietro Asinari, Eliodoro Chiavazzo

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

If you use this data in your research, please cite the follwoing pre-print:

```bibtex
@article{alghamdi2025screening,
  title={Screening novel cathode materials from the Energy-GNoME database using MACE machine learning force field and DFT},
  author={Alghamdi, Nada and de Angelis, Paolo and Asinari, Pietro and Chiavazzo, Eliodoro},
  journal={arXiv preprint arXiv:2511.22504},
  year={2025}
}
```

If you use the Energy-GNoME data, please follow citaiton indeicate by [Energy-GNoME](https://paolodeangelis.github.io/Energy-GNoME/)
