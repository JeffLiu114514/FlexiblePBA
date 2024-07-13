# FPBA: Flexible Percentile-Based Allocation

This repository contains the code and resources related to our research paper titled "FPBA: Flexible Percentile-Based Allocation for Multiple-Bits-Per-Cell RRAM", submitted as part of the ASP-DAC 2025 submission.

## Run Instructions

The major test script is experiments/comprehensive_test.py. Please follow the comments in main to decide what tests are run. Results are saved under experiments/all_tests/

## Directory Structure

The repository is organized as follows:\
├── experiments/ # Contains all the scripts used in the project\
│ ├── all_tests/ # Temporary folder for test results\
│ ├── coding/ # Contains scripts for gray coding\
│ ├── ecc/ # Necessary files for ECC overhead calculation\
│ ├── ember_capacity/ # Temporary files\
│ ├── comprehensive_test.py # Main test file. Tailor what to test in main.\
│ └── ...\
├── model/ # Contains the relaxation models. The experiments input.\
├── settings/ # Contains different bits-per-cell(bpc) settings.\
├── tests/ # Contains our test results archived.\
├── README.md # This README file\
└── FPBA-appendix-vF.pdf # The Appendix of our submitted paper.

## Acknowledgements

We would like to thank the authors of the citation[29] in the FPBA paper for making their code and results public, which greatly assisted our research. 
- Specific scripts/folders we used:
  - [relaxation model data](https://github.com/JeffLiu114514/FlexiblePBA/tree/main/model)
  - [bpc settings](https://github.com/JeffLiu114514/FlexiblePBA/tree/main/settings)
  - [dala.py](https://github.com/JeffLiu114514/FlexiblePBA/blob/main/experiments/dala.py), [dala_genmatrix.py](https://github.com/JeffLiu114514/FlexiblePBA/blob/main/experiments/dala_genmatrix.py), [ecc.py](https://github.com/JeffLiu114514/FlexiblePBA/blob/main/experiments/ecc.py), [trans.py](https://github.com/JeffLiu114514/FlexiblePBA/blob/main/experiments/trans.py) under experiments/
