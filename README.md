# FPBA: Flexible Percentile-Based Allocation

This repository contains the code and resources related to our research paper titled "FPBA: Flexible Percentile-Based Allocation for Multiple-Bits-Per-Cell RRAM", published at the ASP-DAC 2025.

## Run Instructions

The major test script is experiments/comprehensive_test.py. Please follow the comments in main to decide what tests are run. Results are saved under experiments/all_tests/. Results are generally saved in JSON format with each key being the method tested, values[0] being raw Bit-Error-Rate(BER), values[1] being ECC overhead, and values[2] being minimum gamma(if applicable).

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

We would like to thank the authors of the citation[29] in the FPBA paper ([PBA: Percentile-Based Level Allocation for Multiple-Bits-Per-Cell RRAM](https://github.com/Anjiang-Wei/PBA)) for making their code and results public, which greatly assisted our research. 
- Specific scripts/folders we used:
  - [relaxation model data](https://github.com/JeffLiu114514/FlexiblePBA/tree/main/model)
  - [bpc settings](https://github.com/JeffLiu114514/FlexiblePBA/tree/main/settings)
  - [dala.py](https://github.com/JeffLiu114514/FlexiblePBA/blob/main/experiments/dala.py), [dala_genmatrix.py](https://github.com/JeffLiu114514/FlexiblePBA/blob/main/experiments/dala_genmatrix.py), [ecc.py](https://github.com/JeffLiu114514/FlexiblePBA/blob/main/experiments/ecc.py), [trans.py](https://github.com/JeffLiu114514/FlexiblePBA/blob/main/experiments/trans.py) under experiments/

## Citations

```TeX
@inproceedings{10.1145/3658617.3697569,
author = {Liu, Junfei and Kahng, Anson},
title = {FPBA: Flexible Percentile-Based Allocation for Multiple-Bits-Per-Cell RRAM},
year = {2025},
isbn = {9798400706356},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
url = {https://doi.org/10.1145/3658617.3697569},
doi = {10.1145/3658617.3697569},
booktitle = {Proceedings of the 30th Asia and South Pacific Design Automation Conference},
pages = {1202–1208},
numpages = {7},
location = {Tokyo, Japan},
series = {ASPDAC '25}
}
```
