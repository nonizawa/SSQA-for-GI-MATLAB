# SSQA_MATLAB (Stochastic simulated quantum annelinag for graph isomophism)

![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)

Under construction...


SSQA is a MATLAB project that utilizes Stochastic Simulated Quantum Annelinag (SSQA) for Graph Isomorphism (GI) tasks. The project aims to develop a robust algorithm to assess the quality, by evaluating their performance in solving Graph Isomorphism problems. The GI problem is a classical computational problem in which two given graphs, G and H, are determined to be isomorphic if there exists a bijective mapping between their vertex sets that preserves adjacency.

This repository contains an implementation of the SSQA algorithm to tackle the GI problem as described in the research paper ["Stochastic Quantum Monte Carlo Algorithm for Large-Scale Combinatorial Optimization Problems"](https://arxiv.org/abs/2302.12454).

## Installation

### Clone the Repository

```sh
git clone https://github.com/nonizawa/SSQA_MATLAB.git
cd SSQA_MATLAB
```

## Structure

- `SA.m`: This is the MATLAB code that runs the SA for GI algorithm.
- `SSA.m`: This is the MATLAB code that runs the SSA for GI algorithm.
- `SSQA.m`: This is the MATLAB code that runs the SSQA for GI algorithm.
- `GI_dataset.mat` and `GI_dataset2.mat`: This is dataset of GI problems where edges between nodes exist with 50\%.

## Run

There are three algorithms: traditional SA, stochastic simulated annealing (SSA), and SSQA.

### Parameters for SA

- `N`: Number of nodes in a GI problem

- `Mcycle`: Number of cycles for 1 trial

- `run_c`: Number of trials to evaluate the performance on average

- `T_ini`:  A pseudo temperature starts at T_ini

- `T_end`: A pseudo temperature ends at T_end

- `tau`:  A pseudo temperature is decreased every tau cycle


## Contact

For any questions, issues, or inquiries, feel free to create an issue in the repository or contact the repository owner [@nonizawa](https://github.com/nonizawa).

## Citation

If you use this code in your research, please cite the original paper:
```bibtex
@misc{onizawa2023stochastic,
      title={Stochastic Quantum Monte Carlo Algorithm for Large-Scale Combinatorial Optimization Problems}, 
      author={Naoya Onizawa and Ryoma Sasaki and Duckgyu Shin and Warren J. Gross and Takahiro Hanyu},
      year={2023},
      eprint={2302.12454},
      archivePrefix={arXiv},
      primaryClass={quant-ph}
}
```


## License

This project is licensed under the MIT License.
