# SSQA_MATLAB (Stochastic simulated quantum annelinag for graph isomophism)

![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)

SSQA is a MATLAB project that utilizes Stochastic Simulated Quantum Annelinag (SSQA) for Graph Isomorphism (GI) tasks. The project aims to develop a robust algorithm to assess the quality, by evaluating their performance in solving Graph Isomorphism problems. The GI problem is a classical computational problem in which two given graphs, G and H, are determined to be isomorphic if there exists a bijective mapping between their vertex sets that preserves adjacency.

This repository contains an implementation of the SSQA algorithm to tackle the GI problem as described in the research paper ["Stochastic Simualted Quantum Annealing for Fast Solution of Combinatorial Optimization Problems"](https://arxiv.org/abs/2302.12454).

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


### Parameters for SSA

- `N`: Number of nodes in a GI problem

- `Mcycle`: Number of cycles for 1 trial

- `run_c`: Number of trials to evaluate the performance on average

- `I0_min`:  A pseudo inverse temperature starts at I0_min

- `I0_max`: A pseudo inverse temperature ends at I0_max

- `tau` and `beta`:  A pseudo inverse temperature is increased by `beta` every tau cycle

- `nrnd`: A noise magnitude for each spin

### Parameters for SSQA

- `N`: Number of nodes in a GI problem

- `Mcycle`: Number of cycles for 1 trial

- `run_c`: Number of trials to evaluate the performance on average

- `I0`:  A pseudo inverse temperature for each spin

- `M`: Number of layers (number of replicas of spins)

- `T`: Timing of inter-layer interations

- `Q_max`: The maximum inter-layer interaction

- `tau` and `beta`:  An inter-layer interaction is increased by `beta` every tau cycle

- `NI`: Number of iterations. In each iteration, an inter-layer interaction starts from 0 to Q_max.

- `nrnd`: A noise magnitude for each spin

## Contact

For any questions, issues, or inquiries, feel free to create an issue in the repository or contact the repository owner [@nonizawa](https://github.com/nonizawa).

## Citation

If you use this code in your research, please cite the original paper:
```bibtex
@misc{onizawa2023stochastic,
      title={Stochastic Simulated Quantum Annealing for Fast Solution of Combinatorial Optimization Problems}, 
      author={Naoya Onizawa and Ryoma Sasaki and Duckgyu Shin and Warren J. Gross and Takahiro Hanyu},
      year={2023},
      eprint={2302.12454},
      archivePrefix={arXiv},
      primaryClass={quant-ph}
}
```


## License

This project is licensed under the MIT License.
