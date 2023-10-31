# Dandelion 

An efficient full-auto exploration of chemical compound space near the transition state utilizing tight binding, single-ended growing string method, and nudged elastic band.

## Contact

**Minhyeok Lee**  
Email: [mlee@yonsei.ac.kr](mailto:mlee@yonsei.ac.kr)  
Research Group: [TCCL Yonsei](http://tccl.yonsei.ac.kr/)

## Base References

- Grambow, Colin A. et al., _Scientific Data_, 7.1, 137 (2020)
- Schreiner, Mathias et al., _Scientific Data_, 9.1, 779 (2022)


## Setup

### 1. Conda Setup
```bash
conda env create -f environment.yml
conda activate ts
```
### 2. Install pyGSM
Visit the [pyGSM repository](https://github.com/ZimmermanGroup/pyGSM) for instructions or simply use:
```bash
pip install -e git+https://github.com/ZimmermanGroup/pyGSM.git#egg=pyGSM
```

### 3. Install Dandelion
```bash
cd Dandelion
pip install -e .
```
### 4. Install ORCA
Visit the [ORCA forum](https://orcaforum.kofo.mpg.de/) for installation instructions.  

### 5. Setup Environment Variables in `.bashrc`

Add ORCA to your path:
```bash
export PATH="/path/to/your/orca_directory:$PATH"
```

Set the `PYTHONPATH` for pyGSM:
```bash
export PYTHONPATH="/path/to/your/pyGSM_directory:$PYTHONPATH"
```

Adjust the `OMP_NUM_THREADS`:
```bash
export OMP_NUM_THREADS=1
```

Finally, set the xtb configuration:
```bash
export OMP_STACKSIZE=16G
ulimit -s unlimited
```