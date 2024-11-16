# Dandelion 
Dandelion is a code for generating datasets that contain both equilibrium and reactive regions of potential energy surfaces, using automated and efficient sampling of chemical reaction space.

Documentation : <https://jjy1031.github.io>

## References 
- Grambow, Colin A. et al., Scientific Data, 7.1, 137 (2020)
- Schreiner, Mathias et al., Scientific Data, 9.1, 779 (2022)

## Installation
- Prerequisites
  - conda
  - openmpi


### Download Dandelion

```python
git clone https://github.com/mhyeok1/dand.git
cd dand
```

### Setup conda environment

This creates a new Conda environment using the specifications provided in `environment.yml` file.

```python
conda env create -f environment.yml
conda activate ts
pip install -e .
```

### Install pyGSM

Visit the [pyGSM repository](https://github.com/ZimmermanGroup/pyGSM) for instructions or simply use:

```python
git clone https://github.com/ZimmermanGroup/pyGSM
pip install -e .
```
By executing `gsm` in terminal, you can verify that the program has been successfully installed.

### Install Orca

You can install ORCA [here](https://orcaforum.kofo.mpg.de/app.php/portal).
This extracts a tar.xz file.

```python
tar -xf orca.tar.xz
```
By executing `orca` in terminal, you can verify that the program has been successfully downloaded.

### Setup environment variables in `.bashrc`

You can open your `.bashrc` file using:
```python
vi ~/.bashrc
```

Add the following lines to your `.bashrc` file:

1. Add ORCA to your path :
```python
export PATH="/path/to/your/orca/directory:$PATH" \
```

2. Set the `PYTHONPATH` for pyGSM:
```python
export PYTHONPATH=/path/to/your/pyGSM/directory:$PYTHONPATH \
```

3. Adjust the `OMP_NUM_THREADS`:
```python
export OMP_NUM_THREADS=1
```

4. Set the xtb configuration:
```python
export OMP_STACKSIZE=16G \
ulimit -s unlimited\
```

Apply the changes:
```python
source ~/.bashrc
``` 


## Citation
If you are using Dandelion in your research, please cite:
