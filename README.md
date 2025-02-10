# Dandelion 
Dandelion is a code for generating datasets that contain both equilibrium and reactive regions of potential energy surfaces, using automated and efficient sampling of chemical reaction space.

Documentation : <https://mhyeok1.github.io/dand_docs/>



## Citation
If you find this work useful for your research, please consider citing:
- Lee et al. *Adv. Sci.* 202409009 (2025) [LINK](https://doi.org/10.1002/advs.202409009)

This work builds upon two pioneering papers that should also be cited:
- Grambow et al. *Sci. Data* **7**, 137 (2020) [LINK](https://doi.org/10.1038/s41597-020-0460-4)
- Schreiner et al. *Sci. Data* **9**, 779 (2022) [LINK](https://doi.org/10.1038/s41597-022-01870-w)

## Supporting Information
<https://doi.org/10.5281/zenodo.14020916>

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


