# Changelog

## [Unreleased] - 2023-10-12

### Added
- run_gsm.py: added test if gsm command is available
- filter_gsm.py: added step where optimize product with xTB and if RMSD change is big, filter out the reaction


### Added
### Changed 
### Removed 
### Fixed

## [0.7.2] - 2024-09-02
### Added
- Utility that handle db and h5 files are added to utils/db_h5_tools.
  1. db_to_h5.py
  2. h5_to_db.py
  3. make_db_from_xyzs.py
  4. merge_db.py
  5. merge_h5.py

## [0.7.1] - 2024-09-01

### Added
- Normal mode sampling codes are added to utils/nms.
  1. normal_mode_sampling.py
  2. refine_forces_nms.py

### Fixed
- All code now assert the type of the input_path (dir or file)

## [0.7.0] - 2024-08-31

### Added
- Sampling iso/conformers is included as a preparatory step in dandelion.
  1. smiles_to_isoconfs.py
  2. geom_opt.py
  3. dandelion_prep.py
   
- cli.py: to invoke dandelion like 'dand prep -i ./a.smi -n 40' in cli.

### Changed 
- dandelion is shortend as 'dand' in cli.
- dandelion_sample.py: default argument '0_mothers' changed to '0_reactants'
- print_separator, merge_args_with_defaults are moved to init.py
## [0.6.2] - 2024-07-08

### Fixed
- dandelion.py: name changed to dandelion_refine.py



## [0.6.1] - 2024-01-14

### Changed
- compile_refined.py: bug fixed when atomrow doesn't have 'energy' and 'forces'


## [0.6.0] - 2023-11-21

### Added
- filter_neb.py: added function is_valid_reaction to filter out weird rxn

## [0.5.6] - 2023-11-21

### Changed
- refine_forces.py: suppress error in force calculation, save to orca_error.log
- refine_forces.py: now save samples in batch
- refine_forces.py: open .db file with statement

## [0.5.5] - 2023-11-14

### Added
- run_neb.py: added argument fmax_threshold (default=0.1ev/A)

### Fixed
- refine_forces.py: added NoTrah for the orca command


## [0.5.4] - 2023-11-7

### Fixed
- compile_neb.py: fixed argparser that had no required=True


## [0.5.3] - 2023-11-2

### Fixed
- dandelion_refine.py: awesome ascii art


## [0.5.2] - 2023-11-2

### Fixed
- compile_refined.py: sorting the rows in the right order


## [0.5.1] - 2023-10-17

### Added
- opt_mothers.py: optimize crude structures using xTB


## [0.5.0] - 2023-10-12

### Added
- filter_neb.py: xTB normal mode TS validation: is_transition_state


## [0.4.1] - 2023-10-11

### Added
- Added \__init__.py have variable \__version__

### Fixed
- Basis set 6-31g(d) for Br atom in orca was handled thanks to https://github.com/ZimmermanGroup/ORCA-Basis-Sets


## [0.4.0] - 2023-10-10

### Added
- dandelion_refine.py that run refine processes

  
## [0.3.1] - 2023-10-10

### Added
- setup.py, README.md, CHANGELOG.md, LICENSE added 


## [0.2.0] - 2023-09-30

### Added
- dandelion.py that run through neb, refine
- Codes refactored

### Fixed
- Issues with absolute import fixed


## [0.1.0] - 2023-09-10

### Added
- Initial release with features neb, refine, segsm
