# Changelog

## [Unreleased] - 2023-10-12

### Added
### Changed 
### Removed 
### Fixed


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
