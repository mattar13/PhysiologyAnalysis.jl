# PhysiologyAnalysis.jl

[![License][license-img]](LICENSE)

[license-img]: http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat-square

This is a project for opening electrophysiology/engineering data in Julia. 

This package will utilize many different sources for it's data. Each file format will gets its own package. 

### ROADMAP
Version 1.0: Can open, plot, and analyze datafiles from .ABF files. 

Notebooks are included to help analyze some electroretinography data.
To install notebooks use this link: https://github.com/mattar13/PhysiologyInterface.jl


To Do list: 
- [ ] (v0.3.0) Allow for saving .abf files and modifying
- [ ] (v0.2.0) Open .mat files (For use with MatLab and Symphony)
- [ ] (v0.2.0) Open .idata files (For use with MatLab and IrisData https://github.com/sampath-lab-ucla/IrisDVA)
- [ ] (v0.1.0) Update some of the data analysis functions and expand analysis  
- [ ] (v0.1.0) Open .csv files (Some formats are saved as CSV files, especially from LabView products)

Completed Tasks: 
- [x] (v0.1.0) Make a Pluto.jl data entry suite to use as analysis GUI 
- [x] (< v0.1.0) Some basic datasheet manipulation
- [x] (< v0.1.0)Experiments can be plotted
- [x] (< v0.1.0) Experiment struct can be added, subtracted, multiplied, and divided
- [x] (< v0.1.0) Open .abf files (ABFReader.jl)
