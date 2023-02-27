# ePhys


This is a project for opening electrophysiology/engineering data in Julia. 

This package will utilize many different sources for it's data. Each file format will gets its own package. 

### ROADMAP
Version 1.0: Can open, plot, and analyze datafiles from .ABF files. 

To Do list: 
- [ ] (v1.3.0) Allow for saving .abf files and modifying
- [ ] (v1.2.0) Open .mat files (For use with MatLab and Symphony)
- [ ] (v1.2.0) Open .idata files (For use with MatLab and IrisData https://github.com/sampath-lab-ucla/IrisDVA)
- [ ] (v1.1.1) Update some of the data analysis functions and expand analysis  
- [ ] (v1.1.1) Open .csv files (Some formats are saved as CSV files, especially from LabView products)

Completed Tasks: 
- [x] (v0.9.0) Make a Pluto.jl data entry suite to use as analysis GUI 
- [x] (< v1.0.0) Some basic datasheet manipulation
- [x] (< v1.0.0)Experiments can be plotted
- [x] (< v1.0.0) Experiment struct can be added, subtracted, multiplied, and divided
- [x] (< v1.0.0) Open .abf files (ABFReader.jl)