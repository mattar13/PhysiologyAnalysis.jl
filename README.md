# ePhys


This is a project for opening electrophysiology/engineering data in Julia. 

This package will utilize many different sources for it's data. Each file format will gets its own package. 

To Do list: 
- [ ] Make a Pluto.jl data entry suite to use as analysis GUI
- [ ] Open .mat files (For use with MatLab and Symphony)
- [ ] Open .idata files (For use with MatLab and IrisData https://github.com/sampath-lab-ucla/IrisDVA)
- [ ] Open .csv files (Some formats are saved as CSV files, especially from LabView products)
- [ ] Allow for saving .abf files and modifying

Completed Tasks: 
- [x] Experiment struct can be added, subtracted, multiplied, and divided
- [x] Open .abf files (ABFReader.jl)
- [x] Some basic datasheet manipulation
- [x] Experiments can eb plotted

