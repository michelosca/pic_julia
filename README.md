# PIC-Julia
Particle-In-Cell (PIC) model for low temperature plasmas

##  Run PIC-Julia
 1. Load 'PIC_Julia' module
 2. Run 'run_pic(input_file::String)' function, where 'input_file' is a string with the path to the input-deck file
 3. Output data in HDF5® format
 
## Input description
 - Input data is structured per 'blocks'
 - Each block starts with 'begin:<block_name>' and ends with 'end:<block_name>'
