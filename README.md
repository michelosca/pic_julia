# **PIC-Julia**
**Particle-In-Cell (PIC)** (1D) model with **Monte Carlo collisions (MCC)** for low temperature plasmas. 
<!--
MCC null-collision method ([Vahedi *et al*, 1995](https://doi.org/10.1016/0010-4655(94)00171-W)) and  cross-section data is found in *collisions* folder
-->

## Content
<!-- TOC start -->
- [**Run PIC-Julia**](#run-pic-julia)
- [**Input description**](#input-description)
  * [**Input block: SYSTEM**](#input-block-system)
  * [**Input block: SPECIES**](#input-block-species)
<!-- TOC end -->

## **Run PIC-Julia**
 1. Include PIC-Julia folder to ``LOAD_PATH`` and load PIC_Julia module

    ```Julia
    push!(LOAD_PATH, /path/to/PIC-Julia/)

    using PIC_Julia
    ```

 2. Run ``run_pic(input_file::String)`` function, where 'input_file' is a string with the path to the input-deck file
    ```Julia
    run_pic("path/to/input.deck")
    ```
 3. Output data in HDF5® format
 4. Simulation log file is generated
 

## **Input description**
<div id='id-input_deck'/>
 - Input data is structured in *blocks*
 - Each block starts with ``begin:<block_name>`` and ends with ``end:<block_name>``
 - The following blocks are available
    - system: simulation parameters
    - species: charged and neutrals species to be simulated
    - output: defines physical parameters to be output 
    - MCC: Monte Carlo collisions
    - waveform: potential waveforms applied at the boundaries
- SI units apply for physical parameters

### **Input block: SYSTEM**
The system block gathers the simulation parameters required

- *dx* or *cell_width*: cell with
- *x_min*: left boundary location in space, default 0.0
- *x_max*: right boundary location in space, default 0.0
- *Lx*: simulation domain length, default 0.0. In case both x_min/x_max and Lx are defined, the latter prevails.
- *ncells* / *cells*: number of grid cells
- *integration_shape* / *interpolation_shape*: super-particle shape
    - Possible variables: *triangle*
- *particle_bc_min* / *part_bc_min*: particle boundary condition in the left boundary
    - Possible variables: *open*, *periodic*
- *particle_bc_max* / *part_bc_max*: particle boundary condition in the right boundary
    - Possible variables: *open*, *periodic*
- *particle_bc* / *part_bc*: sets particle boundary conditions in both left and right boundaries
    - Possible variables: *open*, *periodic*
- *field_bc_min*: field boundary condition in the left boundary
    - Possible variables: *open*, *periodic*
- *field_bc_max*: field boundary condition in the right boundary
    - Possible variables: *open*, *periodic*
- *field_bc*: field boundary conditions in both left and right boundaries
    - Possible variables: *open*, *periodic*
- *dt*: time step, default 0.0
- *t_end*: simulation time, default 0.0
- *step_end*: simulation steps, default 0

### **Input block: SPECIES**
A species block is required for each species to be simulated

- *name*: name of the species
- *charge*: species electric charge, in electron charge units
- *mass*: species particle mass
- *part_per_cell*: initial number of super-particle loaded per grid-cell
- *particles* / *total_particles*: total number of super-particles initially loaded
- *init_dens* / *density* / *dens*: initial number density
- *init_temp* / *T* / *temp*: initial temperature
- *spatial_dist* / *spatial_distribution*: spatial number density distribution function. This function is multiplied by the *density* parameter described above in this block. 
- *weight* / *part_weight* / *part_ratio*: super-particle to real particle ratio
- *background* / *is_background* / *background_species*: in case the species defined is just a static background field (usually for neutral gases)
    - Possible variables: ``true`` / ``false``