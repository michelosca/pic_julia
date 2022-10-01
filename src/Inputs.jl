module Inputs

using SharedData: System, Species, Particle
using Constants: me, kb, e, amu
using Constants: gc
using Constants: c_bc_periodic, c_bc_open


function GetSystemParameters()
    system = System()

    system.ncells = 100
    system.x_min = 0.0
    system.x_max = 0.1
    system.Lx = system.x_max - system.x_min
    system.dx = system.Lx / (system.ncells-1)
    
    system.t_start = 0.0
    system.t_end = 1.e-6
    system.time = system.t_start
    system.dt = 1.e-8
    system.step = 0

    system.V0_min = 0.0
    system.V0_max = 100.0

    system.bc_field= c_bc_open
    system.bc_part= c_bc_open

    if system.bc_field == c_bc_periodic
        system.V0_min = 0.0
        system.V0_max = 0.0
    end

    return system
end

function GetSpeciesList(system::System)

    species_list = Species[]

    n_particles = 1000

    counter = 0
    # Add electrons
    counter += 1
    temp = 2.0 * e/kb
    is_background = false
    electrons = SetNewSpecies(counter, "electrons", me, -e, 1.e8, n_particles, temp, is_background, system)
    push!(species_list, electrons)

    # Add Ar ions
    counter += 1
    temp = 300.0
    is_background = false
    Ar_ions = SetNewSpecies(counter, "Ar+", 40*amu, e, 1.e8, n_particles, temp, is_background, system)
    push!(species_list, Ar_ions)

    # Add Ar background gas 
    counter += 1
    temp = 300.0
    is_background = true 
    Ar_gas = SetNewSpecies(counter, "Ar", 40*amu, e, 0.0, 0, temp, is_background, system)
    push!(species_list, Ar_gas)

    return species_list
end

function InitParticle(system::System, species::Species, temp::Float64)
    part = Particle()
    sigma = sqrt(temp * kb / species.mass)
    part.pos = rand() * system.Lx + system.x_min
    part.vel = randn(Float64, 3) * sigma
    return part
end


function SetNewSpecies(counter::Int64, name::String, mass::Float64,
    charge::Float64, part_weight::Float64, part_count::Int64, temp::Float64,
    background_species::Bool, system::System)

    species = InitSpeciesBlock(system)
    species.id = counter
    species.name = name
    species.mass = mass
    species.charge = charge 
    species.weight = part_weight
    species.particle_count = part_count
    species.is_background_gas = background_species

    if !background_species
        for i in range(1,species.particle_count,step=1)
            part = InitParticle(system, species, temp)
            push!(species.particle_list, part)
        end
        
        species.particle_grid_list = Vector{Particle}[]
        species.particle_grid_list = Vector{Particle}[]
        for i in range(1, system.ncells, step=1)
            push!(species.particle_grid_list, Particle[])
        end
    end
    return species
end

function InitSpeciesBlock(system::System)
    species = Species()
    species.id = 0
    species.name = ""
    species.mass = -1.0
    species.charge = 0.0
    species.weight = -1.0
    species.dens = zeros(system.ncells)
    species.charge_dens = zeros(system.ncells)
    species.particle_list = Particle[]
    species.particle_count = 0
    species.is_background_gas = false
    return species
end



end