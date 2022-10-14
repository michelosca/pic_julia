module Inputs

using SharedData: System, Species, Particle
using Constants: me, kb, e, amu
using Constants: gc_triangle
using Constants: c_bc_periodic, c_bc_open


function GetSystemParameters()
    system = System()

    system.ncells = 100
    system.gc = gc_triangle
    system.cell_min = system.gc + 1
    system.cell_max = system.gc + system.ncells
    system.ncells_total = system.ncells + 2*system.gc
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
    system.V0_max = 1.0

    system.bc_field= c_bc_periodic #open #
    system.bc_part=  c_bc_periodic #open #

    if system.bc_field == c_bc_periodic
        system.V0_min = 0.0
        system.V0_max = 0.0
    end

    return system
end

function GetSpeciesList(system::System)

    species_list = Species[]

    L = system.Lx

    n_particles = 100000

    counter = 0
    # Add electrons
    counter += 1
    temp = 2.0 * e/kb
    is_background = false
    electron_dist = x -> 1 + 0.5*sin(2*pi/L*x)
    electrons = SetNewSpecies(counter, "electrons", me, -e, 1.e8, n_particles,
        temp, is_background, system, electron_dist )
    push!(species_list, electrons)

    # Add Ar ions
    counter += 1
    temp = 300.0
    is_background = false
    Ar_ions = SetNewSpecies(counter, "Ar+", 40*amu, e, 1.e8, n_particles, temp,
        is_background, system, x -> 1)
    push!(species_list, Ar_ions)

    # Add Ar background gas 
    counter += 1
    temp = 300.0
    is_background = true 
    Ar_gas = SetNewSpecies(counter, "Ar", 40*amu, e, 0.0, 0, temp, is_background,
        system, x -> 1)
    push!(species_list, Ar_gas)

    return species_list
end

function InitParticle(species::Species, temp::Float64, part_pos::Float64)
    part = Particle()
    sigma = sqrt(temp * kb / species.mass)

    part.pos = part_pos
    part.vel = randn(Float64, 3) * sigma

    return part
end


function SetNewSpecies(counter::Int64, name::String, mass::Float64,
    charge::Float64, part_weight::Float64, part_count::Int64, temp::Float64,
    background_species::Bool, system::System, spatial_distribution)

    species = InitSpeciesBlock(system)
    species.id = counter
    species.name = name
    species.mass = mass
    species.charge = charge 
    species.weight = part_weight
    species.particle_count = part_count
    species.is_background_gas = background_species

    if !background_species

        # Spatial particle distribution
        x_min = system.x_min
        x_max = system.x_max
        dx = system.dx
        spatial_cdf = non_uniform_particle_distribution(spatial_distribution, system)
        x_grid = range(x_min, x_max, step=dx)

        for i in range(1,species.particle_count,step=1)

            # Set particle's position
            R = rand()
            part_pos = nothing
            for x in x_grid
                if R < spatial_cdf(x)
                    part_pos = x
                    break
                end
            end
            if part_pos === nothing
                print("***ERROR*** Particle has not been located\n")
            end

            part = InitParticle(species, temp, part_pos)
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
    species.dens = zeros(Float64, system.ncells_total)
    species.particle_list = Particle[]
    species.particle_count = 0
    species.is_background_gas = false
    return species
end

function non_uniform_particle_distribution(dist, system::System)

    # System parameters
    dx = system.dx
    x_min = system.x_min
    x_max = system.x_max

    # The simulation edges are trimmed because they count only half a cell each
    x_grid = range(x_min+0.5*dx, x_max-0.5*dx, step=dx)

    # Normalized spatial distribution
    area = sum( map(dist, x_grid) ) * dx
    dist_norm = x -> dist(x) / area * dx

    # Cumulative distribution function
    cdf = x -> sum( map(dist_norm, x_grid[1:round(Int64,(x-x_min)/dx)+1] ) )

    return cdf
end

end