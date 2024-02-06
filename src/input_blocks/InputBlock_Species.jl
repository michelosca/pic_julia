# Copyright (C) 2021 Michel Osca Engelbrecht
#
# This file is part of GM Julia.
#
# GM Julia is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GM Julia is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GM Julia. If not, see <https://www.gnu.org/licenses/>.

module InputBlock_Species

using Constants: K_to_eV, e, me, amu, kb 
using Constants: c_error
using SharedData: Species, System, Particle
using Tools: GetUnits!
using Printf: @sprintf
using EvaluateExpressions: ReplaceExpressionValues
using PrintModule: PrintErrorMessage, PrintWarningMessage
using DataStructures: MutableLinkedList


function StartFile_Species!(read_step::Int64, species_list::Vector{Species},
    system::System)

    errcode = 0

    return errcode
end

function StartSpeciesBlock!(read_step::Int64, species_list::Vector{Species})

    errcode = 0
    if (read_step == 1)
        species = Species()
        species.id = length(species_list)+1 
        species.name = "None"
        species.mass = 0.0
        species.charge = 0.0
        species.weight = 0.0

        species.particle_count = 0
        species.part_per_cell = 0
        species.particle_list = MutableLinkedList{Particle}()

        # Initial condition parameters
        species.init_dens = 0.0
        species.init_temp = 0.0
        species.is_background_species = false
        species.dens_spatial_distribution = 1
        push!(species_list, species)
    end
    return errcode
end

function ReadSpeciesEntry!(name::SubString{String}, var::SubString{String},
    read_step::Int64, species_list::Vector{Species}, system::System)

    errcode = c_error 

    if (read_step == 1)
        units, name = GetUnits!(name)
        species = species_list[end]

        if (name=="name")
            species.name = strip(var)
            return 0
        end

        if (name=="charge")
            species.charge = parse(Int64,var) * e * units
            return 0
        end

        if (name=="mass")
            expr = Meta.parse(var)
            species.mass = eval(expr) * units
            return 0
        end

        if (name=="part_per_cell")
            expr = Meta.parse(var)
            species.part_per_cell = Int64(eval(expr))
            return 0
        end

        if (name=="particles" || name=="total_particles")
            expr = Meta.parse(var)
            species.particle_count = Int64(eval(expr))
            species.dens_spatial_distribution = 1
            return 0
        end

        if (name=="init_dens" || name=="density" || name=="dens")
            species.init_dens = parse(Float64, var) * units 
            return 0
        end

        if (name=="init_temp" || name=="T" || name=="temp")
            if units == e
                message = @sprintf("Be careful! Temperature in eV must be defined like 'temp_TeV' instead of 'temp_eV' as the latter refers to energy units")
                PrintWarningMessage(system, message)
            end
            species.init_temp = parse(Float64, var) * units 
            return 0
        end

        if (name=="spatial_dist" || name=="spatial_distribution")
            species.dens_spatial_distribution = Meta.parse(var)
            return 0
        end

        if (name=="weight" || name=="part_weight" || name=="part_ratio")
            species.weight = parse(Float64, var) * units 
            return 0
        end

        if (name=="background" || name=="is_background" ||
            name=="background_species")
            species.is_background_species = parse(Bool, var)
            return 0
        end
    elseif read_step == 2
        errcode = 0
    end

    return errcode 
end

function EndSpeciesBlock!(read_step::Int64, species_list::Vector{Species},
    system::System)

    errcode = 0 

    if read_step == 1
        species = species_list[end]
        if species.mass <= 0
            err_message = @sprintf("%s mass must be positive", species.name)
            PrintErrorMessage(system, err_message)
            return c_error
        end

        if species.init_temp < 0
            err_message = @sprintf("%s temperature must be positive", species.name)
            PrintErrorMessage(system, err_message)
            return c_error
        end

        if !species.is_background_species
            # Normalize spatial distribution function
            dist = species.dens_spatial_distribution
            norm_factor = GetDistributionNormalizationFactor(dist, system)
            norm_expr = Expr(:call,:/,dist,norm_factor)
            species.dens_spatial_distribution = norm_expr
        end
    end

    return errcode 
end

function EndFile_Species!(read_step::Int64, species_list::Vector{Species},
    system::System)

    errcode = 0
    
    if (read_step == 1)

        # Species check
        for species in species_list
            if !species.is_background_species
                if (species.part_per_cell > 0) && (species.particle_count == 0)
                    species.particle_count = species.part_per_cell * (system.ncells-1)
                end

                if (species.part_per_cell == 0) && (species.particle_count == 0)
                    err_message = @sprintf("Number of particles for %s has not been set", species.name)
                    PrintErrorMessage(system, err_message)
                    return c_error
                end

                if species.init_dens > 0 && species.particle_count > 0
                    real_part = system.Lx * species.init_dens
                    sim_part = Float64(species.particle_count)
                    species.weight = real_part / sim_part
                elseif species.init_dens > 0 && species.weight > 0
                    real_part = system.Lx * species.init_dens
                    species.particle_count = round(Int64, real_part / species.weight)
                elseif species.weight > 0 && species.particle_count > 0
                    real_part = species.particle_count * species.weight
                    species.init_dens = real_part / system.Lx 
                else
                    err_message = @sprintf("Missing declaration of particle count, density and/or particle weight for %s", species.name)
                    PrintErrorMessage(system, err_message)
                    return c_error
                end
            end
        end

        # Load particles after all the checks required were successfull 
        for species in species_list
            errcode = LoadParticles!(species, system)
            if errcode == c_error
                err_message = @sprintf(" While loading %s", species.name)
                PrintErrorMessage(system, err_message)
            end
        end
    else
        # Load secondary lists
        if system.mcc
            for species in species_list
                species.particle_grid_list = MutableLinkedList{Particle}[]
                # Note that particle_grid_list is ncells-1 long!
                # (because cells surrounding the boundaries are not necessary)
                for i in range(1, system.ncells-1, step=1)
                    push!(species.particle_grid_list, MutableLinkedList{Particle}())
                end
            end
        end
    end
    return errcode
end

function LoadParticles!(species::Species, system::System)

    errcode = 0

    # Particle temperature
    temp = species.init_temp

    # Load particle spatial distribution function
    dist = species.dens_spatial_distribution

    # Particle loading only for species that are not a background species
    if !species.is_background_species

        # Spatial grid
        x_min = system.x_min
        x_max = system.x_max
        dx = system.dx
        ncells = system.ncells
        x_grid = range(x_min, x_max, step=dx)
        
        # Total number of particles to be allocated
        part_count = species.particle_count

        # Cumulative distribution function at x_1 = x_min
        cdf_x_1 = GetCDF_at_x(dist, x_grid[1:1], system)
        x_1 = x_grid[1]

        # Loop over following grid cells and calculate the particles allocated on each cell
        species.particle_count = 0
        for i in range(2,ncells,step=1)
            # Get CDF on i-th position
            cdf_x = GetCDF_at_x(dist, x_grid[1:i], system)

            # Particle rate is the difference between i and i-1 CDF
            rate = cdf_x - cdf_x_1

            # Particles to be allocated
            i_part = round(Int64, part_count * rate)

            # Load particles
            for j in range(1, i_part, step=1)
                part_pos = rand() * dx + x_1
                AddNewParticle!(species, temp, part_pos)
            end
            
            # Set CDF at i-1 for next iteration
            cdf_x_1 = cdf_x
            x_1 = x_grid[i]
        end
        
    end

    return errcode
end

function GetCDF_at_x(dist::Union{Int64, Float64, Expr},
    x_grid_range::StepRangeLen{Float64}, system::System)

    cdf_x = 0.0
    for x in x_grid_range[1:end-1]
        cdf_x += ReplaceExpressionValues(dist, system, pos=x)
    end

    return cdf_x
end

function GetDistributionNormalizationFactor(dist::Union{Int64, Float64, Expr}, system::System)
    # Calculates the area under the distribution function
    # -> Note that to actually get the area the returned value should be
    #    multiplied by the cell width (dx)

    # System parameters
    dx = system.dx
    x_min = system.x_min+dx*0.5
    x_max = system.x_max-dx*0.5
    # The simulation edges are trimmed because they count only half a cell each
    s = 0.0
    for x in range(x_min, x_max, step=dx)
        s += ReplaceExpressionValues(dist, system, pos=x)
    end

    return s
end

function AddNewParticle!(species::Species, temp::Float64, part_pos::Float64)

    part = Particle()
    sigma = sqrt(temp * kb / species.mass)

    part.pos = part_pos
    part.vel = randn(Float64, 3) * sigma

    append!(species.particle_list, part)
    species.particle_count = species.particle_list.len

    return 0 
end

end