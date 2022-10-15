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
        species.particle_list = Particle[]

        # Initial condition parameters
        species.init_dens = 0.0
        species.init_temp = 0.0
        species.is_background_species = false
        species.init_spatial_distribution = 1
        push!(species_list, species)
    end
    return errcode
end

function ReadSpeciesEntry!(name::SubString{String}, var::SubString{String},
    read_step::Int64, species_list::Vector{Species})

    errcode = c_error 

    if (read_step == 1)
        units, name = GetUnits!(name)
        species = species_list[end]

        if (name=="name")
            species.name = var
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
            species.init_spatial_distribution = 1
            return 0
        end

        if (name=="init_dens" || name=="density" || name=="dens")
            species.init_dens = parse(Float64, var) * units 
            return 0
        end

        if (name=="init_temp" || name=="T" || name=="temp")
            species.init_temp = parse(Float64, var) * units 
            return 0
        end

        if (name=="spatial_dist" || name=="spatial_distribution")
            species.init_spatial_distribution = Meta.parse(var)
            return 0
        end

        if (name=="weight" || name=="part_weight" || name=="part_ration")
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
    end

    return errcode 
end

function EndFile_Species!(read_step::Int64, species_list::Vector{Species},
    system::System)

    errcode = 0
    
    if (read_step == 1)

        # Species check
        for species in species_list
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

        # Load particles after all the checks required were successfull 
        for species in species_list
            errcode = LoadParticles!(species, system)
            if errcode == c_error
                err_message = @sprintf(" While loading %s", species.name)
                PrintErrorMessage(system, err_message)
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
    dist = species.init_spatial_distribution

    # Particle loading only for species that are not a background species
    if !species.is_background_species

        # Spatial grid
        x_min = system.x_min
        x_max = system.x_max
        dx = system.dx
        ncells = system.ncells
        x_grid = range(x_min, x_max, step=dx)

        # Spatial distribution function: normalization factor
        norm_factor = GetDistributionNormalizationFactor(dist, system)

        # Load particles
        for i in range(1, species.particle_count, step=1)

            # Set particle's position
            R = rand()
            part_pos = nothing
            x_0 = x_grid[1]
            for j in range(2, ncells, step=1)
                x = x_grid[j]
                if R < GetCDF_at_x(dist, norm_factor, x_grid[1:j], system)
                    part_pos = rand() * dx + x_0
                    break
                end
                x_0 = x
            end

            if part_pos === nothing
                PrintErrorMessage(system, "Particle has not been located in space")
                return c_error
            elseif part_pos > x_max
                PrintErrorMessage("Particle loaded at position > x_max")
                return c_error
            elseif part_pos < x_min
                PrintErrorMessage("Particle loaded at position < x_min")
                return c_error
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

    return errcode
end

function GetCDF_at_x(dist::Union{Int64, Float64, Expr}, norm::Float64,
    x_grid_range::StepRangeLen{Float64}, system::System)

    s = 0.0
    for x in x_grid_range[1:end-1]
        s += ReplaceExpressionValues(dist, system, x)
    end
    cdf_x = s/norm

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
        s += ReplaceExpressionValues(dist, system, x)
    end

    return s
end

function InitParticle(species::Species, temp::Float64, part_pos::Float64)
    part = Particle()
    sigma = sqrt(temp * kb / species.mass)

    part.pos = part_pos
    part.vel = randn(Float64, 3) * sigma

    return part
end

end