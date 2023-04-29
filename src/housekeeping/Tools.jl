module Tools

using SharedData: Species, System, Field
using Constants: c_error
using Constants: c_bc_open, c_bc_periodic
using Constants: c_stag_centre, c_stag_right
using Constants: c_field_magnetic, c_field_electric
using Constants: c_bc_x_min, c_bc_x_max
using Constants: K_to_eV, e
using PrintModule: PrintWarningMessage

function InterpolateParticleToGrid!(field::Vector{Float64}, part_pos::Float64, system::System)

    cell_x, gx = ParticleToGrid(part_pos, system, c_stag_centre)

    # Add particle values to field grid
    field[cell_x - 1] += gx[1]
    field[cell_x]     += gx[2]
    field[cell_x + 1] += gx[3]
end

function InterpolateFieldToPoint(field::Field, part_pos::Float64, system::System)

    stag_list = zeros(Int64, 3)

    # Electric field
    # - Ex is centre staggered
    # - Ey and Ez are right staggered
    # Magnetic field
    # - Bx is right staggered
    # - By and Bz are centre staggered

    if field.id == c_field_electric
        stag_list[1] = c_stag_centre
        stag_list[2] = c_stag_right
        stag_list[3] = c_stag_right
    elseif field.id == c_field_magnetic
        stag_list[1] = c_stag_right
        stag_list[2] = c_stag_centre
        stag_list[3] = c_stag_centre
    end

    cell_x, gx = ParticleToGrid(part_pos, system, stag_list[1])
    cell_y, gy = ParticleToGrid(part_pos, system, stag_list[2])
    cell_z, gz = ParticleToGrid(part_pos, system, stag_list[3])

    field_at_point = zeros(Float64, 3) 
    field_at_point[1] += gx[1] * field.x[cell_x - 1]
    field_at_point[1] += gx[2] * field.x[cell_x]
    field_at_point[1] += gx[3] * field.x[cell_x + 1]

    field_at_point[2] += gy[1] * field.y[cell_y - 1]
    field_at_point[2] += gy[2] * field.y[cell_y]
    field_at_point[2] += gy[3] * field.y[cell_y + 1]

    field_at_point[3] += gz[1] * field.z[cell_z - 1]
    field_at_point[3] += gz[2] * field.z[cell_z]
    field_at_point[3] += gz[3] * field.z[cell_z + 1]

    return field_at_point 
end

function ParticleToGrid(part_pos::Float64, system::System, stagger::Int64)

    x_min = system.x_min
    dx = system.dx
    gc = system.gc

    cell_x_r = (part_pos - x_min) / dx
    if stagger == c_stag_centre
        cell_x = floor(Int64, cell_x_r + 0.5)
        cell_frac_x = Float64(cell_x) - cell_x_r
    elseif stagger == c_stag_right
        cell_x = floor(Int64, cell_x_r)
        cell_frac_x = Float64(cell_x) - cell_x_r
    end

    # Cell position in field grid
    cell_x += gc + 1

    # Particle-to- grid (or viceversa) weighting
    cf2 = cell_frac_x^2
    gx = zeros(Float64, 3)
    gx[1] = 0.25 + cf2 + cell_frac_x
    gx[2] = 1.5 - 2.0 * cf2
    gx[3] = 0.25 + cf2 - cell_frac_x
    gx .*= 0.5

    return cell_x, gx
end

function RealocateParticlesToGridList!(species_list::Vector{Species}, system::System)
    x_min = system.x_min
    dx = system.dx
    ncells = system.ncells
    for species in species_list
        if species.is_background_species
            continue
        end

        # Loop until ncells-1 because particle_grid_list is ncells-1 long
        for i in range(1,ncells-1,step=1)

            # Cell boundaries 
            grid_min = (i - 1.0)*dx + x_min
            grid_max = i*dx + x_min 
            indexes = findall(
                x -> (x.pos >= grid_min) & (x.pos <= grid_max),
                species.particle_list)
            species.particle_grid_list[i] = splice!(species.particle_list, indexes)
        end

        # Check that main particle list is empty
        if length(species.particle_list) > 0
            message = "Main particle list of " * species.name * " is not empty"
            PrintWarningMessage(system, message)
        end
    end
end

function RealocateParticlesToMainList!(species_list::Vector{Species})

    for species in species_list
        if !species.is_background_species
            species.particle_list = reduce(vcat, species.particle_grid_list)
            species.particle_count = length(species.particle_list)
        end
    end
end

function GetUnits!(var::Union{String,SubString{String}}; symb='_'::Char)
    # Identify units definition
    units_fact = 1.0
    units_index = findlast(symb, var)
    if !(units_index === nothing)
        units_str = lowercase(var[units_index+1:end])
        match_flag = false
        if (units_str=="tev")
            # this eV refers to TEMPERATURE!
            units_fact = 1.0/K_to_eV
            match_flag = true
        elseif (units_str=="ev")
            # this eV refers to ENERGY 
            units_fact = e
            match_flag = true
        elseif (units_str=="mtorr")
            units_fact = 0.13332237
            match_flag = true
        elseif (units_str=="kw" || units_str=="khz")
            units_fact = 1.e3
            match_flag = true
        elseif (units_str=="mhz" || units_str=="mw")
            units_fact = 1.e6
            match_flag = true
        elseif (units_str=="ghz")
            units_fact = 1.e9
            match_flag = true
        elseif (units_str=="microns")
            units_fact = 1.e-6
            match_flag = true
        elseif (units_str=="nm")
            units_fact = 1.e-9
            match_flag = true
        elseif (units_str=="sccm")
            # The units_fact still needs to be divided by system.V, however,
            # just in case the volume changes, this is done in FunctionTerms.jl
            ns = 2.686780111798444e25 # Standard density at Ps = 101325 Pa and Ts = 273.15 K
            units_fact = 1.e-6 / 60 * ns
            match_flag = true
        end
        if match_flag
            var = strip(var[1:units_index-1])
        end
    end
    return units_fact, var
end

function parse_boundary(var::Union{String,SubString{String}})

    bc_id = c_error
    if var == "x_min" || var == "left"
        bc_id = c_bc_x_min
    elseif var == "x_max" || var == "right"
        bc_id = c_bc_x_max
    end
    return bc_id
end

function linear_interpolation(data::Vector{Float64}, x::Float64)
    x_min = data[1]
    x_max = data[2]
    y_min = data[3]
    y_max = data[4]
    m = (y_max - y_min) / (x_max - x_min)
    n = y_min - x_min * m
    y = x * m + n 
    return y
end

end