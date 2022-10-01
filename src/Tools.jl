module Tools

using SharedData: Species, System, Field
using Constants: c_bc_open, c_bc_periodic

function GetNumberDensity(species::Species, system::System)
    idx = 1.0 / system.dx
    ncells = system.ncells

    dens = zeros(Float64, ncells)
    for part in species.particle_list
        IntegrateParticleToField!(dens, part, system) 
    end
    dens .*= idx * species.weight
    dens[1] *= 2.0
    dens[ncells] *= 2.0

    return dens
end

function IntegrateParticleToField!(field::Vector{Float64}, part::Particle, system::System)

    x_min = system.x_min
    ncells = system.ncells
    bc_part = system.bc_part

    cell_x_r = (part.pos - x_min) * idx
    cell_x = floor(Int64, cell_x_r + 0.5)
    cell_frac_x = Float64(cell_x) - cell_x_r
    gf = ParticleToGrid(cell_frac_x)
    cell_x += 1

    if cell_x == 1
        if bc_part == c_bc_open
            field[cell_x + 1] += gf[1]
        elseif bc_part == c_bc_periodic
            field[ncells] += gf[1]
        end
        field[cell_x]     += gf[2]
        field[cell_x + 1] += gf[3]
    elseif cell_x == ncells
        field[cell_x - 1] += gf[1]
        field[cell_x]     += gf[2]
        if bc_part == c_bc_open
            field[cell_x - 1] += gf[3]
        elseif bc_part == c_bc_periodic
            field[1] += gf[3]
        end
    else
        field[cell_x - 1] += gf[1]
        field[cell_x]     += gf[2]
        field[cell_x + 1] += gf[3]
    end
end

function ParticleToGrid(cell_frac_x::Float64)
    cf2 = cell_frac_x^2
    gx = zeros(Float64, 3)
    gx[1] = 0.25 + cf2 + cell_frac_x
    gx[2] = 1.5 - 2.0 * cf2
    gx[3] = 0.25 + cf2 - cell_frac_x
    return gx
end

function GetFieldAtPoint(field::Vector{Float64}, gf::Vector{Float64},
    cell_x::Int64, ncells::Int64, bc_field::Int64)

    field_at_pos = 0.0 
    if cell_x == 1
        if bc_field == c_bc_open
            field_at_pos += gf[1] * field[cell_x + 1]
        elseif bc_field == c_bc_periodic
            field_at_pos += gf[1] * field[ncells]
        end
        field_at_pos += gf[2] * field[cell_x]
        field_at_pos += gf[3] * field[cell_x + 1]
    elseif cell_x == ncells
        field_at_pos += gf[1] * field[cell_x - 1]
        field_at_pos += gf[2] * field[cell_x]
        if bc_field == c_bc_open
            field_at_pos += gf[3] * field[cell_x - 1]
        elseif bc_field == c_bc_periodic
            field_at_pos += gf[3] * field[1]
        end
    else
        field_at_pos += gf[1] * field[cell_x - 1]
        field_at_pos += gf[2] * field[cell_x]
        field_at_pos += gf[3] * field[cell_x + 1]
    end
    return field_at_pos
end

function RealocateParticlesToGridList(species_list::Vector{Species})
    x_min = system.x_min
    dx = system.dx
    ncells = system.ncells
    for species in species_list
        for i in range(1,ncells,step=1)
            grid_min = (i - 0.5)*dx + x_min
            grid_max = (i + 0.5)*dx + x_min 
            indexes = findall(
                x -> (x.pos >= grid_min) & (x.pos <= grid_max),
                species.particle_list)
            species.particle_grid_list[i] = splice!(species.particle_list, indexes)
        end
    end
    # Check that main particle list is empty
    if length(species.particle_list) > 0
        print("***WARNING*** Main particle list of ",species.name," is not empty\n")
    end
end

function RealocateParticlesToMainList(species_list::Vector{Species})

    for species in species_list
        # Check that main particle list is empty
        if length(species.particle_list) > 0
            print("***WARNING*** Main particle list of ",species.name," is not empty\n")
        end

        for part_list in species.particle_grid_list
            push!(species.particle_list, part_list)
        end
    end
end

end