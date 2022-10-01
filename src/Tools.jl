module Tools

using SharedData: Species, System, Field
using Constants: c_bc_open, c_bc_periodic
using Constants: gc

function InterpolateParticleToGrid!(field::Vector{Float64}, part_pos::Float64, system::System)

    cell_x, gx = ParticleToGrid(part_pos, system)

    # Add particle values to field grid
    field[cell_x - 1] += gx[1]
    field[cell_x]     += gx[2]
    field[cell_x + 1] += gx[3]
end

function InterpolateGridToPoint(field::Vector{Float64}, part_pos::Float64, system::System)


    cell_x, gx = ParticleToGrid(part_pos, system)

    value_at_point = 0.0 
    value_at_point += gx[1] * field[cell_x - 1]
    value_at_point += gx[2] * field[cell_x]
    value_at_point += gx[3] * field[cell_x + 1]

    return value_at_point
end

function ParticleToGrid(part_pos::Float64, system::System)

    x_min = system.x_min
    dx = system.dx

    cell_x_r = (part_pos - x_min) / dx
    cell_x = floor(Int64, cell_x_r + 0.5)
    cell_frac_x = Float64(cell_x) - cell_x_r

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