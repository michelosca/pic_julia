module GridData

using SharedData: System, Species
using Constants: c_bc_periodic, c_bc_open
using Tools: InterpolateParticleToGrid! 

function GetNumberDensity!(species::Species, system::System)

    idx = 1.0 / system.dx

    dens = zeros(Float64, system.ncells_total)
    for part in species.particle_list
        InterpolateParticleToGrid!(dens, part.pos, system) 
    end
    ApplyBCtoGridData!(dens, system)
    species.dens = dens * idx * species.weight
end

function GetTotalChargeDensity(species_list::Vector{Species}, system::System)

    charge_density = zeros(Float64, system.ncells_total)
    for species in species_list
        if species.is_background_gas
            continue
        end
        GetNumberDensity!(species, system)
        charge_density .+= species.dens * species.charge
    end

    return charge_density
end

function ApplyBCtoGridData!(grid::Vector{Float64}, system::System)

    cell_min = system.cell_min 
    cell_max = system.cell_max
    gc = system.gc

    bc = system.bc_field
    if bc == c_bc_periodic 
        # Past data at ghost cells to the opposite boundary
        for i in range(1,gc,step=1)
            grid[cell_min+i] += grid[cell_max+i]
            grid[cell_max-i] += grid[cell_min-i]
        end
        grid[cell_min] += grid[cell_max]
        grid[cell_max] = grid[cell_min]
    elseif bc == c_bc_open
        # Fold ghost cells data back into the adjacent boundary 
        for i in range(1,gc-1,step=1)
            grid[cell_min+i] += grid[cell_min-i]
            grid[cell_min-i] = 0.0
            grid[cell_max-i] += grid[cell_max+i]
            grid[cell_max+i] = 0.0
        end

        # Last and first grid have a dx/2 width
        grid[cell_min] *= 2.0
        grid[cell_max] *= 2.0
    end
end


end