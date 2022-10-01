module GridData

using SharedData: System, Species
using Constants: gc, c_bc_periodic, c_bc_open
using Tools: InterpolateParticleToGrid! 

function GetNumberDensity(species::Species, system::System)

    idx = 1.0 / system.dx
    ncells = system.ncells

    dens = zeros(Float64, ncells + gc*2)
    for part in species.particle_list
        InterpolateParticleToGrid!(dens, part.pos, system) 
    end
    ApplyBCtoGridData!(dens, system)
    dens .*= idx * species.weight

    return dens
end

function GetTotalChargeDensity(species_list::Vector{Species}, system::System)

    charge_density = zeros(Float64, system.ncells + gc*2)
    for species in species_list
        if species.is_background_gas
            continue
        end
        dens = GetNumberDensity(species, system)
        charge_density .+= dens .* species.charge
    end

    return charge_density
end

function ApplyBCtoGridData!(grid::Vector{Float64}, system::System)
    bc = system.bc_field
    cell_min = gc+1
    cell_max = system.ncells + gc
    if bc == c_bc_periodic 
        for i in range(1,gc,step=1)
            grid[gc+i] += grid[ncells+i]
            grid[ncells+1-i] = grid[i]
        end
    elseif bc == c_bc_open
        for i in range(1,gc,step=1)
            grid[cell_min+i] += grid[cell_min-i]
            grid[cell_max-i] += grid[cell_max+i]
        end

        # Last and first grid have a dx/2 width
        grid[cell_min] *= 2.0
        grid[cell_max] *= 2.0
    end
end


end