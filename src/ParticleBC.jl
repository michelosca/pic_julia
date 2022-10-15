module ParticleBC

using SharedData: Particle, Species, System
using Constants: c_bc_open, c_bc_periodic

function ParticleBoundaryConditions!(i::Int64, species::Species, system::System) 
    x_min = system.x_min
    x_max = system.x_max

    part = species.particle_list[i]
    if part.pos < x_min
        bc_part = system.bc_part_min
        # Particles leaving through the left boundary
        if bc_part == c_bc_open
            deleteat!(species.particle_list,i)
            species.particle_count -= 1
            i -= 1
        elseif bc_part == c_bc_periodic
            part.pos += system.Lx
        end
    elseif part.pos > x_max
        bc_part = system.bc_part_max
        # Particles leaving through the right boundary
        if bc_part == c_bc_open
            deleteat!(species.particle_list,i)
            species.particle_count -= 1
            i -= 1
        elseif bc_part == c_bc_periodic
            part.pos -= system.Lx
        end
    end

end

end