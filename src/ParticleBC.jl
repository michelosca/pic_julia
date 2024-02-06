module ParticleBC

using SharedData: Particle, Species, System
using Constants: c_bc_open, c_bc_periodic

function ParticleBoundaryConditions!(part::Particle, system::System) 
    x_min = system.x_min
    x_max = system.x_max

    remove_flag = false
    if part.pos < x_min
        bc_part = system.bc_part_min
        # Particles leaving through the left boundary
        if bc_part == c_bc_open
            remove_flag = true
        elseif bc_part == c_bc_periodic
            part.pos += system.Lx
        end
    elseif part.pos > x_max
        bc_part = system.bc_part_max
        # Particles leaving through the right boundary
        if bc_part == c_bc_open
            remove_flag = true
        elseif bc_part == c_bc_periodic
            part.pos -= system.Lx
        end
    end

    return remove_flag
end

end