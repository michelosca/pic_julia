module ParticleIntegrator

using SharedData: Species, System, Field
using LinearAlgebra: cross
using Tools: InterpolateFieldToPoint
using ParticleBC: ParticleBoundaryConditions!

function IntegrateParticlesPhaseSpace!(species_list::Vector{Species},
    system::System, electric_field::Field, magnetic_field::Field)

    dt = system.dt

    for species in species_list
        if species.is_background_species
            continue
        end
        #print("Species ", species.name,"\n")

        mass = species.mass
        charge = species.charge
        dt2 = dt * 0.5
        cm = charge/mass
        cm_dt2 = dt2 * cm

        i = 1
        for part in species.particle_list
            pos = part.pos
            vel = part.vel
            #print(" - Part pos ", pos,"\n")
            #print(" - Part vel ", vel,"\n")

            e_field = InterpolateFieldToPoint(electric_field, pos, system)
            b_field = InterpolateFieldToPoint(magnetic_field, pos, system)
            #print(" - E field ", e_field,"\n")
            #print(" - B field ", b_field,"\n")

            # 1st half acceleration
            vel = vel .+ cm_dt2 * e_field
            #print(" - After 1st half-accel. part vel ", vel,"\n")

            # Rotation operators
            t_rot = b_field * cm_dt2
            t_rot2 = sum(t_rot.^2)
            s_rot = t_rot * 2.0/(1.0 + t_rot2 )
            vel_intermediate = vel .+ cross(vel, t_rot)
            vel = vel .+ cross(vel_intermediate, s_rot)
            #print(" - After rotation part vel ", vel,"\n")

            # 1st half acceleration
            vel = vel .+ cm_dt2 * e_field

            # Update particle's position and velocity
            part.pos += vel[1] * dt
            part.vel = vel
            #print(" - New part pos ", pos,"\n")
            #print(" - New part vel ", vel,"\n\n")

            # Check BC and push i-counter
            remove_flag = ParticleBoundaryConditions!(part, system)
            if remove_flag
                delete!(species.particle_list, i)
            else
                i += 1
            end

        end
    end
end

function ParticlesPhasePushBack!(species_list::Vector{Species},
    system::System, electric_field::Field, magnetic_field::Field)

    dt = -system.dt * 0.5

    for species in species_list
        if species.is_background_species
            continue
        end
        #print("Species ", species.name,"\n")

        mass = species.mass
        charge = species.charge
        dt2 = dt * 0.5
        cm = charge/mass
        cm_dt2 = dt2 * cm

        for part in species.particle_list
            pos = part.pos
            vel = part.vel
            #print(" - Part pos ", pos,"\n")
            #print(" - Part vel ", vel,"\n")

            e_field = InterpolateFieldToPoint(electric_field, pos, system)
            b_field = InterpolateFieldToPoint(magnetic_field, pos, system)
            #print(" - E field ", e_field,"\n")
            #print(" - B field ", b_field,"\n")

            # 1st half acceleration
            vel = vel .+ cm_dt2 * e_field
            #print(" - After 1st half-accel. part vel ", vel,"\n")

            # Rotation operators
            t_rot = b_field * cm_dt2
            t_rot2 = sum(t_rot.^2)
            s_rot = t_rot * 2.0/(1.0 + t_rot2 )
            vel_intermediate = vel .+ cross(vel, t_rot)
            vel = vel .+ cross(vel_intermediate, s_rot)
            #print(" - After rotation part vel ", vel,"\n")

            # 1st half acceleration
            vel = vel .+ cm_dt2 * e_field

            # Update particle's velocity
            part.vel = vel
            #print(" - New part vel ", vel,"\n\n")
        end
    end
end

end