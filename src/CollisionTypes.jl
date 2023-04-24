# Copyright (C) 2023 Michel Osca Engelbrecht
#
# This file is part of PIC Julia.
#
# PIC Julia is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PIC Julia is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GM Julia. If not, see <https://www.gnu.org/licenses/>.

module CollisionTypes

using SharedData: CollisionGroup, Collision, Particle
using Random: rand
using LinearAlgebra: norm

function ChargeExchange!(collgroup::CollisionGroup, coll::Collision,
    part1::Particle, part2::Particle, g::Float64)

    vel1 = part1.vel
    part1.vel = part2.vel
    part2.vel = vel1

#    print("Charge exchange\n")
end


function ElasticScattering!(collgroup::CollisionGroup, coll::Collision,
    part1::Particle, part2::Particle, g::Float64)
    
    # Species 1
    s1 = collgroup.colliding_species[1]
    m1 = s1.mass

    # Species 2
    s2 = collgroup.colliding_species[2]
    m2 = s2.mass

    # Total mass
    m_total = m1 + m2
    # Reduced mass
    mu = collgroup.reduced_mass

    # Centre of mass (CM) velocity
    vel_cm = (part1.vel * m1 + part2.vel * m2) / m_total

    # Isotropic scattering in CM frame of reference
    #  - Unit random vector
    r3 = rand(3)
    r = r3 / norm(r3)
    #  - Post collision g
    g_postcoll = g * r

    # Particle scattering
    part1.vel = vel_cm .+ mu * g_postcoll
    part2.vel = vel_cm .- mu * g_postcoll
#    print("Elastic\n")
end


function InelasticScattering!(collgroup::CollisionGroup, coll::Collision,
    part1::Particle, part2::Particle, g::Float64)
    
    # Species 1
    s1 = collgroup.colliding_species[1]
    m1 = s1.mass

    # Species 2
    s2 = collgroup.colliding_species[2]
    m2 = s2.mass

    # Total mass
    m_total = m1 + m2
    # Reduced mass
    mu = collgroup.reduced_mass

    # Centre of mass (CM) velocity
    vel_cm = (part1.vel * m1 + part2.vel * m2) / m_total

    # Isotropic scattering in CM frame of reference
    #  - Unit random vector
    r3 = rand(3)
    r = r3 / norm(r3)
    #  - Post collision g
    E = coll.energy_threshold
    g_postcoll = sqrt(g*g - 2.0 * E / mu) * r

    # Particle scattering
    part1.vel = vel_cm .+ mu * g_postcoll
    part2.vel = vel_cm .- mu * g_postcoll

#    print("Excitation\n")
end


function Ionization!(collgroup::CollisionGroup, coll::Collision,
    part1::Particle, part2::Particle, g::Float64)

    # Species 1
    s1 = collgroup.colliding_species[1]
    m1 = s1.mass

    # Species 2
    s2 = collgroup.colliding_species[2]
    m2 = s2.mass

    # Total mass
    m_total = m1 + m2
    # Reduced mass
    mu = collgroup.reduced_mass

    # Post-collision energy excess
    E_threshold = coll.energy_threshold
    E_excess = 0.5 * mu * g * g - E_threshold
    u_e_excess = sqrt(2.0 * E_excess / m1)

    # Electron post-coll speed (CM frame of reference)
    R1 = rand()
    u_e1 = sqrt(R1) * u_e_excess
    u_e2 = sqrt(1-R1) * u_e_excess

    # Centre of mass (CM) velocity
    vel_cm = (part1.vel * m1 + part2.vel * m2) / m_total

    # Existing electrons
    r3 = rand(3)
    r = r3 / norm(r3)
    part1.vel = vel_cm .+ u_e1 * r

    # New species 
    #  - Remove neutral
    #  - New electron
    r3 = rand(3)
    r = r3 / norm(r3)
    new_electron = Particle()
    new_electron.pos = neutral.pos
    new_electron.vel = vel_cm .+ u_e2 * r
    #  - New ion
    mass_ratio = m_electron / m_neutral
    new_ion = Particle()
    new_ion.pos = neutral.pos
    new_ion.vel = vel_cm * (1 - mass_ratio) .- m_ratio * (u_e1 .+ u_e2)

end

end