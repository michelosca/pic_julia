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

module PrintModule

using SharedData: System, Species, Collision, CollisionGroup
using SharedData: Waveform
using Constants: c_bc_x_min, c_bc_x_max
using Constants: e
using Printf

global tab_str1 = "  "
global tab_str2 = tab_str1 * tab_str1
global tab_str3 = tab_str1 * tab_str1 * tab_str1

function PrintErrorMessage(system::System, message::String)

    open(system.log_file,"a") do file
        @printf(file, "***ERROR*** %s\n", message)
    end
    @printf("***ERROR*** %s\n", message)
end

function PrintWarningMessage(system::System, message::String)

    open(system.log_file,"a") do file
        @printf(file, "***WARNING*** %s\n", message)
    end
    @printf("***WARNING*** %s\n", message)
end

function PrintMessage(system::System, message::String)

    open(system.log_file,"a") do file
        @printf(file, "%s\n", message)
    end
    @printf("%s\n", message)
end

function PrintSpecies(system::System, species::Species)

    open(system.log_file, "a") do file
        @printf(file, "Species name: %s\n", species.name)
        @printf(file, "%s- id: %i\n", tab_str1, species.id)
        @printf(file, "%s- mass: %g kg\n", tab_str1, species.mass)
        @printf(file, "%s- charge: %g C\n", tab_str1, species.charge)
        @printf(file, "%s- initial density: %g m^-3\n", tab_str1, species.init_dens)
        @printf(file, "%s- initial temp: %g K\n", tab_str1, species.init_temp)
        @printf(file, "%s- part. weight: %g\n", tab_str1, species.weight)
        @printf(file, "%s- particle count: %i\n", tab_str1, species.particle_count)
        @printf(file, "%s- is background: %s\n", tab_str1, species.is_background_species)
    end
end

function PrintWaveform(system::System, w::Waveform)

    open(system.log_file,"a") do file
        @printf(file,"Waveform\n")
        @printf(file, "%s- Amplitude: %g V\n", tab_str1, w.amp)
        @printf(file, "%s- Frequency: %g MHz\n", tab_str1, w.freq/1.e6)
        if w.boundary == c_bc_x_max
            b_str = "x-max"
        elseif w.boundary == c_bc_x_min
            b_str = "x-min"
        end
        @printf(file, "%s- Boundary: %s\n", tab_str1, b_str)
        @printf(file, "%s- Start time: %g s\n", tab_str1, w.t_start)
        @printf(file, "%s- End time: %g s\n", tab_str1, w.t_end)
    end

end


function PrintCollision(system::System, c::Collision)
    open(system.log_file, "a") do file
        @printf(file, "%sCollision name: %s\n",tab_str1, c.name)
        coll_str = ""
        r_len = length(c.reactants)
        for (i,s) in enumerate(c.reactants)
            coll_str *= s.name
            if i < r_len
                coll_str *= " + "
            end
        end
        coll_str *= " -> "
        p_len = length(c.products)
        for (i,s) in enumerate(c.products)
            coll_str *= s.name
            if i < p_len
                coll_str *= " + "
            end
        end
        @printf(file, "%s- Reaction equation: %s\n",tab_str2, coll_str)
        species_str = ""
        s_len = length(c.species)
        for (i,s) in enumerate(c.species)
            species_str *= tab_str3 *"-> " * s.name
            species_str *= @sprintf(" - balance: %i", c.species_balance[i])
            if i < s_len
                species_str *= "\n"
            end
        end
        @printf(file, "%s- Involved species:\n%s\n",tab_str2, species_str)
        
        @printf(file, "%s- Energy threshold %g eV\n",tab_str2, c.energy_threshold / e)
        @printf(file, "%s- Energy data length %i\n",tab_str2, length(c.energy_data))
        @printf(file, "%s- Energy units %g\n",tab_str2, c.energy_units)
        @printf(file, "%s- Cross-section data length %i\n",tab_str2, length(c.cross_section_data))
        @printf(file, "%s- Cross section units %g\n",tab_str2, c.cross_section_units)
        #for (energy, sigma) in zip(c.energy_data, c.cross_section_data)
        #    @printf(file, "%s%10.4e %10.4e\n",tab_str3, tab_str,energy/e,sigma)
        #end
    end
end

function PrintCollisionGroup(system::System, cgroup::CollisionGroup)

    open(system.log_file, "a") do file
        @printf(file, "Collision group species: ")
        for s in cgroup.colliding_species
            @printf(file, "%s ", s.name)
        end
        @printf(file, "\n")

        @printf(file, "%s- Max. super-particle weight: %10.4e\n", tab_str1, cgroup.part_weight_max)
        @printf(file, "%s- Max. g*sigma:               %10.4e\n", tab_str1, cgroup.gsigma_max)
        @printf(file, "%s- Reduces mass:               %10.4e\n", tab_str1, cgroup.reduced_mass)
    end

    for c in cgroup.collision_list
        PrintCollision(system,c)
    end

    open(system.log_file, "a") do file
        @printf(file, "\n")
    end

end

function PrintPICJlabel(system::System, print_to_log_file::Bool)

    if print_to_log_file
        open(system.log_file,"a") do file
            @printf(file, "|        _           _       _ _(_)_     |  Particle-In-Cell code for Low Temperature Plasmas\n")
            @printf(file, "|       (_)         (_)     | (_) (_)    |\n")
            @printf(file, "|  ___   _  ___      _ _   _| |_  __ _   |  Documentation: https://github.com/michelosca/pic_julia.git\n") 
            @printf(file, "| |  _ \\| |/ __|    | | | | | | |/ _` |  |\n")
            @printf(file, "| | |_) | | (__     | | |_| | | | (_| |  |  Developer: M. Osca Engelbrecht, michel.osca@protonmail.com\n") 
            @printf(file, "| |  __/|_|\\___|   _/ |\\__'_|_|_|\\__'_|  |\n")
            @printf(file, "| |_|             |__/                   |\n\n")
        end
    else
        @printf("|        _           _       _ _(_)_     |  Particle-In-Cell code for Low Temperature Plasmas\n")
        @printf("|       (_)         (_)     | (_) (_)    |\n")
        @printf("|  ___   _  ___      _ _   _| |_  __ _   |  Documentation: https://github.com/michelosca/pic_julia.git\n") 
        @printf("| |  _ \\| |/ __|    | | | | | | |/ _` |  |\n")
        @printf("| | |_) | | (__     | | |_| | | | (_| |  |  Developer: M. Osca Engelbrecht, michel.osca@protonmail.com\n") 
        @printf("| |  __/|_|\\___|   _/ |\\__'_|_|_|\\__'_|  |\n")
        @printf("| |_|             |__/                   |\n\n")
    end
end

end