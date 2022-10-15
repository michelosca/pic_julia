# Copyright (C) 2021 Michel Osca Engelbrecht
#
# This file is part of GM Julia.
#
# GM Julia is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GM Julia is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GM Julia. If not, see <https://www.gnu.org/licenses/>.

module PrintModule

using SharedData: System, Species
using Printf


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
        @printf(file, "%s", message)
    end
    @printf("%s", message)
end

function PrintSpecies(species::Species)
    @printf("Species name: %s\n", species.name)
    @printf("  - id: %i\n", species.id)
    @printf("  - mass: %g kg\n", species.mass)
    @printf("  - charge: %g C\n", species.charge)
    @printf("  - part. weight: %g\n", species.weight)
    @printf("  - particle count: %i\n", species.particle_count)
    @printf("  - is background: %s\n", species.is_background_species)
end

end