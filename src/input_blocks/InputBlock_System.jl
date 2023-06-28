# Copyright (C) 2021 Michel Osca Engelbrecht
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
# along with PIC Julia. If not, see <https://www.gnu.org/licenses/>.

module InputBlock_System

using Constants: c_error
using Constants: gc_triangle
using Constants: c_bc_open, c_bc_periodic
using SharedData: System, Species
using Tools: GetUnits!
using Dates
using PrintModule: PrintErrorMessage


function StartFile_System!(read_step::Int64, species_list::Vector{Species}, system::System) 

    errcode = c_error
    
    if read_step == 1
        errcode = 0
    elseif read_step == 2

        # At this point species are already defined and time-step conditions can be set
        errcode = SetTimeStep!(system, species_list)
        if errcode == c_error
            PrintErrorMessage(system, "Setting up time step")
            errcode = c_error
            return errcode 
        end

        system.step_start = round(Int64, (system.t_start)/system.dt)
        system.step_end = round(Int64, (system.t_end)/system.dt) - system.step_start
        system.step = system.step_start

        errcode = 0
    end

    return errcode
end


function StartSystemBlock!(read_step::Int64, system::System)

    errcode = c_error

    if (read_step == 1)
        system.ncells = -1 #101
        system.gc = gc_triangle
        system.cell_min = -1 #system.gc + 1
        system.cell_max = -1 #system.gc + system.ncells
        system.ncells_total = -1 #system.ncells + 2*system.gc
        system.x_min = 0.0
        system.x_max = 0.0 #0.1
        system.Lx = 0.0 #system.x_max - system.x_min
        system.dx = 0.0 #system.Lx / (system.ncells-1)
        
        system.t_start = 0.0
        system.t_end = 0.0 #1.12e-8
        system.time = 0.0 #system.t_start
        system.dt = 0.0 #1.e-10
        system.step = 0
        system.step_start = 0
        system.step_end = 0
    
        system.bc_field_min = -1 #c_bc_periodic #open #
        system.bc_field_max = -1 #c_bc_periodic #open #
        system.bc_part_min = -1 #c_bc_periodic #open #
        system.bc_part_max = -1 #c_bc_periodic #open #

        system.mcc = false

        # Generate a log file name
        now_stamp = Dates.now()
        log_num = Dates.format(now_stamp, "yyyymmddHHMMSS")
        system.log_file = "PICJ_run_" * log_num * ".log"

        errcode = 0
    elseif read_step == 2
        errcode = 0
    end
    return errcode
end


function ReadSystemEntry!(name::SubString{String}, var::SubString{String},
    read_step::Int64, system::System)

    errcode = c_error 

    if (read_step == 1)
        
        units_fact, name = GetUnits!(name)

        # Identify variable
        lname = lowercase(name)
        if (name=="dx" || lname=="cell_width")
            system.dx = parse(Float64,var) * units_fact
            errcode = 0
        elseif lname=="x_min"
            system.x_min = parse(Float64, var) * units_fact
            errcode = 0
        elseif lname=="x_max"
            system.x_max = parse(Float64, var) * units_fact
            errcode = 0
        elseif lname=="lx"
            system.Lx = parse(Float64, var) * units_fact
            errcode = 0
        elseif lname=="ncells" || lname=="cells"
            system.ncells = parse(Int64, var) 
            errcode = 0
        elseif (name=="integration_shape" || lname=="interpolation_shape")
            lvar = lowercase(var)
            if (lvar == "triangle")
                system.gc = gc_triangle 
                errcode = 0
            end
        elseif (lname=="particle_bc_min" || lname=="part_bc_min")
            lvar = lowercase(var)
            if (lvar == "open")
                system.bc_part_min = c_bc_open
            elseif (lvar == "periodic")
                system.bc_part_min = c_bc_periodic
            end
            errcode = 0
        elseif (lname=="particle_bc_max" || lname=="part_bc_max")
            lvar = lowercase(var)
            if (lvar == "open")
                system.bc_part_max = c_bc_open
            elseif (lvar == "periodic")
                system.bc_part_max = c_bc_periodic
            end
            errcode = 0
        elseif (lname=="particle_bc" || lname=="part_bc")
            lvar = lowercase(var)
            if (lvar == "open")
                system.bc_part_min = c_bc_open
                system.bc_part_max = c_bc_open
            elseif (lvar == "periodic")
                system.bc_part_min = c_bc_periodic
                system.bc_part_max = c_bc_periodic
            end
            errcode = 0
        elseif lname=="field_bc_min"
            lvar = lowercase(var)
            if (lvar == "open")
                system.bc_field_min = c_bc_open
            elseif (lvar == "periodic")
                system.bc_field_min = c_bc_periodic
            end
            errcode = 0
        elseif lname=="field_bc_max"
            lvar = lowercase(var)
            if (lvar == "open")
                system.bc_field_max = c_bc_open
            elseif (lvar == "periodic")
                system.bc_field_max = c_bc_periodic
            end
            errcode = 0
        elseif lname=="field_bc"
            lvar = lowercase(var)
            if (lvar == "open")
                system.bc_field_min = c_bc_open
                system.bc_field_max = c_bc_open
            elseif (lvar == "periodic")
                system.bc_field_min = c_bc_periodic
                system.bc_field_max = c_bc_periodic
            end
            errcode = 0
        elseif lname=="dt"
            system.dt = parse(Float64, var) * units_fact
            errcode = 0
        elseif lname=="t_end"
            system.t_end = parse(Float64, var) * units_fact
            errcode = 0
        elseif lname=="step_end"
            system.step_end = parse(Int64, var)
            errcode = 0
        end
    elseif read_step == 2
        errcode = 0
    end
    return errcode 
end


function EndSystemBlock!(read_step::Int64, system::System)

    errcode = 0

    if read_step == 1

        ### Spatial check
        if system.Lx == 0.0
            system.Lx = system.x_max - system.x_min 
            if system.Lx <= 0.0
                PrintErrorMessage(system, "Simulation domain length must be > 0")
                errcode = c_error
                return errcode 
            end
        else
            system.x_min = 0.0
            system.x_max = system.Lx
        end

        if system.ncells == -1
            if system.dx > 0.0
                system.ncells = round(Int64, system.Lx/system.dx) + 1
                system.x_max = (system.ncells-1) * system.dx + system.x_min
                system.Lx = system.x_max - system.x_min 
            else
                PrintErrorMessage(system, "Number of cells must be defined in advance")
                errcode = c_error
                return errcode 
            end
        else
            system.dx = system.Lx / system.ncells 
            system.ncells += 1

            system.cell_min = system.gc + 1
            system.cell_max = system.gc + system.ncells
            system.ncells_total = system.ncells + 2*system.gc
        end

        system.ncells_total = system.ncells + system.gc * 2
        system.cell_min = system.gc + 1
        system.cell_max = system.ncells + system.gc

        ### Temporal check
        if system.t_end == 0
            system.t_end = system.step_end * system.dt
        end

        if system.step_end == 0
            system.step_end = round(Int64, system.t_end / system.dt)
        end

        if system.t_end <= system.t_start
            PrintErrorMessage(system, "Simulation times must be t_end > t_start")
            errcode = c_error
            return errcode
        end

        system.time = system.t_start

    end

    return errcode
end
    

function EndFile_System!(read_step::Int64, system::System)

    errcode = 0 

    return errcode
end

function SetTimeStep!(system::System, species_list::Vector{Species})
    return 0
end

end