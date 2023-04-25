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

module InputBlock_Waveform

using Constants: c_error
using Constants: c_bc_open
using SharedData: System, Waveform
using PrintModule: PrintErrorMessage
using Tools: parse_boundary, GetUnits!
using Printf


function StartFile_Waveform!(read_step::Int64, waveform_list::Vector{Waveform},
    system::System) 

    errcode = c_error
    
    if read_step == 1
        errcode = 0
    elseif read_step == 2
        errcode = 0
    end

    return errcode
end


function StartWaveformBlock!(read_step::Int64, waveform_list::Vector{Waveform},
    system::System)

    errcode = c_error

    if (read_step == 1)
        errcode = 0
    elseif read_step == 2
        if system.bc_field_min == c_bc_open &&
            system.bc_field_max == c_bc_open &&
            system.bc_part_min == c_bc_open &&
            system.bc_part_max == c_bc_open
            waveform = InitWaveform(system)
            push!(waveform_list, waveform)
            errcode = 0
        else
            message = "Waveform requires open field and particle boundary conditions"
            PrintErrorMessage(system, message)
            return errcode
        end
    end
    return errcode
end


function ReadWaveformEntry!(name::SubString{String}, var::SubString{String},
    read_step::Int64, waveform_list::Vector{Waveform})

    errcode = c_error 

    if (read_step == 1)
        errcode = 0
    elseif read_step == 2
        waveform = waveform_list[end]

        units_fact, name = GetUnits!(name)
        if name == "amplitude" || name == "amp"
            waveform.amp = parse(Float64, var) * units_fact
            errcode = 0
        elseif name == "f" || name == "frequency" || name == "freq"
            waveform.freq = parse(Float64, var) * units_fact
            errcode = 0
        elseif name == "waveform" || name == "function" || name == "shape"
            waveform.wavefunction = Meta.parse(var)
            errcode = 0
        elseif name == "boundary"
            waveform.boundary = parse_boundary(var)
            errcode = 0
        elseif name == "t_start"
            waveform.t_start = parse(Float64, var)
            errcode = 0
        elseif name == "t_end"
            waveform.t_end = parse(Float64, var)
            errcode = 0
        end
    end
    return errcode 
end


function EndWaveformBlock!(read_step::Int64, waveform_list::Vector{Waveform},
    system::System)

    errcode = 0

    if read_step == 1
        return errcode
    end

    waveform = waveform_list[end]

    if waveform.boundary == c_error 
        err_message = @sprintf("Waveform boundary location has not been defined")
        PrintErrorMessage(system, err_message)
        return c_error
    end

    if waveform.t_end < waveform.t_start
        err_message = @sprintf("Waveform must have t_start < t_end")
        PrintErrorMessage(system, err_message)
        return c_error
    end

    if waveform.boundary == c_error
        err_message = @sprintf("Waveform boundary not recognized")
        PrintErrorMessage(system, err_message)
        return c_error
    end

    return errcode
end
    

function EndFile_Waveform!(read_step::Int64, waveform_list::Vector{Waveform},
    system::System) 

    errcode = 0 

    return errcode
end


function InitWaveform(system::System)

    waveform = Waveform()
    waveform.amp = 0.0
    waveform.freq = 0.0
    waveform.wavefunction = 1.0
    waveform.boundary = c_error
    waveform.t_start = system.t_start
    waveform.t_end = system.t_end

    return waveform
end

end