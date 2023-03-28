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

module InputBlock_

using Constants: c_error
using Constants: c_bc_open
using PrintModule: PrintErrorMessage


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
        errcode = 0
    end
    return errcode
end


function ReadWaveformEntry!(name::SubString{String}, var::SubString{String},
    read_step::Int64, waveform_list::Vector{Waveform})

    errcode = c_error 

    if (read_step == 1)
        units_fact, name = GetUnits!(name)
        errcode = 0
    elseif read_step == 2
        errcode = 0
    end
    return errcode 
end


function EndWaveformBlock!(read_step::Int64, waveform_list::Vector{Waveform},
    system::System)

    errcode = 0

    return errcode
end
    

function EndFile_Waveform!(read_step::Int64, waveform_list::Vector{Waveform},
    system::System) 

    errcode = 0 

    return errcode
end

end