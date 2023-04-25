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

module InputBlock_Fields

using Constants: c_error
using Constants: c_bc_open, c_bc_periodic
using Tools: GetUnits!
using PrintModule: PrintErrorMessage


function StartFile_Fields(read_step::Int64) 

    errcode = c_error
    
    if read_step == 1
        errcode = 0
    elseif read_step == 2
        errcode = 0
    end

    return errcode
end


function StartFieldsBlock!(read_step::Int64, system::System)

    errcode = c_error

    if (read_step == 1)
        errcode = 0
    elseif read_step == 2
        errcode = 0
    end
    return errcode
end


function ReadFieldsEntry!(name::SubString{String}, var::SubString{String},
    read_step::Int64)

    errcode = c_error 

    if (read_step == 1)
        units_fact, name = GetUnits!(name)
        errcode = 0
    elseif read_step == 2
        errcode = 0
    end
    return errcode 
end


function EndSystemBlock!(read_step::Int64)

    errcode = 0

    return errcode
end
    

function EndFile_System!(read_step::Int64)

    errcode = 0 

    return errcode
end

end