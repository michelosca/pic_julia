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

module InputBlock_Output

using Constants: c_error, c_o_all_species, c_o_none_species
using Constants: c_o_density, c_o_potential, c_o_electric_field, c_o_probe
using Constants: c_o_phase_space, c_o_neutral_collisions
using Constants: c_dir_x, c_dir_y, c_dir_z, c_dir_none
using SharedData: Species, System, OutputBlock, OutputDataStruct
using Tools: GetUnits!
using PrintModule: PrintErrorMessage


function StartFile_Output!(read_step::Int64, output_list::Vector{OutputBlock},
    system::System) 

    errcode = 0

    return errcode
end

function StartOutputBlock!(read_step::Int64, output_list::Vector{OutputBlock},
    system::System)

    errcode = c_error

    if read_step == 1
        errcode = 0
    elseif (read_step == 2)
        o_block = OutputBlock()

        o_block.name = "stdout"
        o_block.t_start = system.t_start
        o_block.t_end = system.t_end
        o_block.dt = -1.0

        o_block.step_start = 0 
        o_block.step_end = system.step_end
        o_block.step_jump = -1

        o_block.averaged = false
        o_block.dt_av = -1.0
        o_block.step_av = -1

        o_block.param_list = OutputDataStruct[]

        o_block.file_id = 0
        o_block.zero_pad = 5
        o_block.time_dump= -1.0
        o_block.step_dump = -1
        o_block.step_av_start = -1
        o_block.step_av_end = -1
        o_block.time_av_start = -1.0
        o_block.time_av_end = -1.0

        push!(output_list, o_block)

        errcode = 0
    end
    return errcode
end

function ReadOutputEntry!(name::SubString{String}, var::SubString{String},
    read_step::Int64, output_list::Vector{OutputBlock},
    species_list::Vector{Species}, system::System)

    errcode = c_error 

    if read_step == 1
        errcode = 0

    elseif read_step == 2
        o_block = output_list[end]
        
        units_fact, name = GetUnits!(name)

        if name == "name"
            o_block.name = var 
            return 0
        end

        if name == "zero_pad"
            o_block.zero_pad = parse(Int64,var) 
            return 0
        end

        if name == "t_start"
            o_block.t_start = parse(Float64, var) * units_fact
            return 0
        end
        
        if name == "t_end"
            o_block.t_end = parse(Float64, var) * units_fact
            return 0
        end
        
        if name == "dt"
            o_block.dt = parse(Float64, var) * units_fact
            return 0
        end
        
        if name == "step_jump"
            o_block.step_jump = parse(Int64, var)
            return 0
        end
        
        if name == "step_start"
            o_block.step_start = parse(Int64, var)
            return 0
        end
        
        if name == "step_end"
            o_block.step_end = parse(Int64, var)
            return 0
        end
        
        if name == "averaged"
            o_block.averaged = parse(Bool, var)
            return 0
        end
        
        if (name == "dt_av") || (name == "dt_average")
            o_block.dt_av = parse(Float64, var) * units_fact
            return 0
        end
        
        if name == "step_av"
            o_block.step_av = parse(Int64, var)
            return 0
        end

        if (name == "density") || (name == "number_density")
            p_id = c_o_density
            var_list = SplitVariableString(var)

            data_struct = OutputDataStruct() 
            Initialize_OutputDataStruct!(data_struct, name, p_id, var_list,
                o_block, species_list, system)
            push!(o_block.param_list, data_struct)
            return 0
        end

        if (name == "potential") || (name == "electric_potential")
            p_id = c_o_potential
            var_list = SplitVariableString(var)

            data_struct = OutputDataStruct() 
            Initialize_OutputDataStruct!(data_struct, name, p_id, var_list,
                o_block, species_list, system)
            push!(o_block.param_list, data_struct)
            return 0
        end

        if name == "electric_field"
            p_id = c_o_electric_field
            var_list = SplitVariableString(var)

            data_struct = OutputDataStruct() 
            Initialize_OutputDataStruct!(data_struct, name, p_id, var_list,
                o_block, species_list, system)
            push!(o_block.param_list, data_struct)
            return 0
        end

        if name == "phase_space"
            p_id = c_o_phase_space
            var_list = SplitVariableString(var)

            data_struct = OutputDataStruct() 
            Initialize_OutputDataStruct!(data_struct, name, p_id, var_list,
                o_block, species_list, system)
            push!(o_block.param_list, data_struct)
            return 0
        end
    end
    return errcode 
end

function EndOutputBlock!(read_step::Int64, output_list::Vector{OutputBlock},
    system::System)

    errcode = c_error

    if read_step == 1
        errcode = 0
    elseif read_step == 2

        # Check output parameters
        # - The step parameters will be used during main simulation and therefore
        # they are used as reference

        # Check output dump period
        o_block = output_list[end]
        if o_block.step_jump <= 0
            if o_block.dt > 0.0
                o_block.step_jump = round(Int64, o_block.dt / system.dt)
            else
                PrintErrorMessage(system, "Output period (dt or step_jump must be defined > 0")
                return c_error
            end
        else
            o_block.dt = o_block.step_jump * system.dt
        end

        if o_block.dt < system.dt
            PrintWarningMessage(system, "Output period has been set equal to simulation time step")
            o_block.step_jump = 1
            o_block.dt = system.dt
        end


        # Check output start and end times
        if o_block.step_end <= o_block.step_start
            if o_block.t_start < o_block.t_end
                o_block.step_start = round(Int64, o_block.t_start / system.dt)
                o_block.step_end = round(Int64, o_block.t_end / system.dt)
            else
                PrintErrorMessage(system, "Output time ranges are badly defined")
                return c_error
            end
        else
            o_block.t_start = o_block.step_start * system.dt
            o_block.t_end = o_block.step_end * system.dt
        end

        # Check averaged outputs
        if o_block.averaged

            # Check average time range
            if o_block.step_av <= 0
                if o_block.dt_av > 0.0
                    o_block.step_av = round(Int64, o_block.dt_av / system.dt)
                else
                    PrintErrorMessage(system, "Output average time must be defined positive")
                    return c_error
                end
            else
                o_block.dt_av = o_block.step_av * system.dt
            end

            # Set average time parameters
            o_block.step_av_start = o_block.step_start
            o_block.step_av_end = o_block.step_av_start + o_block.step_av
            o_block.step_dump = o_block.step_av_end
            o_block.time_av_start = round(Int64, o_block.step_av_start * system.dt)
            o_block.time_av_end = round(Int64, o_block.step_av_end * system.dt)
            o_block.time_dump = round(Int64, o_block.step_dump * system.dt)
        else
            o_block.step_dump = o_block.step_start
        end

        errcode = 0
    end

    return errcode
end
    
function EndFile_Output!(read_step::Int64, output_list::Vector{OutputBlock},
    system::System)
    errcode = 0
    return errcode
end

function Initialize_OutputDataStruct!(data_struct::OutputDataStruct,
    name::Union{String, SubString{String}}, p_id::Int64, var_list::Vector{String}, o_block::OutputBlock,
    species_list::Vector{Species}, system::System)

    # Identify the output parameter by its ID and name
    data_struct.id = p_id
    data_struct.name = name

    # Check wether parameter has a direction 
    dir_id = GetSpatialDirectionID(var_list)
    data_struct.dir_id = dir_id

    # Check wether parameter is species related
    s_id = GetSpeciesID(var_list, species_list)
    data_struct.species_id = s_id
    if s_id == c_o_all_species
        data_struct.species_name = "all" 
    elseif s_id == c_o_none_species
        data_struct.species_name = "N/A"
    elseif s_id > 0
        data_struct.species_name = species_list[s_id].name 
    end

    # Data array is only used if
    # - averaged data need to be buffed
    # - particle probe data is buffed
    if (o_block.averaged) &&
        (data_struct.parameter_id != c_o_phase_space) &&
        (data_struct.parameter_id != c_o_probe)
        # particle phase-space is the only parameter which is not averaged

        # The length of the data is given by the number of cells
        ncells = system.ncells

        # Depending on whether one species (or none) or many species are buffered
        if data_struct.species_id == c_o_all_species
            n_species = 0
            for s in species_list
                if !s.is_background_species
                    n_species += 1
                end
            end
            data_struct.data = zeros(Float64, ncells, n_species) 
        else
            data_struct.data = zeros(Float64, ncells, 1) 
        end
    end

end

function SplitVariableString(var::Union{String, SubString{String}})

    var = strip(var)
    var_list = String[]
    if length(var) == 0
        return var_list
    end

    comma = ","
    find_comma = true
    while find_comma
        index = findfirst(comma, var)
        if index === nothing
            find_comma = false
            push!(var_list, strip(var))
        else
            ix_end = index[1]-1
            var_item = var[1:ix_end]
            push!(var_list, strip(var_item))

            var = var[ix_end+2:end]
        end
    end
    return var_list
end

function GetSpeciesID(var_list::Union{Vector{String}, Vector{SubString{String}}},
    species_list::Vector{Species})

    s_id = c_o_none_species

    for var in var_list
        for s in species_list
            if s.name == var
                s_id = s.id
            end
        end
        if (var == "all") || (var == "allspecies")
            s_id = c_o_all_species
        end
    end
    return s_id 
end

function GetSpatialDirectionID(var_list::Union{Vector{String}, Vector{SubString{String}}})
    
    dir_id = c_dir_none

    for var in var_list
        if var == "x"
            dir_id = c_dir_x
            return dir_id
        elseif var == "y"
            dir_id = c_dir_y
            return dir_id
        elseif var == "z"
            dir_id = c_dir_z
            return dir_id
        end
    end

    return dir_id
end

end