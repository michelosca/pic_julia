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

module Inputs

using Constants: c_error
using Constants: c_block_system, c_block_species, c_block_output, c_block_mcc
using Constants: c_block_constants, c_block_waveform
using SharedData: Species, System, OutputBlock, Waveform
using SharedData: Collision, CollisionGroup
using InputBlock_System: StartFile_System!, StartSystemBlock!
using InputBlock_System: ReadSystemEntry!, EndSystemBlock!, EndFile_System!
using InputBlock_Species: StartFile_Species!, StartSpeciesBlock!
using InputBlock_Species: ReadSpeciesEntry!, EndSpeciesBlock!, EndFile_Species!
using InputBlock_Output: StartFile_Output!, StartOutputBlock!
using InputBlock_Output: ReadOutputEntry!, EndOutputBlock!, EndFile_Output!
using InputBlock_Waveform: StartFile_Waveform!, StartWaveformBlock!
using InputBlock_Waveform: ReadWaveformEntry!, EndWaveformBlock!, EndFile_Waveform!
using InputBlock_MCC: StartFile_MCC!, StartMCCBlock!
using InputBlock_MCC: ReadMCCEntry!, EndMCCBlock!, EndFile_MCC!

# The input deck is read twice
# 1ST READ
## -> Read SYSTEM paramters: ncells, dx, x_min, x_max, etc.
##    - System parameters are checked after block is read so that these values
##    should be ready to use at the end of the first read of the input deck,
##    i.e. at EndFile_<block> at 1st read system parameters can be used
##
## -> Read SPECIES blocks and parameters: names, id, mass, part.count, etc.
##    - Species parameters are check at  EndFile_Species, therefore species
##    data should be fully available for the second read
## -> Read ELECTRIC/MAGNETIC FIELDS and ELECTRIC POTENTIAL.
##    - The data is not spatially allocated (this could be done at 2nd read) but parameters
##      describing the fields, e.g. magnitude values, etc. are ready to use at 2nd read

# 2ND READ
## -> TIME STEP CONDITIONS at SYSTEM are set at StartFile_System
## -> Read OUTPUT parameters: when system and species data is already available
##    - All output parameters are defined at the end of each block
## -> Read WAVEFORM parameters
##    - Waveform functions for potential
## -> Read COLLISION parameters
##    - Read and setup MCC

# INPUT BLOCK IDs
global block_id = 0

function SetupInputData!(filename::String
    , species_list::Vector{Species}
    , system::System
    , output_list::Vector{OutputBlock}
    , waveform_list::Vector{Waveform}
    , collisiongroup_list::Vector{CollisionGroup}
    )

    errcode = 0

    collision_list = Collision[]

    # Opens the file given in filename and reads each line
    for read_step in 1:2
        print("Reading $read_step of the input deck...\n")
        errcode = StartFile!(read_step, species_list, system, output_list,
            waveform_list, collision_list)
        if (errcode == c_error) return errcode end

        errcode = ReadFile!(filename, read_step, species_list, system,
            output_list, waveform_list, collision_list)
        if (errcode == c_error) return errcode end

        errcode = EndFile!(read_step, species_list, system, output_list,
            waveform_list, collision_list, collisiongroup_list)
        if (errcode == c_error) return errcode end
        print("End of input deck reading\n\n")
    end
    GetInputFolder!(system, filename)
    
    return errcode
end


function ReadFile!(filename::String, read_step::Int64,
    species_list::Vector{Species},
    system::System, output_list::Vector{OutputBlock},
    waveform_list::Vector{Waveform},
    collision_list::Vector{Collision}
    )
    
    errcode = 0

    open(filename,"r") do f
        constants = Tuple{SubString{String}, SubString{String}}[]
        line = 1
        while ! eof(f)
            s = readline(f)
            # ReadLine identifies each line on filename
            errcode = ReadLine!(s, read_step, species_list,
                system, output_list, constants, waveform_list,
                collision_list)
            if (errcode == c_error)
                print("***ERROR*** Stop reading at file line ", line,"\n")
                return errcode  
            end
            line += 1
        end
    end
    return errcode
end


function ReadLine!(str::String, read_step::Int64,
    species_list::Vector{Species},
    system::System, output_list::Vector{OutputBlock},
    constants::Vector{Tuple{SubString{String},SubString{String}}},
    waveform_list::Vector{Waveform},
    collision_list::Vector{Collision}
    )

    errcode = c_error 

    # Trimm out commets
    i_comment = findfirst("#", str)
    if !(i_comment === nothing)
        i_comment = i_comment[1]
        str = str[begin:i_comment-1]
    end
    str = strip(str)

    # Signal wether this is a begin/end:block 
    i_block = findfirst(":", str)
    # Signla whether it is a "name = var" line
    i_eq = findfirst("=", str)

    # Check line
    if !(i_block === nothing)
        i_block = i_block[1]
        global block_name = str[i_block+1:end]
        if (occursin("begin", str))
            errcode = StartBlock!(block_name, read_step, species_list,
                system, output_list, waveform_list, collision_list)
            if (errcode == c_error)
                print("***ERROR*** Something went wrong starting the ",
                    block_name, " block\n")
                return c_error
            end
        elseif (occursin("end", str))
            errcode = EndBlock!(block_name, read_step, species_list,
                system, output_list, waveform_list, collision_list)
            if (errcode == c_error)
                print("***ERROR*** Something went wrong ending the ",
                    block_name, " block\n")
                return c_error
            end
        end
    elseif !(i_eq === nothing)
        i_eq = i_eq[1]
        name = strip(str[begin:i_eq-1])
        var = strip(str[i_eq+1:end])
        errcode = ReadInputDeckEntry!(name, var, read_step,
            species_list, system, output_list, constants, waveform_list,
            collision_list)
        if (errcode == c_error)
            print("***WARNING*** Entry in ", block_name,
                "-block has not been located\n")
            print("  - Input entry: ", name ," = ",var ,"\n")
        end
    else
        errcode = 0
    end
    return errcode 
end


function StartFile!(read_step::Int64, species_list::Vector{Species},
    system::System, output_list::Vector{OutputBlock},
    waveform_list::Vector{Waveform},
    collision_list::Vector{Collision}
    )

    errcode = StartFile_Species!(read_step, species_list, system) 
    if (errcode == c_error)
        print("***ERROR*** While initializing the input species block")
    end

    errcode = StartFile_System!(read_step, species_list, system) 
    if (errcode == c_error)
        print("***ERROR*** While initializing the input system block")
    end
    
    errcode = StartFile_Output!(read_step, output_list, system) 
    if (errcode == c_error)
        print("***ERROR*** While initializing the input output block")
    end

    errcode = StartFile_Waveform!(read_step, waveform_list, system) 
    if (errcode == c_error)
        print("***ERROR*** While initializing the input waveform block")
    end

    errcode = StartFile_MCC!(read_step, collision_list, system) 
    if (errcode == c_error)
        print("***ERROR*** While initializing the input MCC block")
    end
    
    return errcode
end


function EndFile!(read_step::Int64, species_list::Vector{Species},
    system::System, output_list::Vector{OutputBlock},
    waveform_list::Vector{Waveform},
    collision_list::Vector{Collision},
    collisiongroup_list::Vector{CollisionGroup}
    )

    errcode = EndFile_Species!(read_step, species_list, system)
    if (errcode == c_error)
        print("***ERROR*** While finalizing the input species block\n")
        return errcode
    end

    errcode = EndFile_System!(read_step, system) 
    if (errcode == c_error)
        print("***ERROR*** While finalizing the input system block\n")
        return errcode
    end
    
    errcode = EndFile_Output!(read_step, output_list, system) 
    if (errcode == c_error)
        print("***ERROR*** While finalizing the input output block\n")
        return errcode
    end
    
    errcode = EndFile_Waveform!(read_step, waveform_list, system) 
    if (errcode == c_error)
        print("***ERROR*** While finalizing the input waveform block\n")
        return errcode
    end
    
    errcode = EndFile_MCC!(read_step, collisiongroup_list,
        collision_list, species_list, system) 
    if (errcode == c_error)
        print("***ERROR*** While finalizing the input MCC block\n")
        return errcode
    end
    
    return errcode
end


function StartBlock!(name::SubString{String}, read_step::Int64,
    species_list::Vector{Species},
    system::System, output_list::Vector{OutputBlock},
    waveform_list::Vector{Waveform},
    collision_list::Vector{Collision}
    )

    errcode = c_error
    if (occursin("system",name))
        global block_id = c_block_system
        errcode = StartSystemBlock!(read_step, system)
    elseif (occursin("species",name))
        global block_id = c_block_species
        errcode = StartSpeciesBlock!(read_step, species_list)
    elseif (occursin("output",name))
        global block_id = c_block_output
        errcode = StartOutputBlock!(read_step, output_list, system)
    elseif (occursin("waveform",name))
        global block_id = c_block_waveform
        errcode = StartWaveformBlock!(read_step, waveform_list, system)
    elseif (occursin("MCC",name))
        global block_id = c_block_mcc
        errcode = StartMCCBlock!(read_step, collision_list, species_list,
            system)
    elseif (occursin("constants",name))
        global block_id = c_block_constants
        errcode = 0
    end
    return errcode
end


function EndBlock!(name::SubString{String}, read_step::Int64,
    species_list::Vector{Species}, system::System,
    output_list::Vector{OutputBlock}, waveform_list::Vector{Waveform},
    collision_list::Vector{Collision}
    )

    errcode = c_error
    global block_id = 0
    if (occursin("system",name))
        errcode = EndSystemBlock!(read_step, system)
    elseif (occursin("species",name))
        errcode = EndSpeciesBlock!(read_step, species_list, system)
    elseif (occursin("output",name))
        errcode = EndOutputBlock!(read_step, output_list, system)
    elseif (occursin("waveform",name))
        errcode = EndWaveformBlock!(read_step, waveform_list, system)
    elseif (occursin("MCC",name))
        errcode = EndMCCBlock!(read_step, collision_list, system)
    elseif (occursin("constants",name))
        errcode = 0
    end
    return errcode
end


function ReadInputDeckEntry!(name::SubString{String}, var::SubString{String},
    read_step::Int64, species_list::Vector{Species}, system::System,
    output_list::Vector{OutputBlock},
    constants::Vector{Tuple{SubString{String},SubString{String}}},
    waveform_list::Vector{Waveform},
    collision_list::Vector{Collision})

    errcode = c_error 

    var = CheckConstantValues!(var, constants)
    
    if (block_id == c_block_system)
        errcode = ReadSystemEntry!(name, var, read_step, system)
    elseif (block_id == c_block_species)
        errcode = ReadSpeciesEntry!(name, var, read_step, species_list, system)
    elseif (block_id == c_block_output)
        errcode = ReadOutputEntry!(name, var, read_step, output_list,
            species_list, system)
    elseif (block_id == c_block_waveform)
        errcode = ReadWaveformEntry!(name, var, read_step, waveform_list)
    elseif (block_id == c_block_mcc)
        errcode = ReadMCCEntry!(name, var, read_step, collision_list)
    elseif (block_id == c_block_constants)
        push!(constants, (name,var))
        errcode = 0
    end

    return errcode 
end


function CheckConstantValues!(var::SubString{String},
    constants::Vector{Tuple{SubString{String},SubString{String}}})

    for c in constants
        if var == c[1]
            var = c[2]
            break
        end
    end

    return var
end


function GetInputFolder!(system::System, filename::String)

    index = findlast("/", filename)
    if index === nothing
        system.folder = "./"
    else
        system.folder = filename[1:index[1]]
    end
    system.log_file = system.folder * system.log_file
end

end