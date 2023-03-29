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

module InputBlock_MCC

using Constants: c_error
using Constants: c_bc_open
using PrintModule: PrintMessage, PrintErrorMessage, PrintCollision
using SharedData: System, Species, Collision
using Printf
using Tools: GetUnits!


function StartFile_MCC!(read_step::Int64, collision_list::Vector{Collision},
    system::System) 

    errcode = c_error
    
    if read_step == 1
        errcode = 0
    elseif read_step == 2
        errcode = 0
    end

    return errcode
end


function StartMCCBlock!(read_step::Int64, collision_list::Vector{Collision},
    species_list::Vector{Species}, system::System)

    errcode = c_error

    if (read_step == 1)
        errcode = 0
    elseif read_step == 2
        system.mcc = true
        errcode = LoadCollisions!(collision_list, species_list, system)
        if (errcode == c_error)
            message = @sprintf("While loading collisions")
            PrintErrorMessage(system, message)
            return errcode  
        end

        for c in collision_list
            PrintCollision(c)
            print("\n\n")
        end

        errcode = 0
    end
    return errcode
end


function ReadMCCEntry!(name::SubString{String}, var::SubString{String},
    read_step::Int64, collision_list::Vector{Collision})

    errcode = c_error 

    if (read_step == 1)
        errcode = 0
    elseif read_step == 2
        units_fact, name = GetUnits!(name)
        errcode = 0
    end
    return errcode 
end


function EndMCCBlock!(read_step::Int64, collision_list::Vector{Collision},
    system::System)

    errcode = 0

    return errcode
end
    

function EndFile_MCC!(read_step::Int64, collision_list::Vector{Collision},
    system::System) 

    errcode = 0 

    return errcode
end


function LoadCollisions!(collision_list::Vector{Collision}
    , species_list::Vector{Species}
    , system::System)

    errcode = 0

    main = "./src/collisions/"
    for species in species_list

        path = main * species.name
        table_list = ""
        try
            table_list = readdir(path)
        catch
            message = "No collisions for " * species.name
            PrintMessage(system, message)
            continue
        end

        for table in table_list
            
            # Only files with '.table' extension are evaluated
            ix = findlast('.', table)
            if table[ix+1:end] != "table"
                continue
            end

            # For each table generate a Collision structure and load table data into it
            open(path * "/" * table) do file
                line = 1
                while ! eof(file)
                    s = readline(file)
                    # ReadLine identifies each line on filename
                    errcode = ReadTableLine!(s, collision_list, species_list, system)
                    if (errcode == c_error)
                        message = @sprintf("Stop reading at file %s (species %s) at line %i",
                            table, species.name, line)
                        PrintErrorMessage(system, message)
                        return errcode  
                    end
                    line += 1
                end
            end
        end
    end

    return errcode
end

function ReadTableLine!(str::String, collision_list::Vector{Collision},
    species_list::Vector{Species}, system::System)

    errcode = 0

    # Name description
    i_name = findfirst(':', str)
    if i_name === nothing
        return errcode
    end

    # Get name string
    str_name = strip(str[1:i_name-1])

    # Get variable
    var_str = str[i_name + 1:end] 
    var_str = strip(var_str)

    if str_name == "SPECIES"
        # This entry is necessary in order to create a new collision structure
        collision = InitCollision()
        push!(collision_list, collision)
        errcode = ParseSpecies!(var_str, collision_list[end], species_list, system)
    elseif str_name == "PROCESS"
        errcode = ParseProcess!(var_str, collision_list[end], species_list, system)
    elseif str_name == "PARAM."
        errcode = ParseParam!(var_str, collision_list[end], species_list, system)
    #elseif str_name == "COLUMNS"
    end
    return errcode
end

function ParseSpecies!(var::Union{String, SubString{String}},
    collision::Collision, species_list::Vector{Species},
    system::System)

    errcode = 0
    
    species_name_list = SplitString(var,'/')
    for species_name in species_name_list
        s, errcode = IdentifySpecies(species_name, species_list)

        # In case species_name is not found then send an error message
        if errcode == c_error
            message = @sprintf("Species %s in collision table must be declared in input deck", species_name)
            PrintErrorMessage(system, message)
            return errcode
        else
            is_species_already_included = false
            for species in collision.species 
                if species.name == s.name
                    is_species_already_included = true 
                    break
                end
            end

            if !is_species_already_included
                push!(collision.species, s)
            end
        end
        
    end
    
    return errcode
end

function ParseProcess!(var::Union{String, SubString{String}},
    collision::Collision, species_list::Vector{Species},
    system::System)

    errcode = 0
    
    # Split collision reaction from collision type (separator ',')
    collision_reaction_and_type = SplitString(var, ',')
    coll_reaction_str = collision_reaction_and_type[1]

    ### Collision type (name)
    coll_type_str = collision_reaction_and_type[2]
    collision.name = coll_type_str

    # Find two reaction parts separated by "->" 
    reactant_and_product = SplitString(coll_reaction_str, "->")
    if length(reactant_and_product) < 2 
        message = @sprintf("Reaction equation %s has could not be parsed", var)
        PrintErrorMessage(system, message)
        return c_error
    end

    ### Reaction reactants 
    reactant_str = reactant_and_product[1] 
    reactant_species_names = SplitString(reactant_str, " + ")
    for species_name in reactant_species_names
        s, errcode = IdentifySpecies(species_name, species_list)
        errcode = ParseSpecies!(species_name, collision, species_list, system)

        # In case species_name is not found then send an error message
        if errcode == c_error
            message = @sprintf("Reactant %s in collision %s must be declared in input deck", species_name, collision.name)
            PrintErrorMessage(system, message)
            return errcode
        else
            push!(collision.reactants, s)
        end
        
    end

    ### Reaction products
    product_str = reactant_and_product[2] 
    product_species_names = SplitString(product_str, " + ")
    for species_name in product_species_names
        s, errcode = IdentifySpecies(species_name, species_list)
        errcode = ParseSpecies!(species_name, collision, species_list, system)

        # In case species_name is not found then send an error message
        if errcode == c_error
            message = @sprintf("Product %s in collision %s must be declared in input deck", species_name, collision.name)
            PrintErrorMessage(system, message)
            return errcode
        else
            push!(collision.products, s)

            ## In case species is not yet in collision.species list -> include it
            #check_species_list = true 
            #for cs in collision.species
            #    if cs.name == s.name
            #        check_species_list = false 
            #        break
            #    end
            #end
            #if check_species_list
            #    push!(collision.species, s)
            #end
        end
    end

    ### Reaction species balance
    for s in collision.species
        balance = 0

        # Loop over reactant species
        for sr in collision.reactants
            if sr.name == s.name
                balance -= 1
            end
        end
        # Loop over product species
        for sp in collision.products
            if sp.name == s.name
                balance += 1
            end
        end

        push!(collision.species_balance, balance)
    end

    return errcode
end

function ParseParam!(var::Union{String, SubString{String}},
    collision::Collision, species_list::Vector{Species},
    system::System)

    errcode = 0
    
    # Split collision reaction from collision type (separator ',')
    e_threshold_and_others = SplitString(var, ',')
    e_threshold_str = e_threshold_and_others[1]

    ### Collision energy threshold 
    name_eq_value = SplitString(e_threshold_str, '=')
    name = name_eq_value[1]
    if name != "E"
        return errcode
    end
    e_value = name_eq_value[2]
    units_fact, val_str = GetUnits!(e_value, symb=' ')


    collision.energy_threshold = parse(Float64, val_str) * units_fact 

    return errcode
end

function IdentifySpecies(species_name::Union{String, SubString{String}},
    species_list::Vector{Species})

    errcode = c_error
    #If species is not identify the flag an error
    for s in species_list

        # Match species found with declared species
        if lowercase(s.name) == lowercase(species_name)
            errcode = 0
            return s, errcode
        end
    end
    return nothing, errcode
end

function SplitString(var::Union{String, SubString{String}},
    symb::Union{Char, String})
    
    species_name_list = String[]

    # Find any species separated by 'symb'
    ix = findfirst(symb, var)
    while ix !== nothing

        if typeof(symb) == String
            # Identify species
            species_name = strip(var[1:ix[1]-1])
            push!(species_name_list, species_name)

            # Trim the line and look for the next species
            var = var[ix[2]+1:end]
            ix = findfirst(symb, var)

        else
            # Identify species
            species_name = strip(var[1:ix-1])
            push!(species_name_list, species_name)

            # Trim the line and look for the next species
            var = var[ix+1:end]
            ix = findfirst(symb, var)
        end
    end
    
    # Identify last species in line
    species_name = strip(var)
    push!(species_name_list, species_name)

    return species_name_list 
end

function InitCollision()
    collision = Collision()
    collision.name = "collision_name"
    collision.reactants = Species[]
    collision.products = Species[]
    collision.species = Species[]
    collision.species_balance= Int64[]
    collision.energy_threshold = 0.0
    collision.energy_data = Float64[]
    collision.cross_section_data = Float64[]
    collision.diagnostic = Float64[]
    return collision
end

end