module Outputs

using HDF5
using Printf: @sprintf
using Constants: c_error, c_o_all_species, c_o_none_species
using Constants: c_o_density, c_o_potential, c_o_electric_field, c_o_probe
using Constants: c_o_phase_space, c_o_neutral_collisions
using Constants: c_dir_x, c_dir_y, c_dir_z
using Constants: c_bc_open, c_bc_periodic
using SharedData: System, Field, Species, CollisionGroup
using SharedData: OutputBlock, OutputDataStruct
using PrintModule: PrintErrorMessage

function GenerateOutputs!(output_list::Vector{OutputBlock},
    species_list::Vector{Species}, system::System, electric_potential::Vector{Float64},
    electric_field::Field, collgroup_list::Vector{CollisionGroup})

    step = system.step 
    for o_block in output_list
        if o_block.averaged
            av_start = o_block.step_av_start
            av_end = o_block.step_av_end
            buffer_on = (step >= av_start) && (step <= av_end) 
            average_flag = step == av_end

            if buffer_on 
                BuffAveragedData!(o_block, system, species_list, electric_potential,
                electric_field, collgroup_list, average_flag)
            end
        end

        dump_step = o_block.step_dump
        if step >= dump_step

            DumpOutputData(o_block, system, species_list, electric_potential, electric_field, collgroup_list)

            # Update o_block parameters
            o_block.file_id += 1
            o_block.step_dump += o_block.step_jump
            o_block.time_dump += o_block.dt
            if o_block.averaged
                o_block.step_av_start += o_block.step_jump
                o_block.step_av_end   = o_block.step_av_start + o_block.step_av
                o_block.time_av_start += o_block.dt
                o_block.time_av_end   = o_block.time_av_start + o_block.dt_av
            end

        end

    end # Loop over OutputBlocks in output_list

end

function BuffAveragedData!(o_block::OutputBlock, system::System,
    species_list::Vector{Species}, electric_potential::Vector{Float64},
    electric_field::Field, collgroup_list::Vector{CollisionGroup}, average_flag::Bool)

    if average_flag
        step_range = Float64(o_block.step_av_end - o_block.step_av_start+1)
    end

    for param in o_block.param_list

        # DENSITY data
        if param.id == c_o_density
            c_min = system.cell_min
            c_max = system.cell_max

            # Cases: either average all species or just one single species
            n_species = 1
            for s in species_list
                if s.is_background_species
                    continue
                elseif param.species_id == c_o_all_species ||
                    param.species_id == s.id
                    param.data[:, n_species] += s.dens[c_min:c_max]
                    
                    n_species += 1
                end
            end

            # If last buffer step, average over the n steps
            if average_flag
                param.data ./= step_range
            end

        end
    end

end

function DumpOutputData(o_block::OutputBlock, system::System,
    species_list::Vector{Species}, electric_potential::Vector{Float64},
    electric_field::Field, collgroup_list::Vector{CollisionGroup})

    filename = @sprintf("%s/%s%05i.h5", system.folder, o_block.name, o_block.file_id)
    h5open(filename,"w") do fid

        # SYSTEM parameters
        errcode = WriteSystemDataToH5!(fid, system)
        if errcode == c_error
            err_message = @sprintf("While writing system data in output block %s", o_block.name)
            PrintErrorMessage(system, err_message)
            return errcode
        end

        for param in o_block.param_list

            # DENSITY data
            if param.id == c_o_density
                errcode = WriteDensityToH5!(fid, species_list, param, system, o_block.averaged)
                if errcode == c_error
                    err_message = @sprintf("While writing density data in output block %s", o_block.name)
                    PrintErrorMessage(system, err_message)
                    return errcode
                end
            end

            # ELECTRIC POTENTIAL data
            if param.id == c_o_potential
                errcode = WriteElectricPotentialToH5!(fid, electric_potential, param, system)
                if errcode == c_error
                    err_message = @sprintf("While writing electric potential data in output block %s", o_block.name)
                    PrintErrorMessage(system, err_message)
                    return errcode
                end
            end

            # ELECTRIC field data
            if param.id == c_o_electric_field
                errcode = WriteElectricFieldToH5!(fid, electric_field, param, system)
                if errcode == c_error
                    err_message = @sprintf("While writing electric field data in output block %s", o_block.name)
                    PrintErrorMessage(system, err_message)
                    return errcode
                end
            end

            # Phase-Space data
            if param.id == c_o_phase_space
                for species in species_list
                    if species.is_background_species
                        continue
                    elseif (param.species_id == species.id) ||
                        (param.species_id == c_o_all_species)
                        errcode = WritePhaseSpaceToH5!(fid, species, param)
                        if errcode == c_error
                            err_message = @sprintf("While writing %s phase-space data in output block %s",
                                species.name, o_block.name)
                            PrintErrorMessage(system, err_message)
                            return errcode
                        end
                    end
                end
            end

            # COLLISIONS data
            if param.id == c_o_neutral_collisions
                errcode = WriteNeutralCollisionsToH5!(fid, collgroup_list, param, system)
                if errcode == c_error
                    err_message = @sprintf("While writing neutral collisions in output block %s", o_block.name)
                    PrintErrorMessage(system, err_message)
                    return errcode
                end
            end


        end # Loop over OutputBlockStruct in current OutputBlock:o_block
    end # HDF5 file is closed
end

function WriteDensityToH5!(fid::HDF5.File, species_list::Vector{Species},
    param::OutputDataStruct, system::System, averaged::Bool)

    errcode = c_error

    ncells = system.ncells
    c_min = system.cell_min
    c_max = system.cell_max
    g = create_group(fid, "Number_Density")
    n_species = 1
    for s in species_list
        if s.is_background_species
            continue
        elseif (s.id == param.species_id) ||
            (param.species_id == c_o_all_species)
            dset = create_dataset(g, s.name, Float64,(ncells,))

            if averaged
                write(dset, param.data[:,n_species])
                param.data[:,n_species] .= s.dens[c_min:c_max]
                n_species += 1
            else
                write(dset, s.dens[c_min:c_max])
            end
            errcode = 0
        end
    end

    return errcode
end

function WriteSystemDataToH5!(fid::HDF5.File, system::System)

    errcode = c_error

    ncells = system.ncells
    x_min = system.x_min
    x_max = system.x_max

    # Create System group
    g = create_group(fid,"System")

    # Write spatial grid layout
    dset = create_dataset(g,"Grid",Float64,(ncells,))
    x = map(x -> x, LinRange(x_min,x_max,ncells) )
    write(dset, x)
        
    # Write system attributes
    attributes(g)["ncells"] = ncells 
    attributes(g)["dx"] = system.dx 
    attributes(g)["dt"] = system.dt 
    attributes(g)["time"] = system.time
    attributes(g)["step"] = system.step

    if system.bc_field_max == c_bc_open
        bc_str = "open"
    elseif system.bc_field_max == c_bc_periodic
        bc_str = "periodic"
    end
    attributes(g)["bc_field_max"] = bc_str 

    if system.bc_field_min == c_bc_open
        bc_str = "open"
    elseif system.bc_field_min == c_bc_periodic
        bc_str = "periodic"
    end
    attributes(g)["bc_field_min"] = bc_str 

    if system.bc_part_max == c_bc_open
        bc_str = "open"
    elseif system.bc_part_max == c_bc_periodic
        bc_str = "periodic"
    end
    attributes(g)["bc_part_max"] = bc_str 

    if system.bc_part_min == c_bc_open
        bc_str = "open"
    elseif system.bc_part_min == c_bc_periodic
        bc_str = "periodic"
    end
    attributes(g)["bc_part_min"] = bc_str

    errcode = 0

    return errcode
end

function WriteElectricPotentialToH5!(fid::HDF5.File, pot::Vector{Float64},
    param::OutputDataStruct, system::System)

    errcode = c_error

    ncells = system.ncells
    c_min = system.cell_min
    c_max = system.cell_max
    g = create_group(fid, "Electric_Potential")
    if param.dir_id == c_dir_x
        dset = create_dataset(g, "Vx", Float64,(ncells,))
        write(dset, pot[c_min:c_max])
        errcode = 0
    end


    return errcode
end

function WriteElectricFieldToH5!(fid::HDF5.File, efield::Field,
    param::OutputDataStruct, system::System)

    errcode = c_error

    ncells = system.ncells
    c_min = system.cell_min
    c_max = system.cell_max

    # Open E-field group
    g_name = "Electric_Field"
    if haskey(fid, g_name)
        g = fid[g_name]
    else
        g = create_group(fid, g_name)
    end

    # Write data
    if param.dir_id == c_dir_x
        dset = create_dataset(g, "Ex", Float64,(ncells,))
        write(dset, efield.x[c_min:c_max])
        errcode = 0
    elseif param.dir_id == c_dir_y
        dset = create_dataset(g, "Ey", Float64,(ncells,))
        write(dset, efield.y[c_min:c_max])
        errcode = 0
    elseif param.dir_id == c_dir_z
        dset = create_dataset(g, "Ez", Float64,(ncells,))
        write(dset, efield.z[c_min:c_max])
        errcode = 0
    end

    return errcode
end

function WritePhaseSpaceToH5!(fid::HDF5.File, species::Species,
    param::OutputDataStruct)

    errcode = c_error

    # Open particle group
    g_name = "Particles"
    if haskey(fid, g_name)
        g = fid[g_name]
    else
        g = create_group(fid, g_name)
    end
    
    # Open SPECIES group
    if haskey(g, species.name)
        g_species = g[species.name]
    else
        g_species = create_group(g, species.name)
        errcode = WriteParticleAttributesToH5!(g_species, species)
        if errcode == c_error
            err_message = @sprintf("While writing %s attributes data", species.name)
            PrintErrorMessage(system, err_message)
            return errcode
        end
    end

    # Number of particles to be writen
    nparts = species.particle_count

    # Write particles position 
    pos_label = "x"
    if !haskey(g_species, pos_label)
        dset = create_dataset(g_species, pos_label, datatype(Float64), dataspace(nparts,1))
        pos_x_array = map(x -> x.pos, species.particle_list) 
        write(dset, pos_x_array)
    end

    # Write particle velocity
    if param.dir_id == c_dir_x
        dset = create_dataset(g_species, "vx", datatype(Float64), dataspace(nparts,1))
        v_array = map(x -> x.vel[1], species.particle_list) 
        write(dset, v_array)
        errcode = 0
    elseif param.dir_id == c_dir_y
        dset = create_dataset(g_species, "vy", datatype(Float64), dataspace(nparts, 1))
        v_array = map(x -> x.vel[2], species.particle_list) 
        write(dset, v_array)
        errcode = 0
    elseif param.dir_id == c_dir_z
        dset = create_dataset(g_species, "vz", datatype(Float64), dataspace(nparts, 1))
        v_array = map(x -> x.vel[3], species.particle_list) 
        write(dset, v_array)
        errcode = 0
    end

    return errcode
end

function WriteNeutralCollisionsToH5!(fid::HDF5.File, collgroup_list::Vector{CollisionGroup},
    param::OutputDataStruct, system::System)

    errcode = c_error

    ncells = system.ncells - 1

    g_name = "Neutral_Collisions"
    if haskey(fid, g_name)
        g = fid[g_name]
    else
        g = create_group(fid, "Neutral_Collisions")
        # Add specific grid because NC use a cell less
        dset = create_dataset(g, "NC_Grid", Float64,(ncells,))
        dx = system.dx
        x_min = system.x_min + 0.5*dx
        x_max = system.x_max - 0.5*dx
        x = map(x -> x, LinRange(x_min,x_max,ncells) )
        write(dset, x )
    end
    
    # Load NC to H5 file
    for cgroup in collgroup_list
        species_combination = ""
        for (i,s) in enumerate(cgroup.colliding_species)
            if i == length(cgroup.colliding_species)
                species_combination *= s.name
            else
                species_combination *= s.name * ":"
            end
        end
        g1 = create_group(g, species_combination)

        for coll in cgroup.collision_list
            dset = create_dataset(g1, coll.name, Float64,(ncells,))

            write(dset, coll.diagnostic)
            errcode = 0
        end
    end

    errcode = 0
    return errcode
end

function WriteParticleAttributesToH5!(g::HDF5.Group, s::Species)

    errcode = c_error

    attributes(g)["weight"] = s.weight 
    attributes(g)["mass"] = s.mass 
    attributes(g)["charge"] = s.charge 
    attributes(g)["part_count"] = s.particle_count 
    attributes(g)["background"] = s.is_background_species

    errcode = 0
    return errcode
end

end