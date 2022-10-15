module Outputs

using HDF5
using Printf: @sprintf
using SharedData: System, Field, Species

function WriteDataToHDF5(species_list::Vector{Species}, system::System,
    pot::Vector{Float64}, Efield::Field)

    ncells = system.ncells
    c_min = system.cell_min
    c_max = system.cell_max

    filename = @sprintf("%s/%05i.h5",system.folder, system.step)
    h5open(filename,"w") do fid  #preserving existing "cw"

        g = create_group(fid,"System")
        dset = create_dataset(g,"Grid",Float64,(ncells,))
        x_min = system.x_min
        x_max = system.x_max
        x = map(x -> x, LinRange(x_min,x_max,ncells) )
        write(dset, x)
        attributes(g)["ncells"] = ncells 
        attributes(g)["dx"] = system.dx 
        attributes(g)["dt"] = system.dt 
        attributes(g)["time"] = system.time

        g = create_group(fid,"Electric_Potential")
        dset = create_dataset(g,"Vx",Float64,(ncells,))
        write(dset, pot[c_min:c_max])

        g = create_group(fid,"Electric_Field")
        dset = create_dataset(g,"Ex",Float64,(ncells,))
        write(dset, Efield.x[c_min:c_max])

        g = create_group(fid,"Number_Density")
        for s in species_list
            if !s.is_background_species
                dset = create_dataset(g, s.name, Float64, (ncells,))
                write(dset, s.dens[c_min:c_max])
            end
        end
    end

end

end