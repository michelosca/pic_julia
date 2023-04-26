# Copyright (C) 2023 Michel Osca Engelbrecht
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
# along with GM Julia. If not, see <https://www.gnu.org/licenses/>.

module NeutralCollisions

using SharedData: System, Species, CollisionGroup, Particle
using Constants: c_error, kb
using EvaluateExpressions: ReplaceExpressionValues
using Tools: linear_interpolation
using PrintModule: PrintErrorMessage

using Random: rand, randn, shuffle!
using LinearAlgebra: norm

function NeutralCollisions!(coll_list::Vector{CollisionGroup}, species_list::Vector{Species}, system::System)

    dt = system.dt
    dx = system.dx
    idx = 1.0/dx

    for coll_group in coll_list

        gsigma_max = coll_group.gsigma_max
        collision_list = coll_group.collision_list
        iweight_max = 1.0 / coll_group.part_weight_max

#        print(coll_group.colliding_species[1].name," - ")
#        print(coll_group.colliding_species[2].name,"\n")

        # Loop over ncells-1 because particle_grid_list is ncells-1 long
        for cell in range(1,system.ncells-1,step=1)
#            print("Cell ", cell,"\n")

            # Center position of current cell
            x_pos = dx * (Float64(cell) - 0.5) + system.x_min

            # 1.- SELECT NUMBER OF PARTICLES
            #   1.1.- List with particle count
            n_part_list = Int64[]
            for s in coll_group.colliding_species
                if !s.is_background_species
                    push!(n_part_list, length(s.particle_grid_list[cell]))
#                    print("  - ",s.name," ", length(s.particle_grid_list[cell]),"\n")
                end
            end
            #   1.2.- Set number of particles
            n_part_min = minimum(n_part_list)
            if n_part_min == 0
                # If no particles, move to next cell
                continue
            end
            
            # 2.- COMPUTE MAXIMUM NUMBER OF COLLISIONS
            #   2.1.- Compute density product
            dens_prod = 1.0
            for s in coll_group.colliding_species 
                if s.is_background_species
                    n = ReplaceExpressionValues(s.dens_spatial_distribution, system, pos=x_pos)
                    n *= s.init_dens
#                    print("  - Spatial dist eval", n,"\n")
                else
                    # Number of particles of species s in cell
                    n_parts = length(s.particle_grid_list[cell])
                    # Density
                    n = n_parts * s.weight * idx
                end
#                print("  - dens ", s.name," ",n,"\n")
                # Product
                dens_prod *= n
            end
#            print("  - gsigma-max ", gsigma_max,"\n")
#            print("  - dt ", dt ,"\n")
            #   2.2.- Max. collision events (float) 
            N_coll_max_real = dt * dens_prod * gsigma_max * iweight_max * dx 
            #   2.3.- Max. collision events (int)
            N_coll_max = ceil(Int64, N_coll_max_real)
#            print("  - N_coll_max_real ", N_coll_max_real, " ", N_coll_max,"\n")
            if N_coll_max > n_part_min
                PrintErrorMessage(system, "More collisions than particles")
                return c_error
            end
            #   2.4.- Correction factor (Boyd) - USED LATER
            f_B = N_coll_max_real / N_coll_max

            # 3.- COMPUTE N_coll_max COLLISIONS
            #   3.1.- Set species and point to particle lists
            species1 = coll_group.colliding_species[1]
            if species1.is_background_species
                grid_list1 = GenerateBackgroundPartList(species1, N_coll_max, x_pos)
            else
                shuffle!(species1.particle_grid_list[cell])
                grid_list1 = species1.particle_grid_list[cell]
            end
            species2 = coll_group.colliding_species[2]
            if species2.is_background_species
                grid_list2 = GenerateBackgroundPartList(species2, N_coll_max, x_pos)
            else
                shuffle!(species2.particle_grid_list[cell])
                grid_list2 = species2.particle_grid_list[cell]
            end

            #   3.3.- Loop N_coll_max collisions
            for (item, (part1, part2)) in enumerate(zip(grid_list1, grid_list2)) 

                # 4.- COMPUTE COLLISION TYPE
                #   4.1.- Velocity difference
                g_vec = part1.vel .- part2.vel
                g = norm(g_vec)
#                print("  - g ", g,"\n")

                #   4.2.- Interpolate cross-section data
                cross_section_data = Float64[]
                for colltype in collision_list
                    # Cross-section vs. g data tables
                    s_data = colltype.cross_section_data
                    g_data = colltype.energy_data

                    # g value is at extrema
#                    print("  - ", colltype.name, " g_min ", g_data[1]," ", g < g_data[1],"\n")
                    if g < g_data[1] 
                        # Case where g is beyond lowest energy data
                        push!(cross_section_data, s_data[1])
                        continue
                    end
#                    print("  - ", colltype.name, " g_max ", g_data[end]," ", g > g_data[end],"\n")
                    if g > g_data[end] 
                        # Case where g is beyond highest energy data
                        push!(cross_section_data, s_data[end])
                        continue
                    end

                    # g withing data range. Found surrounding g data
                    ix_min = findlast(x -> x <= g, g_data)
                    ix_max = ix_min + 1#findlast( x -> x > g, g_data)
#                    print("    - index ", ix_min,"  ", ix_max,"\n")
                    g_min = g_data[ix_min]
                    g_max = g_data[ix_max]
                    s_min = s_data[ix_min]
                    s_max = s_data[ix_max]
                    interp_list = [g_min , g_max, s_min, s_max]
#                    print("    - Interp limints ",interp_list,"\n")
                    cross_section = linear_interpolation(interp_list, g)
#                    print("    - SIGMA: ",cross_section,"\n")
                    push!(cross_section_data, cross_section)

                end
#                print("  - Cross sections ", cross_section_data,"\n")
                #cross_section_total = sum(cross_section_data)

                #   4.4.- Compute maximum collision energy
                P_data = g * cross_section_data / gsigma_max * f_B
#                print("  - P_max ", sum(P_data),"\n")

                #   4.5.- Select collision 
                R1 = rand()
#                print("  - Rand 1 ", R1,"\n")
                P_min = 0.0
                P_max = 0.0
                for (i, colltype) in enumerate(collision_list)
                    P_max += P_data[i]
#                    print("  - ", colltype.name," ", P_min," - ", P_max,"\n")
                    if R1 >= P_min && R1 < P_max
                        # Run this collision
                        colltype.collfunction(coll_group, colltype, part1,
                            part2, g, item, cell)
                        # Add diagnostic
                        colltype.diagnostic[cell] += 1
                        break
                    end
                    P_min = P_max
                end  # Collision types
#                print("\n\n")
            end  # N_max collisions
        end  # Cells
    end # Collision group

    return 0
end


function GenerateBackgroundPartList(species::Species, N_coll_max::Int64, cell_pos::Float64)
    part_list = Particle[]
    Tn = species.init_temp
    m = species.mass
    dev = sqrt(kb * Tn / m)
    for i in range(1,N_coll_max,step=1)
        part = Particle()
        part.pos = cell_pos
        part.vel = randn(3) * dev 
        push!(part_list, part)
    end
    return part_list
end

end