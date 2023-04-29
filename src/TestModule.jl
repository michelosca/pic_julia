module TestModule

using Plots
using HDF5
using Printf

using SharedData: Species, System
using Constants: c_field_electric, c_field_pot, c_field_rho
using Constants: c_bc_open
using Constants: epsilon_0, e, me, kb
using Tools: linear_interpolation

function test_species_densities(species_list::Vector{Species}, system::System)

    x = zeros(Float64, system.ncells)
    dx = system.dx
    x_min = system.x_min
    c_min = system.cell_min
    c_max = system.cell_max
    for i in range(1,system.ncells,step=1)
        x[i] = dx * (i - 1) - x_min 
    end
    p = plot(title="Species densities")
    for s in species_list
        if !s.is_background_species
            plot!(p, x, s.dens[c_min:c_max], label=s.name)
        end
    end
    plot!(p
        , ylabel="Particle number density [m-3]"
        , xlabel="Position [m]"
    )

    png("species_number_densities")
    return p
end

function test_plot_field(data::Vector{Float64}, system::System, id::Int64)

    x_min = system.x_min
    x_max = system.x_max
    c_min = system.cell_min
    c_max = system.cell_max
    dx = system.dx
    x = range(x_min, x_max, step=dx)
    if id == c_field_rho
        plot(x, data[c_min:c_max]
            , xlabel = "Position [m]"
            , ylabel = "C/m^3"
            , title = "Charge density"
        )
        png("charge_density_total")
    elseif id == c_field_pot
        plot(x, data[c_min:c_max]
            , label = "Simulation"
            , xlabel = "Position [m]"
            , ylabel = "V"
            , title = "Electric potential"
        )

        #pot_theory = potential_theory(x, system)
        #plot!(x, pot_theory
        #    , label = "Theory"
        #    , legend = :bottom
        #)
        png("electric_potential")
    elseif id == c_field_electric
        plot(x, data[c_min:c_max]
            , label = "Simulation"
            , xlabel = "Position [m]"
            , ylabel = "V/m"
            , title = "Electric field"
        )

        #ex_th = ex_theory(x, system)
        #plot!(x, ex_th
        #    , label = "Theory"
        #    , legend = :bottom
        #)
        png("electric_field")
    end
end

function potential_theory(x, system::System)

    x_min = system.x_min
    x_max = system.x_max
    V0_min = system.V0_min
    V0_max = system.V0_max

    rho = 5
    # Derivation
    # d^2 phi / dx^2 = -rho / eps0
    # d phi / dx = -rho/eps0*x + C1
    # phi = -x^2rho/eps0/2 + C1*x + C2
    # x  = x_min -> phi = V0_min; V0_min = -x_min^2*rho/eps0/2 + C1*x_min + C2
    # x  = x_max -> phi = V0_max; V0_max = -x_max^2*rho/eps0/2 + C1*x_max + C2
    #                     V0_max - V0_min = -(x_max^2 - x_min^2)*rho/eps0/2 + C1(x_max - x_min)
    # C1 = ( (V0_max - V0_min)+(x_max^2 - x_min^2)*rho/eps0/2 ) / (x_max - x_min)
    # C2 = V0_min + x_min^2*rho/eps0/2 - C1*x_min

    C1 = ( (V0_max - V0_min)+(x_max^2 - x_min^2)*rho/epsilon_0/2 ) / (x_max - x_min)
    C2 = V0_min + x_min^2*rho/epsilon_0/2 - C1*x_min
    pot_theory = -x.^2*rho/epsilon_0/2 .+ C1*x .+ C2

    return pot_theory
end

function ex_theory(x, system::System)

    x_min = system.x_min
    x_max = system.x_max
    V0_min = system.V0_min
    V0_max = system.V0_max

    rho = 5
    # Derivation
    # d^2 phi / dx^2 = -rho / eps0
    # Ex = -d phi / dx = rho/eps0*x - C1
    # C1 = ( (V0_max - V0_min)+(x_max^2 - x_min^2)*rho/eps0/2 ) / (x_max - x_min)

    C1 = ( (V0_max - V0_min)+(x_max^2 - x_min^2)*rho/epsilon_0/2 ) / (x_max - x_min)
    ex_theory = rho/epsilon_0 * x .- C1

    return ex_theory
end

function plasma_electron_frequency_cold(dens::Float64)
    return sqrt(dens * e * e / epsilon_0 / me)
end

function plasma_electron_frequency_warm(dens::Float64, temp::Float64, k::Float64)
    omega_e = plasma_electron_frequency_cold(dens)
    thermal_speed = sqrt(kb * temp / me)
    return omega_e*omega_e + 3.0 * (k * thermal_speed)^2
end


function debye_length(dens::Float64, temp_eV::Union{Int64, Float64})
    charge2 = e * e
    eV_to_K = e/kb
    temp_K = temp_eV * eV_to_K
    return sqrt(epsilon_0 * kb * temp_K / dens / charge2)
end


function ColdLangmuirOscillationTest()

    for i in 0:50
        filename = @sprintf("stdout%05i.h5",i)
        h5open("sim/"* filename,"r") do fid
            system = read(fid, "System")
            x = system["Grid"]

            # Plot plasma potential
            ### Simulation rasults
            pot_struct = read(fid, "Electric_Potential")
            pot = pot_struct["Vx"]
            p = plot(x,pot
                , frame = :box
                , linewidth = 2
                , xlabel = "Position / m"
                , ylabel = "Plasma potential / V"
                , ylims = (-1,1)
                , xlims = (x[1], x[end])
                , legend = false
            )
            ### Theory results
            wavelength = 0.1
            k = 2*pi/wavelength
            n0 = 1.e12
            dens_amp = 0.1 * n0
            rho_amp = dens_amp * e
            pot_amp = - rho_amp / epsilon_0 / k / k
            omega_e = plasma_electron_frequency_cold(n0)
            system = fid["System"]
            time = attrs(system)["time"]
            pot_theory = map(x -> pot_amp * sin.(k * x) .* cos(omega_e * time), x)
            plot!(p, x, pot_theory
                , linewidth = 2
                , linestyle = :dot
            )
            fig_name = @sprintf("potential_%05i.png",i)
            savefig(p, "sim/pot/"*fig_name)

            # Plot densities
            dens_struct = read(fid, "Number_Density")
            p = plot(
                frame = :box
                , xlabel = "Position / m"
                , ylabel = "Number density / m^{-3}"
                , ylims = (0.75e12,1.25e12)
                , xlims = (x[1], x[end])
                , legend = true 
            )
            for (species, dens) in dens_struct
                plot!(p, x, dens
                    , linewidth = 2
                    , label = species
                )
            end
            fig_name = @sprintf("dens_%05i.png",i)
            savefig(p, "sim/dens/"*fig_name)
        end
    end
end


function RF_potential_wave()

    for i in 0:200
        filename = @sprintf("stdout%05i.h5",i)
        h5open("sim/"* filename,"r") do fid
            system = read(fid, "System")
            x = system["Grid"]

            # Plot plasma potential
            ### Simulation rasults
            pot_struct = read(fid, "Electric_Potential")
            pot = pot_struct["Vx"]
            p = plot(x,pot
                , frame = :box
                , linewidth = 2
                , xlabel = "Position / m"
                , ylabel = "Plasma potential / V"
                , ylims = (-210,210)
                , xlims = (x[1], x[end])
                , legend = false
            )

            # Analytic solution
            system = fid["System"]
            time = attrs(system)["time"]
            an_pot_left  = 100.0 * sin(4*pi*13.56e6*time)
            an_pot_right = 100.0 * (sin(2.0*pi*13.56e6*time))
            plot!(p, [x[1], x[50]], [an_pot_left, an_pot_left], c = :black, lw = 2)
            plot!(p, [x[50], x[end]], [an_pot_right, an_pot_right], c = :green, lw = 2)
            fig_name = @sprintf("potential_%05i.png",i)
            savefig(p, "sim/pot/"*fig_name)
        end
    end
end

function MCC_test()

    p = plot(
        frame = :box
        , xlabel = "Position / m"
        , ylabel = "Number density / m^{-3}"
        , legend = true 
    )
    ls_list = [:solid, :dot]
    for i in 1:1
        filename = @sprintf("stdout%05i.h5",i)
        h5open("sim/"* filename,"r") do fid

            # NC grid distribution
            NC_struct = read(fid, "Neutral Collisions")
            x = NC_struct["NC_Grid"]

            # Plot collision rate 
            for (speciescomb, set) in NC_struct
                print(speciescomb,"\n")
                if speciescomb == "NC_Grid"
                    continue
                end
                for colltype in keys(set)
                    print("  - ", colltype,"\n")
                    plot!(x, set[colltype]
                        , label = speciescomb * "-" *colltype
                    )
                end

            end
            #for (species, dens) in dens_struct
            #    if species == "electrons"
            #        plot!(p, x, dens
            #            , linewidth = 2
            #            , label = species
            #            , xlims = (x[1], x[end])
            #            , ls = ls_list[i+1]
            #        )
            #    end
            #end
        end
    end
    #fig_name = @sprintf("dens_mcc.png")
    #savefig(p, "sim/"*fig_name)
    return p
end

function RateCoefficient(energy_data::Vector{Float64},
    cross_section_data::Vector{Float64})
    # From Lieberman and Lichtenberg (2005). Chapter 3. Section 3.5. Equation 3.5.2
    # Cross section convolution with Maxwellian distribution
    me = 9.10938188e-31
    kb = 1.380649e-23
    qe = 1.602176634e-19
    eV_to_K = qe/kb

    Te_eV_data= 10.0 .^ (-1:0.1:3)
    K_data = Float64[]
    for Te_eV in Te_eV_data

        # Generate velocity range for Maxwelllian DF
        # Set velocity limits based on current Te_eV
        v_max = sqrt(2 * (Te_eV*qe) / me) * 5
        dv = v_max / 100
        Te = Te_eV * eV_to_K
        v_list = 0:dv:v_max

        # Convolution sigma and Maxwellian DF
        K = 0.0
        for v_point in v_list[2:end]
            v = v_point - 0.5*dv
            # Interpolate cross-section data
            e = 0.5 * me * v * v
            if e >= energy_data[end]
                sigma = cross_section_data[end]
            elseif e <= energy_data[1]
                sigma = cross_section_data[1]
            else
                ix_min = findlast(x -> x <= e, energy_data)
                ix_max = ix_min + 1#findlast( x -> x > g, g_data)
                e_min = energy_data[ix_min]
                e_max = energy_data[ix_max]
                s_min = cross_section_data[ix_min]
                s_max = cross_section_data[ix_max]
                interp_list = [e_min , e_max, s_min, s_max]
                sigma = linear_interpolation(interp_list, e)
            end

            # Integration bit
            K += sigma * v^3 * exp(- me * v * v / 2 / kb / Te) 
        end
        K *= (me / 2 / pi / kb / Te)^(3/2) * 4 * pi * dv
        push!(K_data, K)
    end


    return Te_eV_data, K_data
end

function ReadCrossSectionData(filepath::String, energy_units::Float64,
    sigma_units::Float64)

    energy_data = Float64[]
    sigma_data = Float64[]
    open(filepath) do file
        while ! eof(file)
            str = strip(readline(file))
            ix = findfirst('\t', str)
            energy_str = strip(str[1:ix-1])
            sigma_str = strip(str[ix+1:end])

            energy = parse(Float64, energy_str) * energy_units
            push!(energy_data, energy)
    
            # Cross section
            sigma = parse(Float64, sigma_str) * sigma_units
            push!(sigma_data, sigma)
        end
    end

    return energy_data, sigma_data
end

function MCC_RateCoefficients()
    main = "/home/moe505/Documents/pic_julia/tests/rate_coefficient_Ar_test/"
    fe_elast = "e_Ar_elastic.table"
    fe_excit = "e_Ar_excitation.table"
    fe_ioniz = "e_Ar_ionization.table"

    data_list = [fe_elast, fe_excit, fe_ioniz]

    p = plot()

    # Analytic results
    for data_path in data_list
        e_data, s_data = ReadCrossSectionData(main * data_path, e, 1.0)
        Te_list, K_list = RateCoefficient(e_data, s_data)
        plot!(p, Te_list, K_list
            , lw = 2
            , label = data_path
        )
    end

    # Simulation results
    Lx = 0.025
    W_e = 7.e8
    ncells = 400
    npart_per_cell = 100
    Lx = 0.025
    e_dens = W_e * npart_per_cell * ncells / Lx
    Ar_dens = 2.069e21
    # 1.- Gather data
    K_elastic = 0.0
    K_excitation = 0.0
    K_ionization = 0.0
    n_files = 101
    for i in 0:(n_files-1)
        filename = @sprintf("stdout%05i.h5",i)
        h5open("sim/"* filename,"r") do fid
            nc = read(fid, "Neutral_Collisions")
            system = fid["System"]
            dt = attrs(system)["dt"]
            dx = attrs(system)["dx"]

            # Electron super-part. weight
            #part = fid["Particles"]
            #part_e = part["e"]
            #W_e = attrs(part_e)["weight"]

            nc_e = nc["e:Ar"]
            # Elastic scattering
            nc_e_elastic    = sum(nc_e["Elastic"])
            K_elastic += nc_e_elastic/Lx / dt * W_e / Ar_dens / e_dens

            nc_e_excitation = sum(nc_e["Excitation"])
            K_excitation += nc_e_excitation/Lx / dt * W_e / Ar_dens / e_dens

            nc_e_ionization = sum(nc_e["Ionization"])
            K_ionization += nc_e_ionization/Lx / dt * W_e / Ar_dens / e_dens
        end
    end
    K_elastic /= n_files
    K_excitation /= n_files
    K_ionization /= n_files
    print("K elastic    ", K_elastic,"\n")
    print("K excitation ", K_excitation,"\n")
    print("K ionization ", K_ionization,"\n")

    Te_eV = 15
    scatter!(p, [Te_eV], [K_elastic]
        , markersize = 4
        , marker = :circle
        , markercolor = 1
        , markerstrokecolor = 1
    )

    scatter!(p, [Te_eV], [K_excitation]
        , markersize = 4
        , marker = :circle
        , markercolor = 2
        , markerstrokecolor = 2
    )

    scatter!(p, [Te_eV], [K_ionization]
        , markersize = 4
        , marker = :circle
        , markercolor = 3
        , markerstrokecolor = 3
    )

    plot!(p
        , size = (500, 300)
        , frame = :box
        , xlims = (1.e-1, 1.e3)
        , ylims = (1.e-22, 1.e-12)
        , yscale=:log10
        , xscale=:log10
        , legend = :bottomright
    )
    return p
end

end