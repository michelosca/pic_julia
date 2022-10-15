module TestModule

using Plots
using SharedData: Species, System
using Constants: c_field_electric, c_field_pot, c_field_rho
using Constants: c_bc_open
using Constants: epsilon_0, e, me, kb

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

end