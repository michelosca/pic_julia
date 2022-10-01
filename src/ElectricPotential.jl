module ElectricPotential

using SharedData: System, Species, Field
using Constants: epsilon_0
using Tools: GetNumberDensity
using SparseArrays: sparse
using LinearSolve: solve, LinearProblem

function GetElectricPotential(charge_density::Vector{Float64}, system::System)

    ncells = system.ncells-2
    dx = system.dx

    mid_diag = (range(1,ncells,step=1), range(1,ncells,step=1))
    upper_diag =  (range(1,ncells-1,step=1), range(2,ncells,step=1))
    bottom_diag = (range(2,ncells,step=1), range(1,ncells-1,step=1))
    matrix = sparse(cat(dims=1,mid_diag[1],upper_diag[1],bottom_diag[1]),
        cat(dims=1, mid_diag[2], upper_diag[2], bottom_diag[2]),
        cat(dims=1, ones(ncells)*-2.0, ones(2*(ncells-1))))

    charge_density_input = charge_density[2:end-1].*-dx*dx/epsilon_0
    charge_density_input[1] -= system.V0_min
    charge_density_input[end] -= system.V0_max
    prob = LinearProblem(matrix, charge_density_input)
    sol = solve(prob)
    electric_potential_output = sol.u
    electric_potential = cat(dims=1,[system.V0_min], electric_potential_output,
        [system.V0_max])
    return electric_potential 
end


function GetElectricField(electric_potential::Vector{Float64},
    charge_dens_min::Float64, charge_dens_max::Float64, system::System)

    dx = system.dx
    idx = 1.0/dx
    ncells = system.ncells
    electric_field = zeros(ncells)
    for i in range(2,ncells-1, step=1)
        electric_field[i] = electric_potential[i-1] - electric_potential[i+1]
    end
    electric_field .*= 0.5*idx
    electric_field[1] = electric_field[2] - charge_dens_min * dx /epsilon_0
    electric_field[end] = electric_field[end-1] + charge_dens_max * dx /epsilon_0
    return electric_field
end

function UpdateElectricField!(e_field::Field,
    electric_potential::Vector{Float64}, rho_min::Float64, rho_max::Float64,
    system::System)
    
    e_field.x = GetElectricField(electric_potential, rho_min, rho_max, system)
    # Changes in time in e_field.y or z should go here, e.g. inductive heating

end

function InitializeElectricField(system::System)

    electric_field = Field()
    electric_field.x = zeros(Float64, system.ncells)
    electric_field.y = zeros(Float64, system.ncells)
    electric_field.z = zeros(Float64, system.ncells)
    return electric_field
end

function InitializeMagneticField(system::System)

    magnetic_field = Field()
    magnetic_field.x = ones(Float64, system.ncells)
    magnetic_field.y = zeros(Float64, system.ncells)
    magnetic_field.z = zeros(Float64, system.ncells)
    return magnetic_field 
end

function GetTotalChargeDensity(species_list::Vector{Species}, system::System)

    charge_density = zeros(Float64, system.ncells)
    for species in species_list
        if species.is_background_gas
            continue
        end
        dens = GetNumberDensity(species, system)
        charge_density .+= dens .* species.charge
    end
    return charge_density
end

end