module ElectricMagneticFields

using SharedData: System, Species, Field, Waveform
using Constants: epsilon_0
using Constants: c_field_electric, c_field_magnetic
using Constants: c_bc_periodic, c_bc_open
using Constants: c_bc_x_min, c_bc_x_max
using SparseArrays: sparse
using LinearSolve: solve, LinearProblem
using EvaluateExpressions: ReplaceExpressionValues


function GetElectricPotential(charge_density::Vector{Float64},
    waveform_list::Vector{Waveform}, system::System)

    V0_min, V0_max = GetPotential_boundaries(waveform_list, system)
    dx = system.dx
    ncells = system.ncells - 2
    cell_start = system.cell_min + 1
    cell_end = system.cell_max - 1

    # Prepare linear solver matrix
    mid_diag = (range(1,ncells,step=1), range(1,ncells,step=1))
    upper_diag =  (range(1,ncells-1,step=1), range(2,ncells,step=1))
    bottom_diag = (range(2,ncells,step=1), range(1,ncells-1,step=1))
    matrix = sparse(cat(dims=1,mid_diag[1],upper_diag[1],bottom_diag[1]),
        cat(dims=1, mid_diag[2], upper_diag[2], bottom_diag[2]),
        cat(dims=1, ones(ncells)*-2.0, ones(2*(ncells-1))))

    # Prepare charge density array
    charge_density_input = charge_density[cell_start:cell_end].*-dx*dx/epsilon_0
    charge_density_input[1] -= V0_min
    charge_density_input[end] -= V0_max

    # Solve for electric potential 
    prob = LinearProblem(matrix, charge_density_input)
    sol = solve(prob)

    # Dump solution into electric potential array
    electric_potential_output = sol.u

    # Extent electric potential array to include  boundary & ghost cells
    gc = system.gc
    gc_left = zeros(Float64, gc+1)
    gc_left[end] = V0_min
    gc_right = zeros(Float64, gc+1)
    gc_right[1] = V0_max
    electric_potential = cat(dims=1, gc_left, electric_potential_output, gc_right)

    # Apply boundary conditions to the electric potential array
    ApplyElectricPotentialBC!(electric_potential, system)

    return electric_potential 
end

function ApplyElectricPotentialBC!(electric_potential::Vector{Float64}, system::System)

    cell_min = system.cell_min
    cell_max = system.cell_max
    gc = system.gc

    if system.bc_field_min == c_bc_periodic
        for i in range(1,gc,step=1)
            electric_potential[cell_min - i] = electric_potential[cell_max - 1]
        end
    end

    if system.bc_field_max == c_bc_periodic
        for i in range(1,gc,step=1)
            electric_potential[cell_max + i] = electric_potential[cell_min + i]
        end
    end

end

function GetElectricField(electric_potential::Vector{Float64},
    system::System)

    dx = system.dx
    idx = 1.0/dx
    cell_min = system.cell_min
    cell_max = system.cell_max
    electric_field = zeros(Float64, system.ncells_total)
    for i in range(cell_min, cell_max, step=1)
        electric_field[i] = electric_potential[i-1] - electric_potential[i+1]
    end
    electric_field .*= 0.5*idx

    return electric_field
end

function ApplyElectricFieldBC!(electric_field::Field, system::System,
    charge_dens_min::Float64, charge_dens_max::Float64)

    ex = electric_field.x
    cell_min = system.cell_min
    cell_max = system.cell_max
    gc = system.gc
    if system.bc_field_min == c_bc_open
        dx = system.dx
        ex[cell_min] = ex[cell_min+1] - charge_dens_min * dx /epsilon_0

    elseif system.bc_field_min == c_bc_periodic
        for i in range(1,gc,step=1)
            ex[cell_min-i] = ex[cell_min+i]
        end
    end

    if system.bc_field_max == c_bc_open
        dx = system.dx
        ex[cell_max] = ex[cell_max-1] + charge_dens_max * dx /epsilon_0

    elseif system.bc_field_max == c_bc_periodic
        for i in range(1,gc,step=1)
            ex[cell_max+i] = ex[cell_max-i]
        end
    end
end

function UpdateElectricField!(e_field::Field,
    electric_potential::Vector{Float64}, rho_min::Float64, rho_max::Float64,
    system::System)
    
    e_field.x = GetElectricField(electric_potential, system)

    # Changes in time in e_field.y or z should go here, e.g. inductive heating

    ApplyElectricFieldBC!(e_field, system, rho_min, rho_max)

end

function InitializeElectricField(system::System)

    electric_field = Field()
    electric_field.id = c_field_electric
    electric_field.x = zeros(Float64, system.ncells_total)
    electric_field.y = zeros(Float64, system.ncells_total)
    electric_field.z = zeros(Float64, system.ncells_total)
    return electric_field
end

function InitializeMagneticField(system::System)

    magnetic_field = Field()
    magnetic_field.id = c_field_magnetic
    magnetic_field.x =  ones(Float64, system.ncells_total)
    magnetic_field.y = zeros(Float64, system.ncells_total)
    magnetic_field.z = zeros(Float64, system.ncells_total)
    return magnetic_field 
end

function GetPotential_boundaries(waveform_list::Vector{Waveform}, system::System)

    V0_min = 0.0
    V0_max = 0.0

    for (i, waveform) in enumerate(waveform_list)
        if waveform.boundary == c_bc_x_min
            V0_min += ReplaceExpressionValues(waveform.wavefunction, system,
                waveform=waveform) * waveform.amp
        elseif waveform.boundary == c_bc_x_max
            V0_max += ReplaceExpressionValues(waveform.wavefunction, system,
                waveform=waveform) * waveform.amp
        end
    end
    return V0_min, V0_max
end

end