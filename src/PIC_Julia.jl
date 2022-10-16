module PIC_Julia 

using SharedData: System, Species, OutputBlock
using Inputs: SetupInputData!
using GridData: GetTotalChargeDensity
using ElectricMagneticFields: UpdateElectricField!, GetElectricPotential
using ElectricMagneticFields: InitializeElectricField, InitializeMagneticField
using ParticleIntegrator: IntegrateParticlesPhaseSpace!
#using Tools: RealocateParticlesToGridList!
#using Tools: RealocateParticlesToMainList!
using Outputs: GenerateOutputs! 

using TestModule: test_species_densities, test_plot_field
using Constants: c_field_electric, c_field_pot, c_field_rho
using Constants: c_error

function run_pic(input_file::String)

    # Create system structure
    system = System()

    # Initialize output block list
    output_list = OutputBlock[]

    # Initialize species data list
    species_list = Species[]

    errcode = SetupInputData!(input_file::String, species_list::Vector{Species}, 
        system::System, output_list::Vector{OutputBlock})
    if errcode == c_error
        return c_error
    end

    # Load electric and magnetic field data
    electric_field = InitializeElectricField(system)
    magnetic_field = InitializeMagneticField(system)

    # 01 - push back particles velocity
    #

    # 02 - Initialize fields: densities, charge densities, potential, E/B-fields
    # 1.- Get charge density
    charge_density = GetTotalChargeDensity(species_list, system) 
    # 2.- Calculate electric potential
    electric_potential = GetElectricPotential(charge_density, system)
    # 3.- Calculate electric field
    charge_dens_min = charge_density[system.cell_min]
    charge_dens_max = charge_density[system.cell_max]
    UpdateElectricField!(electric_field, electric_potential, charge_dens_min,
        charge_dens_max, system)
    #UpdateMagneticField!(magnetic_field, system)

    # 03.- Output initial conditions
    GenerateOutputs!(output_list, species_list, system, electric_potential, electric_field)

    while system.step <= system.step_end
        # 1.- Get charge density
        charge_density = GetTotalChargeDensity(species_list, system) 

        # 2.- Calculate electric potential
        electric_potential = GetElectricPotential(charge_density, system)

        # 3.- Calculate electric field
        charge_dens_min = charge_density[system.cell_min]
        charge_dens_max = charge_density[system.cell_max]
        UpdateElectricField!(electric_field, electric_potential, charge_dens_min,
            charge_dens_max, system)
        #UpdateMagneticField!(magnetic_field)

        # 4.- Update particle pos and vel
        IntegrateParticlesPhaseSpace!(species_list, system, electric_field,
            magnetic_field)

        # 5.- Boundary conditions
        # -> particle BC applied while integrating Phase-Space 

        # 6.- Collisions: 6.1.- gather particles in cells
        #RealocateParticlesToGridList!(species_list)
        #for species in species_list
        #    print(species.name, " secondary list ", length(species.particle_grid_list),"\n")
        #end
        #                 6.2.- collider
        #                 6.3.- 
        #RealocateParticlesToMainList!(species_list)

        # 6.- Update time
        system.time += system.dt
        system.step += 1

        # 7.- Output data
        GenerateOutputs!(output_list, species_list, system, electric_potential, electric_field)

    end

    return
end

end