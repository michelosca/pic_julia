module PIC_Julia 

using Inputs: GetSystemParameters, GetSpeciesList
using PrintModule: PrintSpecies
using GridData: GetTotalChargeDensity
using ElectricMagneticFields: UpdateElectricField!, GetElectricPotential
using ElectricMagneticFields: InitializeElectricField, InitializeMagneticField
using ParticleIntegrator: IntegrateParticlesPhaseSpace!
#using Tools: RealocateParticlesToGridList!
#using Tools: RealocateParticlesToMainList!

function run_pic()
    system = GetSystemParameters()
    species_list = GetSpeciesList(system)
    #for species in species_list
    #    PrintSpecies(species)
    #end
    electric_field = InitializeElectricField(system)
    magnetic_field = InitializeMagneticField(system)

    # [ ] test particles temperatures
    
    # 01 - push back particles velocity
    #while system.time <= system.t_end
        # 1.- Get charge density
        charge_density = GetTotalChargeDensity(species_list, system) 

        # 2.- Calculate electric potential
        electric_potential = GetElectricPotential(charge_density, system)

        ## 3.- Calculate electric field
        charge_dens_min = charge_density[1]
        charge_dens_max = charge_density[end]
        UpdateElectricField!(electric_field, electric_potential, charge_dens_min,
            charge_dens_max, system)
        #UpdateMagneticField!(magnetic_field)

        ## 4.- Update particle pos and vel
        IntegrateParticlesPhaseSpace!(species_list, system, electric_field,
            magnetic_field)

        # 5.- Boundary conditions
        # -> particle BC applied in the integrator

        # 6.- Collisions: 6.1.- gather particles in cells
        #RealocateParticlesToGridList!(species_list)
        #for species in species_list
        #    print(species.name, " secondary list ", length(species.particle_grid_list),"\n")
        #end
        #                 6.2.- collider
        #                 6.3.- 
        #RealocateParticlesToMainList!(species_list)
        # 6.- Update time
        #system.time += system.dt
    #end
    return charge_density, electric_potential, electric_field#, electric_potential, electric_field
end

end