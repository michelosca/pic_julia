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

module PIC_Julia 

using SharedData: System, Species, OutputBlock, Waveform, CollisionGroup
using Inputs: SetupInputData!
using GridData: GetTotalChargeDensity
using ElectricMagneticFields: UpdateElectricField!, GetElectricPotential
using ElectricMagneticFields: InitializeElectricField, InitializeMagneticField
using ParticleIntegrator: IntegrateParticlesPhaseSpace!
using Tools: RealocateParticlesToGridList!
using Tools: RealocateParticlesToMainList!
using Outputs: GenerateOutputs! 
using Constants: c_error
using NeutralCollisions: NeutralCollisions!

using PrintModule: PrintCollisionGroup
using PrintModule: PrintErrorMessage

function run_pic(input_file::String)

    # Create system structure
    system = System()

    # Initialize output block list
    output_list = OutputBlock[]

    # Initialize species data list
    species_list = Species[]

    # Initialize waveform list
    waveform_list = Waveform[]

    # Initialize collision list
    collision_list = CollisionGroup[]

    errcode = SetupInputData!(input_file
        , species_list
        , system
        , output_list
        , waveform_list
        , collision_list
    )
    if errcode == c_error
        PrintErrorMessage(system,"Setting up input data")
        return c_error
    end


    if system.mcc
        for coll_group in collision_list
            PrintCollisionGroup(coll_group)
        end
    end

    # Load electric and magnetic field data
    electric_field = InitializeElectricField(system)
    magnetic_field = InitializeMagneticField(system)

    # 01 - push back particles velocity
    #

    # 02 - Initialize fields: densities, charge densities, potential, E/B-fields
    ### 1.- Get charge density
    charge_density = GetTotalChargeDensity(species_list, system) 
    ### 2.- Calculate electric potential
    electric_potential = GetElectricPotential(charge_density, waveform_list, system)
    ### 3.- Calculate electric field
    charge_dens_min = charge_density[system.cell_min]
    charge_dens_max = charge_density[system.cell_max]
    UpdateElectricField!(electric_field, electric_potential, charge_dens_min,
        charge_dens_max, system)
    #UpdateMagneticField!(magnetic_field, system)

    # 03.- Output initial conditions
    GenerateOutputs!(output_list, species_list, system, electric_potential, electric_field)

    while system.step < system.step_end

        # 1.- Get charge density
        charge_density = GetTotalChargeDensity(species_list, system) 

        # 2.- Calculate electric potential
        electric_potential = GetElectricPotential(charge_density, waveform_list, system)

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
        if system.mcc
            RealocateParticlesToGridList!(species_list, system)
            errcode = NeutralCollisions!(collision_list, species_list, system)
            RealocateParticlesToMainList!(species_list)
        end

        # 7.- Update time
        system.time += system.dt
        system.step += 1

        # 8.- Output data
        GenerateOutputs!(output_list, species_list, system, electric_potential, electric_field)

    end

    return
end

end