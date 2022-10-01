module PrintModule

using Printf
using SharedData: Species

function PrintSpecies(species::Species)
    @printf("Species name: %s\n", species.name)
    @printf("  - id: %i\n", species.id)
    @printf("  - mass: %g kg\n", species.mass)
    @printf("  - charge: %g C\n", species.charge)
    @printf("  - part. weight: %g\n", species.weight)
    @printf("  - particle count: %i\n", species.particle_count)
    @printf("  - is background: %s\n", species.is_background_gas)
end

end