module SharedData 


mutable struct Particle
    pos::Float64
    vel::Vector{Float64}

    Particle() = new()
end

mutable struct Species
    id::Int64
    name::String
    mass::Float64
    charge::Float64
    weight::Float64

    dens::Vector{Float64}
    
    particle_list::Vector{Particle}
    particle_grid_list::Vector{Vector{Particle}}
    particle_count::Int64
    part_per_cell::Int64

    is_background_species::Bool
    init_temp::Float64
    init_dens::Float64
    dens_spatial_distribution::Union{Int64, Float64, Expr}

    Species() = new()
end

mutable struct System

    dt::Float64
    t_start::Float64
    time::Float64
    t_end::Float64
    
    step::Int64
    step_start::Int64
    step_end::Int64
    
    ncells::Int64
    cell_min::Int64
    cell_max::Int64
    gc::Int64
    ncells_total::Int64

    x_min::Float64
    x_max::Float64
    Lx::Float64
    dx::Float64

    bc_field_min::Int64
    bc_field_max::Int64
    bc_part_min::Int64
    bc_part_max::Int64

    mcc::Bool

    folder::String
    log_file::String

    System() = new()
end

mutable struct Field

    id::Int64
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    bc_min::Int64
    bc_max::Int64


    Field() = new()
end

mutable struct OutputDataStruct

    id::Int64
    name::String

    species_id::Int64
    species_name::String

    dir_id::Int64

    data::Array{Float64,}

    OutputDataStruct() = new()
end

mutable struct OutputBlock
    
    name::String

    t_start::Float64
    t_end::Float64
    dt::Float64
    step_start::Int64
    step_end::Int64
    step_jump::Int64

    averaged::Bool
    dt_av::Float64
    step_av::Int64

    param_list::Vector{OutputDataStruct}

    file_id::Int64
    zero_pad::Int64
    time_dump::Float64
    step_dump::Int64
    time_av_start::Int64
    time_av_end::Int64
    step_av_start::Int64
    step_av_end::Int64


    OutputBlock() = new()

end

mutable struct Waveform

    amp::Float64
    freq::Float64
    wavefunction::Union{Int64, Float64, Expr}
    boundary::Int64
    t_start::Float64
    t_end::Float64

    Waveform() = new()
end

mutable struct Collision
    name::String
    reactants::Vector{Species}
    products::Vector{Species}
    species::Vector{Species}
    species_balance::Vector{Int64}
    energy_threshold::Float64

    energy_data::Vector{Float64}
    energy_units::Float64
    cross_section_data::Vector{Float64}
    cross_section_units::Float64

    diagnostic::Vector{Int64}

    Collision() = new()
end

mutable struct CollisionGroup

    collision_list::Vector{Collision}
    colliding_species::Vector{Species}

    gsigma_max::Float64
    part_weight_max::Float64
    reduced_mass::Float64

    CollisionGroup() = new()
end

end