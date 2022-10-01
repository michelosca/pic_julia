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
    charge_dens::Vector{Float64}
    
    particle_list::Vector{Particle}
    particle_grid_list::Vector{Particle}
    particle_count::Int64

    is_background_gas::Bool
    temp::Float64

    Species() = new()
end

mutable struct System

    dt::Float64
    t_start::Float64
    time::Float64
    t_end::Float64
    
    step::Int64
    
    ncells::Int64
    x_min::Float64
    x_max::Float64
    Lx::Float64
    dx::Float64

    bc_field::Int64
    bc_part::Int64

    V0_min::Float64
    V0_max::Float64
    System() = new()
end

mutable struct PlasmaParameters
    total_dens::Vector{Float64}
    total_charge_dens::Vector{Float64}
    electric_potential::Vector{Float64}
    electric_field::Vector{Float64}
    PlasmaParameters() = new()
end

mutable struct Field
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    Field() = new()
end

end