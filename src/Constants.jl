module Constants

const me = 9.1093837015e-31   # kg
const e = 1.602176634e-19    # C
const kb = 1.380649e-23       # J/K
const amu = 1.6605390666e-27  # kg
const epsilon_0 = 8.8541878128e-12 #F m^-1
const K_to_eV = kb / e

const gc_triangle = 2 # ghost cells, this must be changed if particle-to-grid shapes change

const c_error = 1

const c_bc_periodic = 100
const c_bc_open = 101
const c_bc_x_min = 102
const c_bc_x_max = 103

const c_stag_centre = 200 
const c_stag_right = 201

const c_field_electric = 300
const c_field_magnetic = 301
const c_field_rho = 302
const c_field_pot = 303

const c_block_system = 400
const c_block_species = 401
const c_block_output = 402
const c_block_constants = 403
const c_block_waveform = 404

const c_o_all_species = 500
const c_o_none_species = 501
const c_o_phase_space = 502
const c_o_density = 503
const c_o_potential = 504
const c_o_electric_field = 505
const c_o_magnetic_field = 506
const c_o_neutral_collisions = 507
const c_o_probe = 508

const c_dir_none = 600
const c_dir_x = 601
const c_dir_y = 602
const c_dir_z = 603

end