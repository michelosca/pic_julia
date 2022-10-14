module Constants

const me = 9.1093837015e-31   # kg
const e = 1.602176634e-19    # C
const kb = 1.380649e-23       # J/K
const amu = 1.6605390666e-27  # kg
const epsilon_0 = 8.8541878128e-12 #F m^-1

const gc_triangle = 2 # ghost cells, this must be changed if particle-to-grid shapes change

const c_bc_periodic = 100
const c_bc_open = 101

const c_stag_centre = 200 
const c_stag_right = 201

const c_field_electric = 300
const c_field_magnetic = 301
const c_field_rho = 302
const c_field_pot = 303

end