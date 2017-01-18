FOBJECTS += mod_physics.t mod_physics_roe.t mod_physics_hllc.t mod_physics_ppm.t

# FOBJECTS += mod_dust.t
FOBJECTS += mod_thermal_conduction.t

mod_physics.o: mod_physics_roe.o mod_physics_hllc.o mod_physics_ppm.o
mod_thermal_conduction.o: mod_ghostcells_update.o


