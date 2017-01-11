FOBJECTS += mod_rho_phys$F mod_rho_roe$F mod_rho$F
FOBJECTS += mod_hd_phys$F mod_hd_roe$F mod_hd_hllc$F mod_hd_ppm$F mod_hd$F
FOBJECTS += mod_physics$F mod_physics_roe$F mod_physics_hllc$F	\
mod_physics_ppm$F

# FOBJECTS += mod_dust$F
FOBJECTS += mod_thermal_conduction$F

mod_physics.o: mod_physics_roe.o mod_physics_hllc.o mod_physics_ppm.o

mod_rho.o: mod_physics.o mod_rho_phys.o mod_rho_roe.o
mod_rho_roe.o: mod_rho_phys.o mod_physics_roe.o

mod_hd.o: mod_physics.o mod_hd_phys.o mod_hd_hllc.o mod_hd_roe.o mod_hd_ppm.o #mod_dust.o
# mod_hd_phys.o: mod_dust.o
mod_hd_phys.o: mod_thermal_conduction.o
mod_hd_roe.o: mod_hd_phys.o mod_physics_roe.o
mod_hd_ppm.o: mod_hd_phys.o mod_physics_ppm.o
