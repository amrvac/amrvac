FOBJECTS += mod_rho_phys.t mod_rho_roe.t mod_rho.t

mod_rho.o: mod_physics.o mod_rho_phys.o mod_rho_roe.o
mod_rho_phys.o: mod_physics.o
mod_rho_roe.o: mod_rho_phys.o mod_physics_roe.o
