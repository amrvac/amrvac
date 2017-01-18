FOBJECTS += mod_hd_phys.t mod_hd_roe.t mod_hd_hllc.t mod_hd_ppm.t mod_hd.t

mod_hd.o: mod_physics.o mod_hd_phys.o mod_hd_hllc.o mod_hd_roe.o mod_hd_ppm.o #mod_dust.o
mod_hd_phys.o: mod_thermal_conduction.o
mod_hd_roe.o mod_hd_ppm.o mod_hd_hllc.o: mod_hd_phys.o mod_physics_roe.o
