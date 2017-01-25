FOBJECTS += mod_mhd_phys.t mod_mhd_roe.t mod_mhd_hllc.t mod_mhd_ppm.t mod_mhd.t

mod_mhd.o: mod_physics.o mod_mhd_phys.o mod_mhd_hllc.o mod_mhd_roe.o mod_mhd_ppm.o
mod_mhd_phys.o: mod_thermal_conduction.o mod_radiative_cooling.o
mod_mhd_roe.o mod_mhd_ppm.o mod_mhd_hllc.o: mod_mhd_phys.o mod_physics_roe.o
