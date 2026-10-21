FOBJECTS += mod_eos_container.t mod_eos_shared_functions.t mod_eos_interp.t mod_eos_lte_tables.t mod_eos_lte_saha.t mod_eos_lte_entropy.t mod_eos_lte_state.t mod_eos_pi_tables.t mod_eos_fi.t mod_eos_pi.t mod_eos_lte.t mod_eos.t

# Explicit dependencies (otherwise parallel make can race on .mod files)
mod_eos_shared_functions.o: mod_eos_container.mod
mod_eos_interp.o: mod_eos_container.mod
mod_eos_lte_tables.o: mod_eos_container.mod mod_eos_interp.mod
mod_eos_lte_saha.o: mod_eos_container.mod mod_eos_interp.mod mod_eos_lte_tables.mod
mod_eos_lte_entropy.o: mod_eos_container.mod mod_eos_lte_tables.mod
mod_eos_pi_tables.o: mod_comm_lib.mod mod_constants.mod mod_global_parameters.mod
mod_eos_fi.o: mod_eos_container.mod
mod_eos_pi.o: mod_eos_container.mod mod_eos_pi_tables.mod mod_eos_fi.mod
mod_eos_lte.o: mod_eos_container.mod mod_eos_shared_functions.mod mod_eos_interp.mod mod_eos_lte_saha.mod mod_eos_lte_entropy.mod mod_eos_lte_tables.mod
mod_eos.o: mod_eos_container.mod mod_eos_shared_functions.mod mod_eos_interp.mod mod_eos_lte_saha.mod mod_eos_lte_entropy.mod mod_eos_lte_tables.mod mod_eos_fi.mod mod_eos_pi.mod mod_eos_lte.mod
