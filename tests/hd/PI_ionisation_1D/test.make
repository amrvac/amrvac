SETUP_FLAGS := -d=1 -v=1

# Partial-ionisation (eos_type='PI') smoke matrix: each backend table x closure.
#   chromosphere / no-energy   (Rfactor_from_temperature_ionization path)
#   chromosphere / energy      (eint<->T inversion + energy-mode csound2)
#   flare        / energy      (alternate T-only table)
#   prominence   / no-energy   (T,p pressure-bisection path)
# Each parfile is standalone (no scheme overlay).
TESTS := pi_chromo.log pi_chromoE.log pi_flareE.log pi_prom.log

include ../../test_rules.make

$(eval pi_chromo.log:  pi_chromo.par)
$(eval pi_chromoE.log: pi_chromoE.par)
$(eval pi_flareE.log:  pi_flareE.par)
$(eval pi_prom.log:    pi_prom.par)
