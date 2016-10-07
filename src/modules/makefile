#########################################################################
# Here we have some auxilary modules
#########################################################################
FOBJECTS += mod_constants$F   mod_timing$F   mod_particles$F  \
mod_gridvars$F  mod_odeint$F  mod_oneblock$F mod_lu$F mod_cooling$F

#########################################################################
# Dependencies:
#########################################################################
amrvacusr.o : mod_particles.o mod_odeint.o mod_cooling.o
mod_particles.o : mod_constants.o mod_gridvars.o
advance.o : mod_particles.o
amrvac.o : mod_gridvars.o
mod_metric.o : mod_lu.o
