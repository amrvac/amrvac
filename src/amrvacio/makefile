#########################################################################
# Everything IO-related goes here:
#########################################################################
INCLUDES += mod_config.t
FOBJECTS += mod_input_output.t convert.t kracken.t analysis.t collapse.t	\
slice.t

#########################################################################
# Dependencies:
#########################################################################
mod_input_output.o : kracken.o mod_forest.o
analysis.o : mod_constants.o
