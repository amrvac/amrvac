#########################################################################
# Everything IO-related goes here:
#########################################################################
FOBJECTS += amrio$F  convert$F  kracken$F  analysis$F collapse$F slice$F
OBJECTS += amrio.o convert.o kracken.o  analysis.o collapse.o slice.o

#########################################################################
# Dependencies:
#########################################################################
amrio.o : kracken.o mod_forest.o
analysis.o : mod_constants.o