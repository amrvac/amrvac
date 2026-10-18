SETUP_FLAGS := -d=2 -v=2

# beam is vacuum transport across the periodic seam, slab an absorbing slab that exercises the
# heating blend and taper, and table the Rosseland table opacity on a stratified atmosphere in cgs.
# gather is the replicated reference, wave the ray-by-ray KBA and pipe the pipelined KBA
CASES := beam_gather beam_wave beam_pipe slab_gather slab_pipe table_pipe

TESTS := $(CASES:%=sc_rt_2d_%.log)

include ../../test_rules.make

$(foreach c, $(CASES), $(eval sc_rt_2d_$(c).log: sc_rt_2d.par $(c).par))
