using FloquetExpansions
using JET: JET

jumps = Dict(
  FloquetExpansions.ck_jump_vertex(1, -1) => 1.0,
  FloquetExpansions.ck_jump_vertex(1, 0) => 2.0,
  FloquetExpansions.ck_jump_vertex(1, 1) => 3.0,
)
drifts = Dict(FloquetExpansions.ck_drift_vertex(0) => -0.5)
A1 = FloquetExpansions.ck_kernel_generator(jumps, 0.0)
A2 = FloquetExpansions.ck_kernel_generator(drifts, 0.0)
zero_state = FloquetExpansions.ck_kernel_zero(0, 0.0)

JET.@test_opt target_modules=(FloquetExpansions,) FloquetExpansions.evaluate_ck_hori_deprit(
  [A1, A2], 3, zero_state
)
