using FloquetExpansions

const BlochResidualNode = FloquetExpansions.BlochResidualPlanNode
const BlochEvaluationPlanCounts = FloquetExpansions.BlochProjectionPlanCounts
const BlochEvaluationPlan = FloquetExpansions.BlochProjectionPlan

function compile_bloch_evaluation_plan(support_input, order::Int; zero_harmonic=nothing)
  if isnothing(zero_harmonic)
    return FloquetExpansions.compile_bloch_projection_plan(support_input, order)
  end
  return FloquetExpansions.compile_bloch_projection_plan(
    support_input, order, zero_harmonic
  )
end

function evaluate_bloch_evaluation_plan(args...; kwargs...)
  return FloquetExpansions.evaluate_bloch_projection_plan(args...; kwargs...)
end
