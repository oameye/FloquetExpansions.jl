function completion_scalar(x::SQA.CNum)
  value = SQA.to_num(SQA.trigonometric_form(x))
  return completion_scalar(value)
end

function completion_matrix(matrix::KossakowskiMatrix)
  result = completion_matrix_zeros(size(matrix, 1), size(matrix, 2))
  for index in eachindex(matrix)
    result[index] = completion_scalar(matrix[index])
  end
  return result
end

function coefficient_from_completion(value::CompletionScalar)
  return SQA.simplify(convert(SQA.CNum, simplify_scalar(value)))
end

function coefficient_matrix_from_completion(matrix::CompletionMatrix)
  result = coefficient_matrix(size(matrix, 1), size(matrix, 2))
  for index in eachindex(matrix)
    result[index] = coefficient_from_completion(matrix[index])
  end
  return result
end

function condition_coefficients(values::Vector{CompletionScalar})
  result = SQA.CNum[]
  sizehint!(result, length(values))
  for value in values
    push!(result, coefficient_from_completion(value))
  end
  return result
end

function completion_series(matrices::Vector{KossakowskiMatrix})
  result = CompletionMatrix[]
  sizehint!(result, length(matrices))
  for matrix in matrices
    push!(result, completion_matrix(matrix))
  end
  return result
end

function retained_gksl_data(expansion::FloquetExpansion, frame::SQA.DissipativeFrame)
  components = getfield(expansion, :effective_components)
  drive_frequency = getfield(expansion, :generator).wd
  coherent = zero(SQA.QAdd)
  matrices = KossakowskiMatrix[]
  sizehint!(matrices, length(components))
  for index in eachindex(components)
    hamiltonian, matrix = SQA.gksl_coordinates(components[index], frame)
    coherent = coherent + reattach(hamiltonian, drive_frequency, index - 1)
    push!(matrices, matrix)
  end
  return RetainedGKSLData(SQA.simplify(coherent), matrices)
end

function physical_kossakowski_series(matrices::Vector{KossakowskiMatrix}, wd::Symbolics.Num)
  result = KossakowskiMatrix[]
  sizehint!(result, length(matrices))
  for index in eachindex(matrices)
    matrix = matrices[index]
    grade = index - 1
    scale = if iszero(grade)
      coefficient_one()
    else
      convert(SQA.CNum, inverse_drive_power(wd, grade))
    end
    physical = coefficient_matrix(size(matrix, 1), size(matrix, 2))
    for entry in eachindex(matrix)
      physical[entry] = simplify_coefficient(scale * matrix[entry])
    end
    push!(result, physical)
  end
  return result
end

function seed_completion_conditions!(conditions::CompletionConditions, ::NoProvenance)
  return conditions
end

function seed_completion_conditions!(
  conditions::CompletionConditions, provenance::MicroscopicProvenance
)
  for assumption in provenance.rate_assumptions
    require_positivity!(conditions, completion_scalar(assumption.rate))
  end
  return conditions
end
