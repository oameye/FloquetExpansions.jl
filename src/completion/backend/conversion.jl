function completion_scalar(x::SQA.CNum)::CompletionScalar
  value = SQA.to_num(x)::Complex{Symbolics.Num}
  return completion_scalar(value)
end

function completion_matrix(matrix::KossakowskiMatrix)::CompletionMatrix
  result = completion_matrix_zeros(size(matrix, 1), size(matrix, 2))
  for index in eachindex(matrix)
    result[index] = completion_scalar(matrix[index])
  end
  return result
end

function coefficient_from_completion(value::CompletionScalar)::SQA.CNum
  simplified = simplify_scalar(value)
  exact_value = complex(
    Symbolics.simplify(real(simplified))::Symbolics.Num,
    Symbolics.simplify(imag(simplified))::Symbolics.Num,
  )
  return SQA.simplify(convert(SQA.CNum, exact_value))::SQA.CNum
end

function coefficient_matrix_from_completion(matrix::CompletionMatrix)::KossakowskiMatrix
  result = coefficient_matrix(size(matrix, 1), size(matrix, 2))
  for index in eachindex(matrix)
    result[index] = coefficient_from_completion(matrix[index])
  end
  return result
end

function condition_coefficients(values::Vector{CompletionScalar})::Vector{SQA.CNum}
  result = SQA.CNum[]
  sizehint!(result, length(values))
  for value in values
    push!(result, coefficient_from_completion(value))
  end
  return result
end

function completion_series(matrices::Vector{KossakowskiMatrix})::MatrixSeries
  result = CompletionMatrix[]
  sizehint!(result, length(matrices))
  for matrix in matrices
    push!(result, completion_matrix(matrix))
  end
  return result
end

@inline function inverse_drive_power(wd::Symbolics.Num, grade::Int)::Symbolics.Num
  return iszero(grade) ? Symbolics.Num(1) : (wd^(-grade))::Symbolics.Num
end

function retained_gksl_data(
  expansion::FloquetExpansion, frame::DissipativeFrame
)::RetainedGKSLData
  components = getfield(expansion, :effective_components)
  drive_frequency = getfield(expansion, :generator).wd
  coherent = zero(SQA.QAdd)
  matrices = KossakowskiMatrix[]
  sizehint!(matrices, length(components))
  for index in eachindex(components)
    hamiltonian, matrix = extract_gksl(components[index], frame)
    coherent = (coherent + reattach(hamiltonian, drive_frequency, index - 1))::SQA.QAdd
    push!(matrices, matrix)
  end
  return RetainedGKSLData(SQA.simplify(coherent)::SQA.QAdd, matrices)
end

function physical_kossakowski_series(
  matrices::Vector{KossakowskiMatrix}, wd::Symbolics.Num
)::Vector{KossakowskiMatrix}
  result = KossakowskiMatrix[]
  sizehint!(result, length(matrices))
  for index in eachindex(matrices)
    matrix = matrices[index]
    grade = index - 1
    scale = if iszero(grade)
      coefficient_one()
    else
      convert(SQA.CNum, inverse_drive_power(wd, grade))::SQA.CNum
    end
    physical = coefficient_matrix(size(matrix, 1), size(matrix, 2))
    for entry in eachindex(matrix)
      physical[entry] = simplify_coefficient((scale * matrix[entry])::SQA.CNum)
    end
    push!(result, physical)
  end
  return result
end

function seed_completion_conditions!(
  conditions::CompletionConditions, ::NoProvenance
)::CompletionConditions
  return conditions
end

function seed_completion_conditions!(
  conditions::CompletionConditions, provenance::MicroscopicProvenance
)::CompletionConditions
  for assumption in provenance.rate_assumptions
    require_positivity!(conditions, completion_scalar(assumption.rate))
  end
  return conditions
end
