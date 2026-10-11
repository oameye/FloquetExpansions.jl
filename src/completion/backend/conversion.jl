function completion_scalar(x::SQA.CNum)::CompletionScalar
  value = SQA.to_num(SQA.trigonometric_form(x))::Complex{Symbolics.Num}
  return completion_scalar(value)
end

function completion_matrix(matrix::KossakowskiMatrix)::CompletionMatrix
  result = completion_matrix_zeros(size(matrix, 1), size(matrix, 2))
  for index in eachindex(matrix)
    result[index] = completion_scalar(matrix[index])
  end
  return result
end

function exact_numeric_radical_replacement(rational::Rational{Int})::Symbolics.Num
  reciprocal_root = Symbolics.Num(
    Symbolics.term(sqrt, Symbolics.unwrap(Symbolics.Num(inv(rational))); type=Real)
  )
  return Symbolics.Num(rational) * reciprocal_root
end

function collect_exact_numeric_radicals!(
  replacements::Dict{Symbolics.Num,Symbolics.Num}, node
)::Dict{Symbolics.Num,Symbolics.Num}
  Symbolics.iscall(node) || return replacements
  arguments = Symbolics.arguments(node)
  if Symbolics.operation(node) === sqrt && length(arguments) == 1
    radicand = Symbolics.unwrap_const(only(arguments))
    if radicand isa Union{Int,Rational{Int}} && radicand > 1
      replacements[Symbolics.Num(node)] = exact_numeric_radical_replacement(
        Rational{Int}(radicand)
      )
      return replacements
    end
  end
  for argument in arguments
    collect_exact_numeric_radicals!(replacements, argument)
  end
  return replacements
end

function preserve_exact_numeric_radicals(value::Symbolics.Num)::Symbolics.Num
  simplified = Symbolics.simplify(value)
  replacements = Dict{Symbolics.Num,Symbolics.Num}()
  collect_exact_numeric_radicals!(replacements, Symbolics.unwrap(simplified))
  isempty(replacements) && return simplified
  return Symbolics.substitute(simplified, replacements)
end

function coefficient_from_completion(value::CompletionScalar)::SQA.CNum
  simplified = simplify_scalar(value)
  exact_value = complex(
    preserve_exact_numeric_radicals(real(simplified)),
    preserve_exact_numeric_radicals(imag(simplified)),
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
    coherent = coherent + reattach(hamiltonian, drive_frequency, index - 1)
    push!(matrices, matrix)
  end
  return RetainedGKSLData(SQA.simplify(coherent), matrices)
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
      physical[entry] = simplify_coefficient(scale * matrix[entry])
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
