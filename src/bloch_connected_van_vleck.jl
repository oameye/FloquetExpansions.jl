struct BlochLyndonBracketNode
  left::Int
  right::Int
end

struct BlochLyndonOutputTerm{C}
  node::Int
  coefficient::C
end

struct BlochConnectedLogPlan{H,C}
  leaves::Vector{H}
  brackets::Vector{BlochLyndonBracketNode}
  outputs::Vector{Dict{H,Vector{BlochLyndonOutputTerm{C}}}}
end

struct BlochStaticExpDependency{H}
  left_order::Int
  left_harmonic::H
  right_order::Int
  right_harmonic::H
  right_node::Int
end

struct BlochStaticExpNode{H}
  power::Int
  order::Int
  harmonic::H
  dependencies::Vector{BlochStaticExpDependency{H}}
end

struct BlochStaticExpPlan{H}
  zero_harmonic::H
  nodes::Vector{BlochStaticExpNode{H}}
  targets::Vector{Vector{Int}}
  product_count::Int
end

struct BlochConnectedVanVleckPlan{H,C}
  log::BlochConnectedLogPlan{H,C}
  static::BlochStaticExpPlan{H}
end

struct BlochConnectedVanVleckResult{H,T}
  static_factor::Vector{T}
  inverse_static_factor::Vector{T}
  log_embedding::Vector{Dict{H,T}}
  effective::Vector{T}
  counts::BlochVanVleckCounts
end

function bloch_word_accumulate!(terms::Dict{Tuple,C}, word::Tuple, coefficient) where {C}
  value = get(terms, word, zero(C)) + convert(C, coefficient)
  if iszero(value)
    haskey(terms, word) && delete!(terms, word)
  else
    terms[word] = value
  end
  return terms
end

function bloch_word_terms!(embedding::Dict{H,Dict{Tuple,C}}, harmonic::H) where {H,C}
  return get!(embedding, harmonic) do
    return Dict{Tuple,C}()
  end
end

function bloch_word_add_terms!(
  destination::Dict{Tuple,C}, source::Dict{Tuple,C}, weight::C=one(C)
) where {C}
  for (word, coefficient) in source
    bloch_word_accumulate!(destination, word, weight * coefficient)
  end
  return destination
end

function bloch_word_copy_embedding(embedding::Dict{H,Dict{Tuple,C}}) where {H,C}
  return Dict(harmonic => copy(terms) for (harmonic, terms) in embedding)
end

function bloch_word_add_embedding!(
  destination::Dict{H,Dict{Tuple,C}},
  source::Dict{H,Dict{Tuple,C}},
  weight::C=one(C),
) where {H,C}
  for (harmonic, terms) in source
    bloch_word_add_terms!(bloch_word_terms!(destination, harmonic), terms, weight)
  end
  return destination
end

function bloch_word_periodic_product(
  left::Dict{H,Dict{Tuple,C}}, right::Dict{H,Dict{Tuple,C}}
) where {H,C}
  result = Dict{H,Dict{Tuple,C}}()
  for (left_harmonic, left_terms) in left,
    (right_harmonic, right_terms) in right

    output = bloch_word_terms!(result, left_harmonic + right_harmonic)
    for (left_word, left_coefficient) in left_terms,
      (right_word, right_coefficient) in right_terms

      bloch_word_accumulate!(
        output,
        (left_word..., right_word...),
        left_coefficient * right_coefficient,
      )
    end
  end
  return result
end

function bloch_word_right_static_product(
  periodic::Dict{H,Dict{Tuple,C}}, static::Dict{Tuple,C}
) where {H,C}
  result = Dict{H,Dict{Tuple,C}}()
  for (harmonic, periodic_terms) in periodic
    output = bloch_word_terms!(result, harmonic)
    for (periodic_word, periodic_coefficient) in periodic_terms,
      (static_word, static_coefficient) in static

      bloch_word_accumulate!(
        output,
        (periodic_word..., static_word...),
        periodic_coefficient * static_coefficient,
      )
    end
  end
  return result
end

function bloch_word_scale_embedding!(
  embedding::Dict{H,Dict{Tuple,C}}, weight::C
) where {H,C}
  for terms in values(embedding)
    for word in collect(keys(terms))
      coefficient = weight * terms[word]
      if iszero(coefficient)
        delete!(terms, word)
      else
        terms[word] = coefficient
      end
    end
  end
  return embedding
end

function bloch_compile_connected_log_words(
  support_input, order::Int, zero_harmonic::H
) where {H}
  order >= 1 || throw(ArgumentError("order must be >= 1"))
  support = sort(unique(H[harmonic for harmonic in support_input]))
  isempty(support) && throw(ArgumentError("harmonic support must not be empty"))

  C = Rational{Int}
  one_coefficient = one(C)
  effective = Vector{Dict{Tuple,C}}()
  B0 = Dict{Tuple,C}()
  if zero_harmonic in support
    B0[(zero_harmonic,)] = one_coefficient
  end
  push!(effective, B0)

  wave = Vector{Dict{H,Dict{Tuple,C}}}()
  if order > 1
    X1 = Dict{H,Dict{Tuple,C}}()
    for harmonic in support
      harmonic == zero_harmonic && continue
      bloch_word_terms!(X1, harmonic)[(harmonic,)] = 1 // harmonic
    end
    push!(wave, X1)
  end

  for n in 1:(order - 1)
    residual = Dict{H,Dict{Tuple,C}}()

    for generator_harmonic in support, (wave_harmonic, wave_terms) in wave[n]
      output = bloch_word_terms!(residual, generator_harmonic + wave_harmonic)
      for (wave_word, wave_coefficient) in wave_terms
        bloch_word_accumulate!(
          output,
          (generator_harmonic, wave_word...),
          wave_coefficient,
        )
      end
    end

    for wave_order in 1:n
      effective_order = n - wave_order
      B = effective[effective_order + 1]
      isempty(B) && continue
      for (harmonic, wave_terms) in wave[wave_order]
        output = bloch_word_terms!(residual, harmonic)
        for (wave_word, wave_coefficient) in wave_terms,
          (effective_word, effective_coefficient) in B

          bloch_word_accumulate!(
            output,
            (wave_word..., effective_word...),
            -wave_coefficient * effective_coefficient,
          )
        end
      end
    end

    push!(effective, copy(get(residual, zero_harmonic, Dict{Tuple,C}())))

    if n < order - 1
      Xnext = Dict{H,Dict{Tuple,C}}()
      for (harmonic, residual_terms) in residual
        harmonic == zero_harmonic && continue
        output = bloch_word_terms!(Xnext, harmonic)
        weight = 1 // harmonic
        for (word, coefficient) in residual_terms
          bloch_word_accumulate!(output, word, weight * coefficient)
        end
      end
      push!(wave, Xnext)
    end
  end

  canonical_order = order - 1
  static_factor = Dict{Tuple,C}[Dict(() => one_coefficient)]
  normalized = Vector{Dict{H,Dict{Tuple,C}}}()
  log_embedding = Vector{Dict{H,Dict{Tuple,C}}}()
  powers = [
    [Dict{H,Dict{Tuple,C}}() for _ in 1:max(canonical_order, 1)] for
    _ in 1:max(canonical_order, 1)
  ]

  for n in 1:canonical_order
    prefactor = bloch_word_copy_embedding(wave[n])
    for j in 1:(n - 1)
      correction = bloch_word_right_static_product(wave[j], static_factor[n - j + 1])
      bloch_word_add_embedding!(prefactor, correction)
    end

    nonlinear_log = Dict{H,Dict{Tuple,C}}()
    for power in 2:n
      power_coefficient = Dict{H,Dict{Tuple,C}}()
      for k in 1:(n - power + 1)
        contribution = bloch_word_periodic_product(normalized[k], powers[power - 1][n - k])
        bloch_word_add_embedding!(power_coefficient, contribution)
      end
      powers[power][n] = power_coefficient
      weight = convert(C, (-1)^(power + 1) * (1 // power))
      contribution = bloch_word_copy_embedding(power_coefficient)
      bloch_word_scale_embedding!(contribution, weight)
      bloch_word_add_embedding!(nonlinear_log, contribution)
    end

    candidate = bloch_word_copy_embedding(prefactor)
    bloch_word_add_embedding!(candidate, nonlinear_log)
    static_n = Dict{Tuple,C}()
    for (word, coefficient) in get(candidate, zero_harmonic, Dict{Tuple,C}())
      bloch_word_accumulate!(static_n, word, -coefficient)
    end
    push!(static_factor, static_n)

    normalized_n = bloch_word_copy_embedding(prefactor)
    bloch_word_add_terms!(bloch_word_terms!(normalized_n, zero_harmonic), static_n)
    push!(normalized, normalized_n)
    powers[1][n] = normalized_n

    log_n = bloch_word_copy_embedding(normalized_n)
    bloch_word_add_embedding!(log_n, nonlinear_log)
    push!(log_embedding, log_n)
  end

  return log_embedding
end

function bloch_is_lyndon_word(word::Tuple)
  isempty(word) && return false
  return all(isless(word, word[index:end]) for index in 2:length(word))
end

function bloch_lyndon_standard_factorization(word::Tuple)
  length(word) > 1 || throw(ArgumentError("a Lyndon leaf has no standard factorization"))
  for split in 2:length(word)
    suffix = word[split:end]
    bloch_is_lyndon_word(suffix) && return word[1:(split - 1)], suffix
  end
  return throw(ArgumentError("word is not Lyndon"))
end

function bloch_word_leaf(harmonic, ::Type{C}) where {C}
  return Dict{Tuple,C}((harmonic,) => one(C))
end

function bloch_word_commutator(left::Dict{Tuple,C}, right::Dict{Tuple,C}) where {C}
  result = Dict{Tuple,C}()
  for (left_word, left_coefficient) in left,
    (right_word, right_coefficient) in right

    bloch_word_accumulate!(
      result, (left_word..., right_word...), left_coefficient * right_coefficient
    )
    bloch_word_accumulate!(
      result, (right_word..., left_word...), -right_coefficient * left_coefficient
    )
  end
  return result
end

function bloch_lyndon_bracket_words!(
  cache::Dict{Tuple,Dict{Tuple,C}}, word::Tuple, ::Type{C}
) where {C}
  haskey(cache, word) && return cache[word]
  result = if length(word) == 1
    bloch_word_leaf(first(word), C)
  else
    left_word, right_word = bloch_lyndon_standard_factorization(word)
    left = bloch_lyndon_bracket_words!(cache, left_word, C)
    right = bloch_lyndon_bracket_words!(cache, right_word, C)
    bloch_word_commutator(left, right)
  end
  cache[word] = result
  return result
end

function bloch_lyndon_decomposition(
  word_terms::Dict{Tuple,C}, bracket_cache::Dict{Tuple,Dict{Tuple,C}}
) where {C}
  residual = copy(word_terms)
  coefficients = Dict{Tuple,C}()
  while !isempty(residual)
    word = minimum(keys(residual))
    bloch_is_lyndon_word(word) ||
      throw(ArgumentError("primitive polynomial has a non-Lyndon leading word"))
    coefficient = residual[word]
    coefficients[word] = get(coefficients, word, zero(C)) + coefficient
    bracket = bloch_lyndon_bracket_words!(bracket_cache, word, C)
    for (term_word, term_coefficient) in bracket
      updated = get(residual, term_word, zero(C)) - coefficient * term_coefficient
      if iszero(updated)
        delete!(residual, term_word)
      else
        residual[term_word] = updated
      end
    end
  end
  return coefficients
end

function bloch_lyndon_ensure_node!(
  word::Tuple,
  leaves::Vector{H},
  leaf_ids::Dict{H,Int},
  word_ids::Dict{Tuple,Int},
  brackets::Vector{BlochLyndonBracketNode},
) where {H}
  length(word) == 1 && return leaf_ids[first(word)]
  haskey(word_ids, word) && return word_ids[word]
  left_word, right_word = bloch_lyndon_standard_factorization(word)
  left = bloch_lyndon_ensure_node!(left_word, leaves, leaf_ids, word_ids, brackets)
  right = bloch_lyndon_ensure_node!(right_word, leaves, leaf_ids, word_ids, brackets)
  push!(brackets, BlochLyndonBracketNode(left, right))
  node = length(leaves) + length(brackets)
  word_ids[word] = node
  return node
end

function compile_bloch_connected_log_plan(
  support_input, order::Int, zero_harmonic::H
) where {H}
  C = Rational{Int}
  leaves = sort(unique(H[harmonic for harmonic in support_input]))
  isempty(leaves) && throw(ArgumentError("harmonic support must not be empty"))
  leaf_ids = Dict(harmonic => index for (index, harmonic) in enumerate(leaves))
  word_ids = Dict{Tuple,Int}()
  brackets = BlochLyndonBracketNode[]
  bracket_cache = Dict{Tuple,Dict{Tuple,C}}()
  log_embedding = bloch_compile_connected_log_words(leaves, order, zero_harmonic)
  outputs = Vector{Dict{H,Vector{BlochLyndonOutputTerm{C}}}}()

  for embedding in log_embedding
    compiled = Dict{H,Vector{BlochLyndonOutputTerm{C}}}()
    for (harmonic, word_terms) in embedding
      terms = BlochLyndonOutputTerm{C}[]
      for (word, coefficient) in bloch_lyndon_decomposition(word_terms, bracket_cache)
        node = bloch_lyndon_ensure_node!(word, leaves, leaf_ids, word_ids, brackets)
        push!(terms, BlochLyndonOutputTerm(node, coefficient))
      end
      compiled[harmonic] = terms
    end
    push!(outputs, compiled)
  end

  return BlochConnectedLogPlan(leaves, brackets, outputs)
end

function evaluate_bloch_connected_log_plan(
  plan::BlochConnectedLogPlan{H,C},
  components::AbstractDict{H,T};
  product,
  zero_component::T,
  simplifier=identity,
) where {H,C,T}
  nleaves = length(plan.leaves)
  values = Vector{T}(undef, nleaves + length(plan.brackets))
  for (index, harmonic) in enumerate(plan.leaves)
    values[index] = components[harmonic]
  end

  for (index, bracket) in enumerate(plan.brackets)
    left = values[bracket.left]
    right = values[bracket.right]
    values[nleaves + index] = simplifier(product(left, right) - product(right, left))::T
  end

  output = Vector{Dict{H,T}}()
  for embedding in plan.outputs
    evaluated = Dict{H,T}()
    for (harmonic, terms) in embedding
      value = zero_component
      for term in terms
        value += term.coefficient * values[term.node]
      end
      component = simplifier(value)::T
      iszero(component) || (evaluated[harmonic] = component)
    end
    push!(output, evaluated)
  end
  return output
end

bloch_connected_log_products(plan::BlochConnectedLogPlan) = 2 * length(plan.brackets)

function bloch_static_power_support(plan::BlochConnectedLogPlan{H}) where {H}
  order = length(plan.outputs)
  supports = [[Set{H}() for _ in 1:order] for _ in 1:order]
  for n in 1:order
    union!(supports[1][n], keys(plan.outputs[n]))
  end
  for power in 2:order, n in power:order
    for k in 1:(n - power + 1)
      for left in supports[1][k], right in supports[power - 1][n - k]
        push!(supports[power][n], left + right)
      end
    end
  end
  return supports
end

function compile_bloch_static_exp_plan(
  log_plan::BlochConnectedLogPlan{H}, zero_harmonic::H
) where {H}
  order = length(log_plan.outputs)
  supports = bloch_static_power_support(log_plan)
  nodes = BlochStaticExpNode{H}[]
  node_ids = Dict{Tuple{Int,Int,H},Int}()
  product_count = Ref(0)

  function ensure_node(power::Int, n::Int, harmonic::H)
    state = (power, n, harmonic)
    haskey(node_ids, state) && return node_ids[state]
    dependencies = BlochStaticExpDependency{H}[]
    for k in 1:(n - power + 1)
      right_order = n - k
      for left_harmonic in supports[1][k]
        for right_harmonic in supports[power - 1][right_order]
          left_harmonic + right_harmonic == harmonic || continue
          right_node =
            power == 2 ? 0 : ensure_node(power - 1, right_order, right_harmonic)
          push!(
            dependencies,
            BlochStaticExpDependency(
              k, left_harmonic, right_order, right_harmonic, right_node
            ),
          )
          product_count[] += 1
        end
      end
    end
    push!(nodes, BlochStaticExpNode(power, n, harmonic, dependencies))
    node = length(nodes)
    node_ids[state] = node
    return node
  end

  targets = [zeros(Int, order) for _ in 1:order]
  for n in 2:order, power in 2:n
    zero_harmonic in supports[power][n] || continue
    targets[n][power] = ensure_node(power, n, zero_harmonic)
  end
  return BlochStaticExpPlan(zero_harmonic, nodes, targets, product_count[])
end

function evaluate_bloch_static_exp_plan(
  plan::BlochStaticExpPlan{H},
  log_embedding::Vector{Dict{H,T}},
  identity_component::T,
  zero_component::T;
  product,
  simplifier=identity,
) where {H,T}
  values = Vector{T}(undef, length(plan.nodes))
  for (node_id, node) in enumerate(plan.nodes)
    value = zero_component
    for dependency in node.dependencies
      left = get(
        log_embedding[dependency.left_order], dependency.left_harmonic, zero_component
      )
      right = if dependency.right_node == 0
        get(log_embedding[dependency.right_order], dependency.right_harmonic, zero_component)
      else
        values[dependency.right_node]
      end
      value += product(left, right)
    end
    values[node_id] = simplifier(value)::T
  end

  static_factor = T[identity_component]
  for n in eachindex(log_embedding)
    value = get(log_embedding[n], plan.zero_harmonic, zero_component)
    for power in 2:n
      node = plan.targets[n][power]
      node == 0 && continue
      value += (1 // factorial(power)) * values[node]
    end
    push!(static_factor, simplifier(value)::T)
  end
  return static_factor
end

function compile_bloch_connected_van_vleck_plan(plan::BlochProjectionPlan{H}) where {H}
  log_plan = compile_bloch_connected_log_plan(plan.input_support, plan.order, plan.zero_harmonic)
  static_plan = compile_bloch_static_exp_plan(log_plan, plan.zero_harmonic)
  return BlochConnectedVanVleckPlan(log_plan, static_plan)
end

function evaluate_bloch_connected_van_vleck(
  projection_plan::BlochProjectionPlan{H},
  bloch::BlochProjectionResult{H,T},
  components::AbstractDict{H,T},
  connected_plan::BlochConnectedVanVleckPlan{H};
  product,
  zero_component::T,
  simplifier=identity,
) where {H,T}
  order = length(bloch.effective) - 1
  projection_plan.order == order + 1 ||
    throw(ArgumentError("Bloch plan/result truncations are inconsistent"))
  length(connected_plan.log.outputs) == order ||
    throw(ArgumentError("connected-log plan/result truncations are inconsistent"))
  connected_plan.static.zero_harmonic == projection_plan.zero_harmonic ||
    throw(ArgumentError("static-exponential and Bloch zero harmonics are inconsistent"))

  log_embedding = evaluate_bloch_connected_log_plan(
    connected_plan.log, components; product, zero_component, simplifier
  )
  identity_component = one(first(bloch.effective))
  static_factor = evaluate_bloch_static_exp_plan(
    connected_plan.static,
    log_embedding,
    identity_component,
    zero_component;
    product,
    simplifier,
  )
  counts = BlochVanVleckCounts()
  counts.harmonic_products =
    bloch_connected_log_products(connected_plan.log) + connected_plan.static.product_count
  inverse_static_factor = bloch_vv_static_series_inverse(
    static_factor, order, product, counts; simplifier
  )
  right_transformed = bloch_vv_static_series_product(
    bloch.effective, static_factor, order, product, counts; simplifier
  )
  effective = bloch_vv_static_series_product(
    inverse_static_factor, right_transformed, order, product, counts; simplifier
  )

  return BlochConnectedVanVleckResult(
    static_factor, inverse_static_factor, log_embedding, effective, counts
  )
end
