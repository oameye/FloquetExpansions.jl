using Test

const CKBernoulliExact = Complex{Rational{Int}}
const ck_bernoulli_im = CKBernoulliExact(0 // 1, 1 // 1)

struct CKBernoulliKernelTerm
  coefficient::CKBernoulliExact
  period_power::Int
  coordinate_power::Int
end

function ck_bernoulli_accumulate!(
  target::Dict{Int,CKBernoulliExact}, power::Int, value::CKBernoulliExact
)
  updated = get(target, power, zero(CKBernoulliExact)) + value
  if iszero(updated)
    haskey(target, power) && delete!(target, power)
  else
    target[power] = updated
  end
  return target
end

function ck_periodic_monomial_fourier_integral(degree::Int, harmonic::Int)
  degree >= 0 || throw(ArgumentError("degree must be nonnegative"))
  iszero(harmonic) && throw(ArgumentError("harmonic must be nonzero"))
  degree == 0 && return Dict{Int,CKBernoulliExact}()

  previous = ck_periodic_monomial_fourier_integral(degree - 1, harmonic)
  result = Dict{Int,CKBernoulliExact}()
  inverse = ck_bernoulli_im / harmonic
  ck_bernoulli_accumulate!(result, degree, inverse)
  for (period_power, coefficient) in previous
    ck_bernoulli_accumulate!(
      result, period_power, -(degree * inverse) * coefficient
    )
  end
  return result
end

function ck_bernoulli_kernel(order::Int)
  if order == 1
    return [
      CKBernoulliKernelTerm(ck_bernoulli_im / 2, 1, 0),
      CKBernoulliKernelTerm(-ck_bernoulli_im, 0, 1),
    ]
  elseif order == 2
    return [
      CKBernoulliKernelTerm(CKBernoulliExact(1 // 12), 2, 0),
      CKBernoulliKernelTerm(CKBernoulliExact(-1 // 2), 1, 1),
      CKBernoulliKernelTerm(CKBernoulliExact(1 // 2), 0, 2),
    ]
  elseif order == 3
    return [
      CKBernoulliKernelTerm(ck_bernoulli_im / 12, 2, 1),
      CKBernoulliKernelTerm(-ck_bernoulli_im / 4, 1, 2),
      CKBernoulliKernelTerm(ck_bernoulli_im / 6, 0, 3),
    ]
  elseif order == 4
    return [
      CKBernoulliKernelTerm(CKBernoulliExact(1 // 720), 4, 0),
      CKBernoulliKernelTerm(CKBernoulliExact(-1 // 24), 2, 2),
      CKBernoulliKernelTerm(CKBernoulliExact(1 // 12), 1, 3),
      CKBernoulliKernelTerm(CKBernoulliExact(-1 // 24), 0, 4),
    ]
  end
  throw(ArgumentError("reference kernels are tabulated through order four"))
end

function ck_bernoulli_fourier_coefficient(terms, harmonic::Int)
  result = Dict{Int,CKBernoulliExact}()
  for term in terms
    integral = ck_periodic_monomial_fourier_integral(term.coordinate_power, harmonic)
    for (integral_period_power, coefficient) in integral
      output_period_power = term.period_power + integral_period_power - 1
      ck_bernoulli_accumulate!(
        result, output_period_power, term.coefficient * coefficient
      )
    end
  end
  return result
end

function ck_bernoulli_average(terms)
  result = Dict{Int,CKBernoulliExact}()
  for term in terms
    output_period_power = term.period_power + term.coordinate_power
    coefficient = term.coefficient / (term.coordinate_power + 1)
    ck_bernoulli_accumulate!(result, output_period_power, coefficient)
  end
  return result
end

function ck_bernoulli_term_dictionary(terms)
  result = Dict{Tuple{Int,Int},CKBernoulliExact}()
  for term in terms
    key = (term.period_power, term.coordinate_power)
    updated = get(result, key, zero(CKBernoulliExact)) + term.coefficient
    if iszero(updated)
      haskey(result, key) && delete!(result, key)
    else
      result[key] = updated
    end
  end
  return result
end

function ck_bernoulli_derivative(terms)
  return CKBernoulliKernelTerm[
    CKBernoulliKernelTerm(
      term.coordinate_power * term.coefficient,
      term.period_power,
      term.coordinate_power - 1,
    ) for term in terms if term.coordinate_power > 0
  ]
end

function ck_bernoulli_scale(scale, terms)
  return CKBernoulliKernelTerm[
    CKBernoulliKernelTerm(scale * term.coefficient, term.period_power, term.coordinate_power) for
    term in terms
  ]
end

@testset "periodic Bernoulli kernels exactly quotient continuum homological denominators" begin
  harmonics = vcat(collect(-7:-1), collect(1:7))
  for order in 1:4
    kernel = ck_bernoulli_kernel(order)
    @test isempty(ck_bernoulli_average(kernel))

    for harmonic in harmonics
      coefficient = ck_bernoulli_fourier_coefficient(kernel, harmonic)
      expected = CKBernoulliExact(1 // (harmonic^order))
      @test coefficient == Dict(0 => expected)
    end
  end
end

@testset "repeated continuum homological solves remain finite polynomial kernels" begin
  for order in 1:3
    derivative = ck_bernoulli_term_dictionary(
      ck_bernoulli_derivative(ck_bernoulli_kernel(order + 1))
    )
    expected = ck_bernoulli_term_dictionary(
      ck_bernoulli_scale(ck_bernoulli_im, ck_bernoulli_kernel(order))
    )
    @test derivative == expected
  end

  first = ck_bernoulli_term_dictionary(ck_bernoulli_kernel(1))
  @test first == Dict((1, 0) => ck_bernoulli_im / 2, (0, 1) => -ck_bernoulli_im)
end
