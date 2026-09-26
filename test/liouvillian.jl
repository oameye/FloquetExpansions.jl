using Test
using FloquetExpansions
using Symbolics: Num, Symbolics, @variables
using LinearAlgebra: I, kron

include(joinpath(@__DIR__, "helpers", "shared.jl"))

space = FockSpace(:liouvillian)
a = Destroy(space, :a)
@variables γ::Real δ::Real ω::Real t::Real t₀::Real

function coeff_to_complex(coeff, substitutions::AbstractDict)
  dummy = coeff * one(SQA.QAdd)
  return tomatrix(dummy, 1, substitutions)[1, 1]
end

function superoperator_matrix(L::Liouvillian, space, d::Int, substitutions::AbstractDict)
  matrix = zeros(ComplexF64, d^2, d^2)
  for (left, right, coeff) in terms(L)
    left_mat = tomatrix(left, d, substitutions)
    right_mat = tomatrix(right, d, substitutions)
    c = coeff_to_complex(coeff, substitutions)
    iszero(c) && continue
    matrix .+= c .* kron(transpose(right_mat), left_mat)
  end
  return matrix
end

@testset "Liouvillian terms expose semantic triples" begin
  L = hamiltonian_action(a)
  observed = collect(terms(L))

  @test @inferred(collect(terms(L))) isa Vector{Tuple{SQA.QAdd,SQA.QAdd,SQA.CNum}}
  @test length(observed) == 2
  @test all(length(term) == 3 for term in observed)
  a_q = a + zero(one(a))
  minus_im = convert(SQA.CNum, -im)
  plus_im = convert(SQA.CNum, im)
  @test Set(observed) == Set(((a_q, one(a), minus_im), (one(a), a_q, plus_im)))
end

function independent_dissipator_matrix(C::Matrix{ComplexF64})
  norm = C' * C
  identity = Matrix{ComplexF64}(I, size(C, 1), size(C, 2))
  return kron(conj(C), C) - 0.5 * (kron(identity, norm) + kron(transpose(norm), identity))
end

function independent_liouvillian_matrix(H, C, J, rate)
  identity = Matrix{ComplexF64}(I, size(H, 1), size(H, 2))
  coherent = -im * (kron(identity, H) - kron(transpose(H), identity))
  return coherent +
         independent_dissipator_matrix(C) +
         rate * independent_dissipator_matrix(J)
end

@testset "Liouvillian channel constructors" begin
  H = a' * a
  coherent = hamiltonian_action(H)
  complete = dissipator(a)
  weighted = γ * dissipator(a')

  L = liouvillian(H; channels=(collapse(a), jump(a', γ)))

  @test L == coherent + complete + weighted
  @test L isa Liouvillian
  @test @inferred(hamiltonian_action(H)) isa Liouvillian
  @test @inferred(dissipator(a)) isa Liouvillian
  @test @inferred(γ * dissipator(a')) isa Liouvillian
  @test @inferred(liouvillian(H; channels=(jump(a, γ),))) isa Liouvillian
end

@testset "channel displays preserve repr and render physical terms" begin
  collapse_channel = collapse(a)
  jump_channel = jump(a, γ)
  one_collapse = [collapse_channel]
  two_collapses = [collapse_channel, collapse_channel]
  two_jumps = [jump_channel, jump_channel]
  mixed_channels = Union{typeof(collapse_channel),typeof(jump_channel)}[
    collapse_channel, jump_channel
  ]
  no_collapses = one_collapse[1:0]
  five_collapses = fill(collapse_channel, 5)

  operator_text = sprint(show, a)
  rate_text = sprint(show, γ)
  collapse_repr = "collapse($operator_text)"
  jump_repr = "jump($operator_text, $rate_text)"
  collapse_text = "𝒟[$operator_text]"
  jump_text = "($rate_text)𝒟[$operator_text]"

  @test sprint(show, collapse_channel) == collapse_repr
  @test sprint(show, jump_channel) == jump_repr
  @test sprint(show, MIME"text/plain"(), collapse_channel) == collapse_text
  @test sprint(show, MIME"text/plain"(), jump_channel) == jump_text
  @test sprint(show, MIME"text/plain"(), one_collapse) ==
    "1 collapse channel:\n  1: $collapse_text"
  @test sprint(show, MIME"text/plain"(), two_collapses) ==
    "2 collapse channels:\n  1: $collapse_text\n  2: $collapse_text"
  @test sprint(show, MIME"text/plain"(), two_jumps) ==
    "2 rate-weighted jump channels:\n  1: $jump_text\n  2: $jump_text"
  @test sprint(show, MIME"text/plain"(), mixed_channels) ==
    "2 channels:\n  1: $collapse_text\n  2: $jump_text"
  @test sprint(show, MIME"text/plain"(), no_collapses) == "0 collapse channels"
  full_five_collapses =
    "5 collapse channels:" * join("\n  $index: $collapse_text" for index in 1:5)
  @test sprint(show, MIME"text/plain"(), five_collapses) == full_five_collapses

  function limited_channel_display(mime, value; rows::Int=24)
    buffer = IOBuffer()
    io = IOContext(buffer, :limit => true, :displaysize => (rows, 80))
    show(io, mime, value)
    return String(take!(buffer))
  end

  limited_five_collapses = join(
    [
      "5 collapse channels:",
      "  1: $collapse_text",
      "  2: $collapse_text",
      "  ⋮ 1 channel omitted",
      "  4: $collapse_text",
      "  5: $collapse_text",
    ],
    "\n",
  )
  @test limited_channel_display(MIME"text/plain"(), five_collapses) ==
    limited_five_collapses
  @test limited_channel_display(MIME"text/plain"(), five_collapses; rows=4) ==
    "5 collapse channels:\n  1: $collapse_text\n  ⋮ 4 channels omitted"

  operator_latex = FloquetExpansions.latex_fragment(a)
  rate_latex = FloquetExpansions.latex_fragment(γ)
  collapse_latex = raw"\mathcal{D}\!\left[" * operator_latex * raw"\right]"
  jump_latex =
    raw"\left(" *
    rate_latex *
    raw"\right)\mathcal{D}\!\left[" *
    operator_latex *
    raw"\right]"

  @test showable(MIME"text/latex"(), collapse_channel)
  @test showable(MIME"text/latex"(), jump_channel)
  @test showable(MIME"text/latex"(), two_collapses)
  @test sprint(show, MIME"text/latex"(), collapse_channel) ==
    raw"\[" * collapse_latex * raw"\]"
  @test sprint(show, MIME"text/latex"(), jump_channel) == raw"\[" * jump_latex * raw"\]"
  @test sprint(show, MIME"text/latex"(), two_collapses) ==
    raw"\[\begin{aligned}\mathcal{L}_{\mathrm{diss}} &= " *
        collapse_latex *
        raw"\\&+ " *
        collapse_latex *
        raw"\end{aligned}\]"
  @test sprint(show, MIME"text/latex"(), two_jumps) ==
    raw"\[\begin{aligned}\mathcal{L}_{\mathrm{diss}} &= " *
        jump_latex *
        raw"\\&+ " *
        jump_latex *
        raw"\end{aligned}\]"
  @test limited_channel_display(MIME"text/latex"(), five_collapses) ==
    raw"\[\begin{aligned}\mathcal{L}_{\mathrm{diss}} &= " *
        collapse_latex *
        raw"\\&+ " *
        collapse_latex *
        raw"\\&\quad\vdots\quad\text{(1 channel omitted)}\\&+ " *
        collapse_latex *
        raw"\\&+ " *
        collapse_latex *
        raw"\end{aligned}\]"
  @test sprint(show, MIME"text/latex"(), no_collapses) ==
    raw"\[\mathcal{L}_{\mathrm{diss}} = 0\]"

  # Nested operator and rate LaTeX must be embedded as delimiterless fragments. A renderer strips
  # only the outermost delimiter pair, so an interior `$` or `equation` environment fails to parse.
  for rendered in (
    sprint(show, MIME"text/latex"(), collapse_channel),
    sprint(show, MIME"text/latex"(), jump_channel),
    sprint(show, MIME"text/latex"(), two_jumps),
    limited_channel_display(MIME"text/latex"(), five_collapses),
  )
    @test !occursin('$', rendered)
    @test !occursin("begin{equation}", rendered)
  end
end

@testset "jump rates are real and nonnegative by physical assumption" begin
  zero_H = zero(SQA.QAdd)

  @test liouvillian(zero_H; channels=(jump(a, 0),)) == zero(Liouvillian)
  @test liouvillian(zero_H; channels=(jump(a, 2),)) == 2 * dissipator(a)
  @test liouvillian(zero_H; channels=(jump(a, 2 + 0im),)) == 2 * dissipator(a)
  @test liouvillian(zero_H; channels=(jump(a, γ),)) == γ * dissipator(a)
  @test liouvillian(zero_H; channels=(jump(a, 1 + cos(ω * t)),)) ==
    (1 + cos(ω * t)) * dissipator(a)

  phase_rate = 1 + expim(ω * t) + expim(-ω * t)
  phase_liouvillian = liouvillian(zero_H; channels=(jump(a, phase_rate),))
  @test phase_liouvillian == phase_rate * dissipator(a)
  @test sort!(collect(keys(harmonics(phase_liouvillian, ω, t)))) == [-1, 0, 1]

  # A symbolic expression is assumed nonnegative as a whole; its factors are not sign-split.
  @test liouvillian(zero_H; channels=(jump(a, -γ),)) == -γ * dissipator(a)

  @test_throws ArgumentError jump(a, -1)
  @test_throws ArgumentError jump(a, -1.0)
  @test_throws ArgumentError jump(a, -2 + 0im)
  @test_throws ArgumentError jump(a, Num(-1))
  @test_throws ArgumentError jump(a, im)
  @test_throws ArgumentError jump(a, γ + im * δ)
  @test_throws ArgumentError jump(a, expim(ω * t₀))
end

@testset "collapse amplitudes and jump rates remain distinct representations" begin
  zero_H = zero(SQA.QAdd)
  collapse_form = liouvillian(zero_H; channels=(collapse(2a),))
  rate_form = liouvillian(zero_H; channels=(jump(a, 4),))
  frame = DissipativeFrame(a)

  collapse_d = @inferred kossakowski(collapse_form, frame)
  rate_d = @inferred kossakowski(rate_form, frame)

  @test collapse_d == rate_d
  @test collapse_d[1, 1] == 4
  @test iszero(hamiltonian(collapse_form, frame))
  @test iszero(hamiltonian(rate_form, frame))
end

@testset "Liouvillian arithmetic collects equal terms" begin
  L = hamiltonian_action(a' * a)

  @test iszero(L - L)
  @test zero(L) + L == L
  @test 2 * L == L + L
  @test L * 2 == L + L
  @test SQA.simplify(L - L) == zero(L)
end

@testset "periodic Liouvillian channels collect equal terms" begin
  operator = a + cos(ω * t) * a'
  single = harmonics(liouvillian(zero(SQA.QAdd); channels=(collapse(operator),)), ω, t)
  doubled = harmonics(
    liouvillian(zero(SQA.QAdd); channels=(collapse(operator), collapse(operator))), ω, t
  )

  @test doubled == 2 * single
  @test harmonics(doubled(t), ω, t) == doubled
end

@testset "zero operator factors produce zero maps" begin
  zero_operator = zero(SQA.QAdd)
  @test iszero(hamiltonian_action(zero_operator))
  @test iszero(dissipator(zero_operator))
  @test iszero(compose(hamiltonian_action(a), hamiltonian_action(zero_operator)))
end

@testset "Liouvillian composition is map composition" begin
  coherent = hamiltonian_action(a' * a)
  dissipative = dissipator(a)
  composed = compose(coherent, dissipative)

  @test composed isa Liouvillian
  @test !iszero(composed)
  @test compose(coherent, zero(dissipative)) == zero(coherent)
  @test compose(zero(coherent), dissipative) == zero(coherent)
  @test SQA.commutator(coherent, dissipative) ==
    compose(coherent, dissipative) - compose(dissipative, coherent)
  @test compose(coherent, dissipative) != compose(dissipative, coherent)
end

@testset "Liouvillian channel adapters" begin
  @test liouvillian(a; channels=(collapse(a),)) == hamiltonian_action(a) + dissipator(a)
  @test liouvillian(a; channels=(jump(a, γ),)) == hamiltonian_action(a) + γ * dissipator(a)
  @test liouvillian(zero(SQA.QAdd); channels=(collapse(2a),)) == dissipator(2a)
  @test_throws MethodError jump(a)
  @test_throws ArgumentError liouvillian(a; channels=(a,))
end

@testset "Liouvillians use the common van Vleck expansion" begin
  static = liouvillian(a' * a; channels=(jump(a, γ),))
  driven = liouvillian(a)
  generator = PeriodicGenerator(Dict(0 => static, 1 => driven, -1 => driven), ω)
  expansion = floquet_expansion(
    generator, VanVleck(; algorithm=HoriDeprit(; complete_positive=Val(false))), 2
  )

  @test expansion isa FloquetExpansion
  @test effective_generator(expansion) isa Liouvillian
  @test micromotion(expansion) isa PeriodicGenerator{Liouvillian}
  @test effective_component(expansion, 0) == static
end

@testset "periodic channel forms lower through the public Liouvillian seam" begin
  H = a' * a + cos(ω * t) * (a + a')
  collapse_operator = expim(-ω * t) * a
  jump_operator = expim(-ω * t) * a
  rate = 1 + cos(ω * t)

  native = liouvillian(H; channels=(collapse(collapse_operator), jump(jump_operator, rate)))
  periodic = harmonics(native, ω, t)

  @test periodic isa PeriodicGenerator{Liouvillian}
  @test @inferred(harmonics(native, ω, t)) isa PeriodicGenerator{Liouvillian}
  @test harmonics(periodic(t), ω, t) == periodic
  @test periodic[0] == liouvillian(a' * a; channels=(collapse(a), jump(a, 1)))
  expected_oscillatory = hamiltonian_action((1 // 2) * (a + a')) + (1 // 2) * dissipator(a)
  @test periodic[1] == expected_oscillatory
  @test periodic[-1] == expected_oscillatory

  direct = floquet_expansion(
    native, ω, t, VanVleck(; algorithm=HoriDeprit(; complete_positive=Val(false))), 1
  )
  keyword = floquet_expansion(
    H,
    ω,
    t,
    VanVleck(; algorithm=HoriDeprit(; complete_positive=Val(false))),
    1;
    channels=(collapse(collapse_operator), jump(jump_operator, rate)),
  )
  explicit = floquet_expansion(
    periodic, VanVleck(; algorithm=HoriDeprit(; complete_positive=Val(false))), 1
  )

  @test direct isa FloquetExpansion
  @test effective_generator(direct) == effective_generator(explicit)
  @test effective_generator(keyword) == effective_generator(explicit)

  complete_only = harmonics(
    liouvillian(zero(SQA.QAdd); channels=(collapse(collapse_operator),)), ω, t
  )
  weighted_only = harmonics(
    liouvillian(zero(SQA.QAdd); channels=(jump(jump_operator, 1),)), ω, t
  )
  @test complete_only == weighted_only
end

@testset "operator and rate harmonics convolve in periodic channels" begin
  collapse_operator = a + cos(ω * t) * a'
  jump_operator = a + expim(-ω * t) * a'
  rate = 1 + cos(ω * t)
  periodic = harmonics(
    liouvillian(
      zero(SQA.QAdd); channels=(collapse(collapse_operator), jump(jump_operator, rate))
    ),
    ω,
    t,
  )

  @test sort!(collect(keys(periodic))) == [-2, -1, 0, 1, 2]
  @test harmonics(periodic(t), ω, t) == periodic
end

@testset "general Liouvillian coefficients keep the shared Fourier convention" begin
  coefficient = expim(ω * t₀) * cos(ω * t)
  @test_throws ArgumentError jump(a, coefficient)

  periodic = harmonics(coefficient * dissipator(a), ω, t)

  @test sort!(collect(keys(periodic))) == [-1, 1]
  @test periodic[-1] == (1 // 2) * expim(ω * t₀) * dissipator(a)
  @test periodic[1] == (1 // 2) * expim(ω * t₀) * dissipator(a)
  @test harmonics(periodic(t), ω, t) == periodic
end

@testset "non-periodic Liouvillian phases are rejected" begin
  coefficient = expim(ω * t^2)
  @test_throws ArgumentError harmonics(coefficient * dissipator(a), ω, t)
end

@testset "a zero periodic Liouvillian keeps its public zero prototype" begin
  periodic = harmonics(zero(Liouvillian), ω, t)

  @test periodic isa PeriodicGenerator{Liouvillian}
  @test iszero(periodic)
  @test time_average(periodic) == zero(Liouvillian)
  @test iszero(
    floquet_expansion(
      zero(Liouvillian),
      ω,
      t,
      VanVleck(; algorithm=HoriDeprit(; complete_positive=Val(false))),
      1,
    ) |> effective_generator,
  )
end

@testset "periodic Liouvillian action agrees with an independent matrix model" begin
  finite_space = NLevelSpace(:finite_validation, 2)
  sigma = Transition(finite_space, :sigma, 1, 2)
  H = sigma + sigma'
  collapse_operator = expim(-ω * t) * sigma
  jump_operator = sigma + expim(ω * t) * sigma'
  rate = 2
  native = liouvillian(H; channels=(collapse(collapse_operator), jump(jump_operator, rate)))
  periodic = harmonics(native, ω, t)
  substitutions = Dict(ω => 3.0, t => 0.37)

  H_matrix = tomatrix(H, 2, substitutions)
  C_matrix = tomatrix(collapse_operator, 2, substitutions)
  J_matrix = tomatrix(jump_operator, 2, substitutions)
  rate_value = 2.0
  expected = independent_liouvillian_matrix(H_matrix, C_matrix, J_matrix, rate_value)
  native_matrix = superoperator_matrix(native, finite_space, 2, substitutions)
  periodic_matrix = superoperator_matrix(periodic(t), finite_space, 2, substitutions)

  @test native_matrix ≈ expected
  @test periodic_matrix ≈ expected
end
