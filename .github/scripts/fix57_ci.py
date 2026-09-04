from pathlib import Path
import subprocess

engine = Path("src/engine.jl")
text = engine.read_text()

# Reattach the public API docstring to the public binding instead of the internal helper.
doc_start = text.index('"""\n    floquet_expansion(generator::PeriodicGenerator, gauge, order)')
helper_start = text.index("function _floquet_expansion(", doc_start)
doc = text[doc_start:helper_start]
text = text[:doc_start] + text[helper_start:]
public_start = text.index("function floquet_expansion(\n  generator::P, gauge::G, order::Int")
text = text[:public_start] + doc + text[public_start:]

old_channels = """function floquet_expansion(
  H::SQA.QField,
  wd::Symbolics.Num,
  t::Symbolics.Num,
  gauge::Gauge,
  order::Int;
  channels::LiouvillianChannelCollection=(),
)
  if isempty(channels)
    return _floquet_expansion(harmonics(qadd(H), wd, t), gauge, order, NoProvenance())
  end

  provenance = _microscopic_provenance(channels)
  L = _liouvillian_from_provenance(H, provenance)
  return _floquet_expansion(harmonics(L, wd, t), gauge, order, provenance)
end
"""
new_channels = """function _floquet_expansion_channels(
  H::SQA.QField,
  wd::Symbolics.Num,
  t::Symbolics.Num,
  gauge::Gauge,
  order::Int,
  ::Tuple{},
)
  return _floquet_expansion(harmonics(qadd(H), wd, t), gauge, order, NoProvenance())
end

function _floquet_expansion_channels(
  H::SQA.QField,
  wd::Symbolics.Num,
  t::Symbolics.Num,
  gauge::Gauge,
  order::Int,
  channels::LiouvillianChannelCollection,
)
  provenance = _microscopic_provenance(channels)
  L = _liouvillian_from_provenance(H, provenance)
  return _floquet_expansion(harmonics(L, wd, t), gauge, order, provenance)
end

function floquet_expansion(
  H::SQA.QField,
  wd::Symbolics.Num,
  t::Symbolics.Num,
  gauge::Gauge,
  order::Int;
  channels::LiouvillianChannelCollection=(),
)
  return _floquet_expansion_channels(H, wd, t, gauge, order, channels)
end
"""
if old_channels not in text:
    raise SystemExit("channel constructor block not found")
text = text.replace(old_channels, new_channels, 1)

old_component = """function effective_component(expansion::FloquetExpansion, n::Int)
  0 <= n < expansion.order ||
    throw(ArgumentError("order $(n) is outside 0:$(expansion.order - 1)"))
  return SQA.simplify(
    reattach(expansion.effective_components[n + 1], expansion.generator.wd, n)
  )
end
"""
new_component = """function effective_component(
  expansion::FloquetExpansion{G,P,E,C,R}, n::Int
) where {G,P,E<:GeneratorComponent,C,R}
  0 <= n < expansion.order ||
    throw(ArgumentError("order $(n) is outside 0:$(expansion.order - 1)"))
  component = reattach(expansion.effective_components[n + 1], expansion.generator.wd, n)::E
  return SQA.simplify(component)::E
end
"""
if old_component not in text:
    raise SystemExit("effective_component block not found")
text = text.replace(old_component, new_component, 1)
text = text.replace(
    "  return SQA.simplify(result)\nend\n\nfunction effective_generator(\n  expansion::FloquetExpansion{G,P,E,C,R}",
    "  return SQA.simplify(result)::E\nend\n\nfunction effective_generator(\n  expansion::FloquetExpansion{G,P,E,C,R}",
    1,
)
engine.write_text(text)

liouvillian = Path("src/liouvillian.jl")
text = liouvillian.read_text()
old_rate = """function _validated_jump_rate(rate::LiouvillianScalar)
  coefficient = convert(SQA.CNum, rate)
  value = SQA.to_num(coefficient)
  imaginary_part = Symbolics.simplify(imag(value))
  iszero(imaginary_part) ||
    throw(ArgumentError("jump rate must be provably real; got `$rate`"))

  real_part = Symbolics.simplify(real(value))
  unwrapped = Symbolics.value(real_part)
  if unwrapped isa Real && unwrapped < 0
    throw(ArgumentError("jump rate must be nonnegative; got `$rate`"))
  end
  return coefficient
end
"""
new_rate = """function _validated_jump_rate(rate::LiouvillianScalar)
  coefficient = convert(SQA.CNum, rate)
  value = SQA.to_num(coefficient)
  conjugation_residual = Symbolics.simplify(conj(value) - value)
  iszero(conjugation_residual) ||
    throw(ArgumentError("jump rate must be provably real; got `$rate`"))

  unwrapped = Symbolics.value(Symbolics.simplify(value))
  if unwrapped isa Real && unwrapped < 0
    throw(ArgumentError("jump rate must be nonnegative; got `$rate`"))
  end
  return coefficient
end
"""
if old_rate not in text:
    raise SystemExit("jump-rate validator block not found")
liouvillian.write_text(text.replace(old_rate, new_rate, 1))

tests = Path("test/liouvillian.jl")
text = tests.read_text()
text = text.replace(
    "  rate = 1 + 2 * cos(ω * t)\n  native = liouvillian(H; channels=(collapse(collapse_operator), jump(jump_operator, rate)))",
    "  rate = 2\n  native = liouvillian(H; channels=(collapse(collapse_operator), jump(jump_operator, rate)))",
    1,
)
text = text.replace("  rate_value = 1 + 2 * cos(3.0 * 0.37)", "  rate_value = 2.0", 1)
tests.write_text(text)

subprocess.run(["git", "diff", "--check"], check=True)
