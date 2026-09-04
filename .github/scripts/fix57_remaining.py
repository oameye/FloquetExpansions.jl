from pathlib import Path

engine = Path('src/engine.jl')
text = engine.read_text()
text = text.replace('\\sum_{n<N]','\\sum_{n<N}',1)
old = '''function micromotion(expansion::FloquetExpansion)
  result = zero(expansion.generator)
  for (order, kick) in enumerate(expansion.kick_components)
    result = result + reattach(kick, order)
  end
  return result
end

function micromotion(expansion::FloquetExpansion, n::Int)
  1 <= n < expansion.order ||
    throw(ArgumentError("order $(n) is outside 1:$(expansion.order - 1)"))
  return SQA.simplify(reattach(expansion.kick_components[n], n))
end
'''
new = '''function micromotion(
  expansion::FloquetExpansion{G,P,E,C,R}
) where {G,P<:PeriodicGenerator,E,C,R}
  result = zero(expansion.generator)::P
  for (order, kick) in enumerate(expansion.kick_components)
    result = result + reattach(kick, order)
  end
  return result::P
end

function micromotion(
  expansion::FloquetExpansion{G,P,E,C,R}, n::Int
) where {G,P<:PeriodicGenerator,E,C,R}
  1 <= n < expansion.order ||
    throw(ArgumentError("order $(n) is outside 1:$(expansion.order - 1)"))
  return SQA.simplify(reattach(expansion.kick_components[n], n))::P
end
'''
if old not in text:
    raise SystemExit('micromotion block not found')
engine.write_text(text.replace(old,new,1))

liouvillian = Path('src/liouvillian.jl')
text = liouvillian.read_text()
old = '''function _validated_jump_rate(rate::LiouvillianScalar)
  if rate isa Real && !(rate isa Symbolics.Num) && rate < 0
    throw(ArgumentError("jump rate must be nonnegative; got `$rate`"))
  end

  coefficient = convert(SQA.CNum, rate)
  value = SQA.to_num(coefficient)
  imaginary_part = imag(value)
  iszero(imaginary_part) ||
    throw(ArgumentError("jump rate must be provably real; got `$rate`"))

  if value isa Real && !(value isa Symbolics.Num) && value < 0
    throw(ArgumentError("jump rate must be nonnegative; got `$rate`"))
  end
  return coefficient
end
'''
new = '''function _validated_jump_rate(rate::LiouvillianScalar)
  if rate isa Symbolics.Num
    return convert(SQA.CNum, rate)
  elseif rate isa Real
    rate < 0 && throw(ArgumentError("jump rate must be nonnegative; got `$rate`"))
    return convert(SQA.CNum, rate)
  elseif rate isa Complex
    iszero(imag(rate)) ||
      throw(ArgumentError("jump rate must be provably real; got `$rate`"))
    real_rate = real(rate)
    if real_rate isa Real && !(real_rate isa Symbolics.Num) && real_rate < 0
      throw(ArgumentError("jump rate must be nonnegative; got `$rate`"))
    end
    return convert(SQA.CNum, rate)
  end

  coefficient = convert(SQA.CNum, rate)
  value = SQA.to_num(coefficient)
  if value isa Symbolics.Num
    return coefficient
  elseif value isa Real
    value < 0 && throw(ArgumentError("jump rate must be nonnegative; got `$rate`"))
    return coefficient
  elseif value isa Complex
    iszero(imag(value)) ||
      throw(ArgumentError("jump rate must be provably real; got `$rate`"))
    real_value = real(value)
    if real_value isa Real && !(real_value isa Symbolics.Num) && real_value < 0
      throw(ArgumentError("jump rate must be nonnegative; got `$rate`"))
    end
    return coefficient
  end

  throw(ArgumentError("jump rate must be provably real; got `$rate`"))
end
'''
if old not in text:
    raise SystemExit('jump validator block not found')
liouvillian.write_text(text.replace(old,new,1))
