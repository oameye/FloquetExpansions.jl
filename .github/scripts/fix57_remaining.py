from pathlib import Path

engine = Path('src/engine.jl')
text = engine.read_text()
text = text.replace(r'\sum_{n<N]', r'\sum_{n<N}', 1)
engine.write_text(text)

liouvillian = Path('src/liouvillian.jl')
text = liouvillian.read_text()
start = text.index('function _validated_jump_rate(rate::LiouvillianScalar)')
stop = text.index('\nend\n\n"""\n    jump(', start) + len('\nend')
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
end'''
liouvillian.write_text(text[:start] + new + text[stop:])
