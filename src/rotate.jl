
function _rotate(input::QMul, a::QSym, ω, t; extra_term=true)::QTerm
    out = 1.0
    for arg in input.args_nc
        if arg isa SQA.Create
            out *= exp(-im * ω * t) * arg
        elseif arg isa SQA.Destroy
            out *= exp(im * ω * t) * arg
        else
            error("Unknown type: $(typeof(arg))")
        end
    end
    if extra_term
        input.arg_c * SymbolicUtils.simplify(out) - ω * a' * a
    else
        input.arg_c * SymbolicUtils.simplify(out)
    end
end
function _rotate(input::QAdd, a::QSym, ω, t; extra_term=true)
    if extra_term
        return sum(arg -> _rotate(arg, a, ω, t; extra_term=false), input.arguments) -
               ω * a' * a
    else
        return sum(arg -> _rotate(arg, a, ω, t; extra_term=false), input.arguments)
    end
end
function _rotate(in::Create, a::QSym, ω, t; extra_term=true)
    extra_term ? exp(-im * ω * t) * in - ω * a' * a : exp(-im * ω * t) * in
end
function _rotate(in::Destroy, a::QSym, ω, t; extra_term=true)
    extra_term ? exp(im * ω * t) * in - ω * a' * a : exp(im * ω * t) * in
end
function _rotate(input, a::QSym, ω, t; extra_term=true)
    extra_term ? input - ω * a' * a : input
end

function rotate(input::Union{QTerm,QSym}, a::QSym, ω, t)
    rot = _rotate(input, a, ω, t; extra_term=false)
    av_rot = SQA.average(rot)
    return SQA.undo_average(QuestBase.trig_to_exp(av_rot)) - ω * a' * a
end
