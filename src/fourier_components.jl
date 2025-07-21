# const Qexpression = Union{QTerm,QSym,Number,SymbolicUtils.BasicSymbolic}

function is_fourier(x, ω, t)
    isexp = is_exp(x)
    hasωt = occursin(ω, x) && occursin(t, x)
    # _iscomplex = iscomplex(first(SymbolicUtils.arguments(x)))
    return isexp && hasωt # && _iscomplex
end
function get_fourier_mul(mul::SymbolicUtils.BasicSymbolic, ω, t)
    # SymbolicUtils.ismul(mul) ? error("Not a Mul: $mul") : nothing
    args = SymbolicUtils.issym(mul) ? [mul] : SymbolicUtils.arguments(mul)
    idxs = findall(x -> is_fourier(x, ω, t), args)
    length(idxs) > 1 ? error("More than one Fourier component: $mul") : nothing
    rest = args[setdiff(1:end, idxs)]
    return args[idxs], prod(rest)
end

function add_components!(dict::Dict{Int,QAdd}, i::Int, components::QMul)
    if haskey(dict, i)
        dict[i] += components
    else
        dict[i] = QAdd([components])
    end
    return nothing
end

function remove_constants!(rot::QTerm)
    filter!(isQMul, SymbolicUtils.arguments(rot))
    return rot
end

function extract_factor_exponential(exp::SymbolicUtils.BasicSymbolic, ω, t)::Int
    @assert is_fourier(exp, ω, t) "Not a Fourier component: $exp"
    args = first(SymbolicUtils.arguments(exp))
    multiple = Int(SymbolicUtils.substitute(args, Dict(ω => -im, t => 1)))
    return multiple
end

function get_fourier_components(rot::QAdd, ω, t)
    components = Dict{Int,QAdd}()
    for term in SymbolicUtils.arguments(rot)
        if term isa QMul
            non_commutative = term.args_nc
            commutative = term.arg_c
            exponential_v, rest = get_fourier_mul(commutative, ω, t)
            if isempty(exponential_v)
                add_components!(components, 0, term)
            else
                exponential = first(exponential_v)
                factor = extract_factor_exponential(exponential, ω, t)
                qmul = QMul(rest, non_commutative)
                add_components!(components, factor, qmul)
            end
        end
    end
    return components
end
