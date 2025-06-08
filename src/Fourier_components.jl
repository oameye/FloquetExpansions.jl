const Qexpression = Union{QTerm,QSym,Number,SymbolicUtils.BasicSymbolic}

function is_fourier(x, ω, t)
    isexp = is_exp(x)
    hasωt = occursin(ω, x) && occursin(t, x)
    return isexp && hasωt
end
function get_fourier_mul(mul::SymbolicUtils.BasicSymbolic, ω, t)
    SymbolicUtils.ismul(mul) ? error("Not a Mul: $mul") : nothing
    args = SymbolicUtils.arguments(mul)
    idxs = findall(x -> is_fourier(x, ω, t), args)
    length(idxs) > 1 ? error("More than one Fourier component: $mul") : nothing
    rest = args[setdiff(1:end, idxs)]
    idxs, args[idxs], prod(rest) # Vector{Int}, Vector{BasicSymbolic}, BasicSymbolic
end

function add_components!(dict::Dict{Int,Qexpression}, i::Int, components)
    i ∈ keys(dict) ? dict[i] += components : dict[i] = components
    nothing
end
function add_components!(dict::Dict{Int,Qexpression}, term, ω, t)
    idxs, exponentials, rest = get_fourier_mul(term, ω, t)
    if isempty(idxs)
        add_components!(dict, 0, SQA.undo_average(rest))
    else
        args = first(SymbolicUtils.arguments(first(exponentials))) # can be dealt with better
        multiple = Int(SymbolicUtils.substitute(args, Dict(ω=>-im, t=>1)))
        add_components!(dict, multiple, SQA.undo_average(rest))
    end
end
function get_fourier_components(rot::SymbolicUtils.BasicSymbolic, ω, t)
    components = Dict{Int,Qexpression}()
    if SymbolicUtils.isadd(rot)
        for term in SymbolicUtils.arguments(rot)
            add_components!(components, term, ω, t)
        end
    elseif SymbolicUtils.ismul(rot)
        add_components!(components, rot, ω, t)
    else
        error("Not a sum or product: $rot")
    end
    return components
end

function remove_constants(rot::SymbolicUtils.BasicSymbolic)
    avgsym = Symbolics.filterchildren(w -> w isa SQA.Average, rot)
    constant = SymbolicUtils.substitute(rot, Dict(avgsym .=> 0))
    rot - constant
end
