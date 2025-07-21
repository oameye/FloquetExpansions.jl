
"Returns true if expr is an exponential"
is_exp(expr) = SymbolicUtils.isterm(expr) && expr.f == exp

isQMul(x) = x isa SQA.QMul

function trigonometric_to_exponential(x::QTerm)
    av = SQA.average(x)
    av_simplified = SymbolicUtils.simplify(SymbolicUtils.expand(QuestBase.trig_to_exp(av)))
    return SQA.undo_average(av_simplified)
end
