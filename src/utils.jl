
"Returns true if expr is an exponential"
is_exp(expr) = SymbolicUtils.isterm(expr) && expr.f == exp
