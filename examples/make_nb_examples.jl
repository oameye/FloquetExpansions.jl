using Literate

using VanVleckRecursion

### Process examples
# Always rerun examples
const EXAMPLES_IN = @__DIR__
const OUTPUT_NB_DIR = @__DIR__

examples = filter!(file -> file[(end - 2):end] == ".jl", readdir(EXAMPLES_IN; join=true))
filter!(file -> !contains(file, "make_nb_examples"), examples)

for example in examples
    Literate.notebook(example, OUTPUT_NB_DIR; documenter=false, execute=true)
end
