import Pkg; Pkg.activate(pwd()); Pkg.instantiate()
using DrWatson, Revise

@quickactivate :SAME

#include(srcdir("SAME.jl"))

using SAME


SAME.greet()
SAME.hello()
