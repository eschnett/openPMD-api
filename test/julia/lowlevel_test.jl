import Base

# pass as first argument the path to libopenPMD.jl.so

module openPMD
using CxxWrap
@wrapmodule(ARGS[1])

function __init__()
    @initcxx
end
end

s = openPMD.CXX_Series(
    "../samples/lowlevel_julia_test.json", openPMD.ACCESS_CREATE)

function Base.getindex(
    cont::Cont,
    index,
) where
{
    A,
    B,
    Cont<:openPMD.CXX_Container{A,B}
}
    return openPMD.cxx_getindex(cont, index)
end

iteration = openPMD.cxx_iterations(s)[100]
mesh = openPMD.cxx_meshes(iteration)
# empty brackets to dereference CxxRef
# automatic dereferencing works only for Int types since CxxRwap.jl defines
# Base.getindex(x::CxxBaseRef, i::Int) = Base.getindex(x[], i)
mesh[]["E"][]["x"]
openPMD.cxx_flush(s, "{}")
