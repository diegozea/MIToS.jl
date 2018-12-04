using Pkg

if !Sys.iswindows()
    if !haskey(Pkg.installed(), "GaussDCA")
        Pkg.add(PackageSpec(url="https://github.com/carlobaldassi/GaussDCA.jl"))
    end

    msa = map(Residue, rand(1:21,100,20))
    dca = gaussdca(msa, min_separation=2)

    @test  isnan(dca[1,1])
    @test  isnan(dca[1,2])
    @test !isnan(dca[1,3])
end
# TODO: Solve error in Windows
#
# ERROR: LoadError: syntax: invalid escape sequence
# Stacktrace:
#  [1] include at .\boot.jl:317 [inlined]
#  [2] include_relative(::Module, ::String) at .\loading.jl:1038
#  [3] include(::Module, ::String) at .\sysimg.jl:29
#  [4] exec_options(::Base.JLOptions) at .\client.jl:239
#  [5] _start() at .\client.jl:432
# in expression starting at C:\Users\appveyor\AppData\Local\Temp\1\jl_4DB2.tmp.jl:3
# Information: Error During Test at C:\projects\mitos-jl\test\runtests.jl:46
#   Got exception LoadError("C:\\projects\\mitos-jl\\test\\Information\\Externals.jl", 8, ErrorException("failed process: Process(`'C:\\julia\\bin\\julia.exe' 'C:\\Users\\appveyor\\AppData\\Local\\Temp\\1\\jl_4DB2.tmp.jl'`, ProcessExited(1)) [1]")) outside of a @test
#   LoadError: failed process: Process(`'C:\julia\bin\julia.exe' 'C:\Users\appveyor\AppData\Local\Temp\1\jl_4DB2.tmp.jl'`, ProcessExited(1)) [1]
#   Stacktrace:
#    [1] error(::String, ::Base.Process, ::String, ::Int64, ::String) at .\error.jl:42
#    [2] pipeline_error at .\process.jl:712 [inlined]
#    [3] #run#509(::Bool, ::Function, ::Cmd) at .\process.jl:670
#    [4] run at .\process.jl:668 [inlined]
#    [5] #gaussdca#39(::String, ::Base.Iterators.Pairs{Symbol,Int32,Tuple{Symbol},NamedTuple{(:min_separation,),Tuple{Int32}}}, ::Function, ::Array{Residue,2}) at C:\projects\mitos-jl\src\Information\Externals.jl:35
