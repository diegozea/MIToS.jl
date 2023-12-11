"""
You need to `include` this file (and module) into your Julia session to run the tests using 
`ReTest`. `ReTest` is not a dependency of `MIToS`, so you need to install it manually.
Then, you can do `MIToSTests.retest()` to setup the tests. After that, you can use 
regular expressions to filter the tests you want to run, e.g. `MIToSTests.retest(r"MSA")`.
"""
module MIToSTests
    using ReTest
    include("tests.jl")
end