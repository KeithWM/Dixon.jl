using ReTest

using Dixon

include("TestDixon.jl")

retest(Dixon, TestDixon)