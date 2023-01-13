module TestDixon
using ReTest

using Dixon

cm(x::Number) = compute(x, Dixon.DixonRecursion, :cm)
sm(x::Number) = compute(x, Dixon.DixonRecursion, :sm)

@testset "Test pi3" begin
    @test isapprox(Dixon.pi3, 5.29991625; rtol=1.e-9)
end

@testset "Test coefficients" begin
    table = Dict(
        0 => (c=1, s=1),
        1 => (c=-1//3, s=-1//6),
        2 => (c=1//18, s=2//63),
        3 => (c=-23//2268, s=-13//2268),
        4 => (c=25//13608, s=23//22113),
        5 => (c=-619//1857492, s=-2803//14859936),
    )
    for (n, expected) in table
        coeffficients = Dixon.Coefficients{Rational, Dixon.DixonElliptic}(n+1)
        @test length(coeffficients.cm) == n + 1
        @test length(coeffficients.sm) == n + 1
        for i in 0:n
            @test coeffficients.cm == [table[j][:c] for j in 0:n]
            @test coeffficients.sm == [table[j][:s] for j in 0:n]
        end
    end
end

@testset "Test specific values" begin
    table = Dict(
        -1//3 * pi3	=> (cm=Inf, sm=Inf),
        -1//6 * pi3	=> (cm=cbrt(2), sm=-1),
        0	=> (cm=1, sm=0),
        1//6 * pi3	=> (cm=1/cbrt(2), sm=1/cbrt(2)),
        1//3 * pi3	=> (cm=0, sm=1),
        1//2 * pi3	=> (cm=-1, sm=1/cbrt(2)),
        2//3 * pi3	=> (cm=Inf, sm=Inf),
    )
    for (z, expected) in table
        @test isapprox(sm(z), expected[:sm]; rtol=1.e-10)
        @test isapprox(cm(z), expected[:cm]; rtol=1.e-10)
    end
end

@testset "Test periodicity" begin
    samples = (0:1:1) * pi3
    for omega in exp.((0:2) .* 2 * pi / 3) .* pi3
        for (x, y) in Iterators.product(samples, samples)
            z = Complex(x, y)
            @test isapprox(sm(z), sm(z + omega); rtol=1.e-10)
            @test isapprox(cm(z), cm(z + omega); rtol=1.e-10)
        end
    end
end

end  # module