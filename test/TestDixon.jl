module TestDixon
using ReTest

using Dixon

@testset "Test Taylor" begin
    testvector = [1, 2]
    T = eltype(testvector)
    @test Dixon.create_taylors(testvector, DixonElliptic, :cm) ==
        Dixon.Taylor1{T}([1, 0, 0, 2], 3)
    @test Dixon.create_taylors(testvector, DixonElliptic, :sm) ==
        Dixon.Taylor1{T}([0, 1, 0, 0, 2], 4)
    taylor = Dixon.create_taylors(testvector, DixonElliptic, :cm)
    @test taylor(2) == 17
end

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
    @testset "Coefficients of degree $n" for (n, expected) in table
        coeffficients = Dixon.Coefficients{Rational,Dixon.DixonElliptic}(n + 1)
        @test length(coeffficients.cm) == 3 * n + 1
        @test length(coeffficients.sm) == 3 * n + 2
        @test coeffficients.cm.coeffs[1:3:end] == [table[j][:c] for j in 0:n]
        @test coeffficients.sm.coeffs[2:3:end] == [table[j][:s] for j in 0:n]
    end
end

COEFFICIENTS = Dixon.Coefficients{Float64,Dixon.DixonElliptic}(100, 0, (1, 0))

cm(x::Number) = compute(x, COEFFICIENTS, Val{:cm}())
sm(x::Number) = compute(x, COEFFICIENTS, Val{:sm}())

@testset "Test specific values" begin
    table = Dict(
        # -1//3	=> (cm=Inf, sm=Inf, broken=true),
        -1//6 => (cm=cbrt(2), sm=-1, broken=false),
        0 => (cm=1, sm=0, broken=false),
        1//6 => (cm=1 / cbrt(2), sm=1 / cbrt(2), broken=false),
        1//3 => (cm=0, sm=1, broken=false),
        1//2 => (cm=-1, sm=cbrt(2), broken=false),
        # 2//3	=> (cm=Inf, sm=Inf, broken=true),
        -1//4 => (
            cm=(1 + sqrt(3) + sqrt(2 * sqrt(3))) / 2,
            sm=(-1 - sqrt(3 + 2 * sqrt(3))) / cbrt(4),
            broken=false,
        ),
        -1//12 => (
            cm=(-1 + sqrt(3) + sqrt(2 * sqrt(3))) / (2 * cbrt(2)),
            sm=(-1 + sqrt(3) - sqrt(2 * sqrt(3))) / (2 * cbrt(2)),
            broken=false,
        ),
        1//12 => (
            cm=(-1 + sqrt(3 + 2 * sqrt(3))) / cbrt(4),
            sm=(1 + sqrt(3) - sqrt(2 * sqrt(3))) / 2,
            broken=false,
        ),
        1//4 => (
            cm=(1 + sqrt(3) - sqrt(2 * sqrt(3))) / 2,
            sm=(-1 + sqrt(3 + 2 * sqrt(3))) / cbrt(4),
            broken=false,
        ),
        5//12 => (
            cm=(-1 + sqrt(3) - sqrt(2 * sqrt(3))) / (2 * cbrt(4)),
            sm=(-1 + sqrt(3) + 2 * sqrt(2 * sqrt(3))) / (2 * cbrt(4)),
            broken=true,
        ),
        7//12 => (
            cm=(-1 - sqrt(3 + 2 * sqrt(3))) / cbrt(4),
            sm=(1 + sqrt(3) + sqrt(2 * sqrt(3))) / 2,
            broken=false,
        ),
    )
    @testset "Value $z * pi3" for (z, expected) in table
        @test isapprox(cm(z * pi3), expected[:cm]; rtol=1.e-10) broken = expected[:broken]
        @test isapprox(sm(z * pi3), expected[:sm]; rtol=1.e-10) broken = expected[:broken]
    end
end

@testset "Test periodicity" begin
    samples = (-1.01:0.1:1) * 6  # avoid using pi3 and 0 to prevent degenerate cases.
    @testset "Testing periodicity for ω=$omega." for omega in
                                                     exp.((0:2) .* 2im * pi / 3) .* pi3
        for (x, y) in Iterators.product(samples, samples)
            z = Complex(x, y)
            @test isapprox(
                Dixon.center(z + omega, DixonElliptic),
                Dixon.center(z, DixonElliptic);
                atol=1.e-10,
            )
            @test isapprox(sm(z), sm(z + omega); rtol=1.e-10)
            @test isapprox(cm(z), cm(z + omega); rtol=1.e-10)
        end
    end
end

end  # module