# For integral(f, geometry, settings, FP) where FP is not Float64, ensure results
# have expected level of accuracy and are produce results in appropriate type

# Base value for atol when integrating with a particular FP type
@testsnippet BaseAtol begin
    baseatol = Dict(
        Float32 => 0.01f0,
        BigFloat => BigFloat(0.001)
    )
end

@testitem "Alternate floating types" setup=[Setup, BaseAtol] begin
    @testset "$FP" for FP in (Float32, BigFloat)
        # typeof @test's are currently broken for Float32, see GitHub Issue 74

        # Rectangular volume with unit integrand
        f = p -> one(FP)
        box1d = Box(Point(fill(zero(FP) * u"m", 1)...), Point(fill(one(FP) * u"m", 1)...))
        box2d = Box(Point(fill(zero(FP) * u"m", 2)...), Point(fill(one(FP) * u"m", 2)...))
        box3d = Box(Point(fill(zero(FP) * u"m", 3)...), Point(fill(one(FP) * u"m", 3)...))

        # Check HCubature integrals (same method invoked for all dimensions)
        int_HC = integral(f, box1d, HAdaptiveCubature(); FP = FP)
        @test int_HC≈one(FP) * u"m" atol=baseatol[FP] * u"m"
        @test typeof(int_HC.val)==FP broken=(FP == Float32)

        # Check Gauss-Legendre integral in 1D
        int_GL_1D = integral(f, box1d, GaussLegendre(100); FP = FP)
        @test int_GL_1D≈one(FP) * u"m" atol=baseatol[FP] * u"m"
        @test typeof(int_GL_1D.val)==FP broken=(FP == Float32)

        # Check Gauss-Legendre integral in 2D
        int_GL_2D = integral(f, box2d, GaussLegendre(100); FP = FP)
        @test int_GL_2D≈one(FP) * u"m^2" atol=2baseatol[FP] * u"m^2"
        @test typeof(int_GL_2D.val)==FP broken=(FP == Float32)

        # Check Gauss-Legendre integral in 3D
        int_GL_3D = integral(f, box3d, GaussLegendre(100); FP = FP)
        @test int_GL_3D≈one(FP) * u"m^3" atol=3baseatol[FP] * u"m^3"
        @test typeof(int_GL_3D.val)==FP broken=(FP == Float32)
    end
end

@testitem "Integral Aliases" setup=[Setup] begin
    f = p -> one(Float32)
    box1d = Box(Point(fill(0.0f0u"m", 1)...), Point(fill(1.0f0u"m", 1)...))
    box2d = Box(Point(fill(0.0f0u"m", 2)...), Point(fill(1.0f0u"m", 2)...))
    box3d = Box(Point(fill(0.0f0u"m", 3)...), Point(fill(1.0f0u"m", 3)...))
    box4d = Box(Point(fill(0.0f0u"m", 4)...), Point(fill(1.0f0u"m", 4)...))

    # Check alias functions for accuracy
    glrule = GaussLegendre(100)
    @test lineintegral(f, box1d, glrule, FP=Float32)≈1.0f0u"m" atol=0.01f0u"m"
    @test surfaceintegral(f, box2d, glrule, FP=Float32)≈1.0f0u"m^2" atol=0.02f0u"m^2"
    @test volumeintegral(f, box3d, glrule, FP=Float32)≈1.0f0u"m^3" atol=0.03f0u"m^3"

    # Check for unsupported use of alias functions
    harule = HAdaptiveCubature()
    @test_throws "not supported" lineintegral(f, box4d, harule, FP=Float32)
    @test_throws "not supported" surfaceintegral(f, box4d, harule, FP=Float32)
    @test_throws "not supported" volumeintegral(f, box4d, harule, FP=Float32)
end
