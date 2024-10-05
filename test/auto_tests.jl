################################################################################
#                          Automatic test generation
################################################################################

@testsnippet AutoTests begin
    struct SupportItem{T, Dim, CRS, G <: Meshes.Geometry{Meshes.𝔼{Dim}, CRS}}
        name::String
        type::Type{T}
        geometry::G
        integral::Bool
        lineintegral::Bool
        surfaceintegral::Bool
        volumeintegral::Bool
        gausslegendre::Bool
        gausskronrod::Bool
        hadaptivecubature::Bool
    end

    # Constructor to explicitly convert Ints (0,1) to Bool values
    SupportItem(name, type, geometry, checkboxes::Vararg{I, 7}) where {I <: Integer} = SupportItem(
        name, type, geometry, Bool.(checkboxes)...)

    # If method is supported, test it on scalar- and vector-valued functions.
    # Otherwise, test that its use throws a MethodError
    function integraltest(intf, geometry, rule, supported, T)
        f(::Point) = T(1)
        fv(::Point) = fill(T(1), 2)

        if supported
            a1 = intf(f, geometry, rule)
            b1 = measure(geometry)
            @test a1 ≈ b1
            @test typeof(a1) == typeof(b1)
            @test intf(fv, geometry, rule) ≈ fill(b1, 2)
        else
            @test_throws "not supported" intf(f, geometry, rule)
        end
    end

    # Generate a @testset for item
    function autotest(item::SupportItem)
        N = (item.type == Float32) ? 1000 : 100
        algorithm_set = [
            (GaussLegendre(N), item.gausslegendre),
            (GaussKronrod(), item.gausskronrod),
            (HAdaptiveCubature(), item.hadaptivecubature)
        ]

        method_set = [
            (integral, item.integral),
            (lineintegral, item.lineintegral),
            (surfaceintegral, item.surfaceintegral),
            (volumeintegral, item.volumeintegral)
        ]

        itemsupport = Iterators.product(method_set, algorithm_set)

        # For each enabled solver type, run the test suite
        @testset "$(item.name)" begin
            for ((method, msupport), (alg, asupport)) in itemsupport
                integraltest(method, item.geometry, alg, msupport && asupport, item.type)
            end
        end
    end
end
