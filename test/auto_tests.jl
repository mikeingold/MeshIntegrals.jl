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
        #@assert item.type == coordtype(item.geometry) "Item type mismatch"

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

@testitem "Integrals" setup=[Setup, AutoTests] begin
    # Spatial descriptors
    origin3d(T) = Point(T(0), T(0), T(0))
    origin2d(T) = Point(T(0), T(0))
    ẑ(T) = Vec(T(0), T(0), T(1))
    plane_xy(T) = Plane(origin3d(T), ẑ(T))

    # Points on xy-plane at unit distance on axes
    pt_n(T) = Point(T(0), T(1), T(0))
    pt_w(T) = Point(T(-1), T(0), T(0))
    pt_e(T) = Point(T(1), T(0), T(0))
    pt_s(T) = Point(T(0), T(-1), T(0))
    pt_z(T) = Point(T(0), T(0), T(1))

    # Test Geometries
    ball2d(T) = Ball(origin2d(T), T(2.0))
    ball3d(T) = Ball(origin3d(T), T(2.0))
    box1d(T) = Box(Point(T(-1)), Point(T(1)))
    box2d(T) = Box(Point(T(-1), T(-1)), Point(T(1), T(1)))
    box3d(T) = Box(Point(T(-1), T(-1), T(-1)), Point(T(1), T(1), T(-1)))
    circle(T) = Circle(plane_xy(T), T(2.5))
    cyl(T) = Cylinder(pt_e(T), pt_w(T), T(2.5))
    cylsurf(T) = CylinderSurface(pt_e(T), pt_w(T), T(2.5))
    disk(T) = Disk(plane_xy(T), T(2.5))
    parab(T) = ParaboloidSurface(origin3d(T), T(2.5), T(4.15))
    sphere2d(T) = Sphere(origin2d(T), T(2.5))
    sphere3d(T) = Sphere(origin3d(T), T(2.5))
    tetra(T) = Tetrahedron(pt_n(T), pt_w(T), pt_e(T), pt_n(T) + ẑ(T))
    triangle(T) = Ngon(pt_e(T), pt_n(T), pt_w(T))
    torus(T) = Torus(origin3d(T), ẑ(T), T(3.5), T(1.25))

    SUPPORT_MATRIX(T) = [
        # Name, T type, example,    integral,line,surface,volume,    GaussLegendre,GaussKronrod,HAdaptiveCubature
        SupportItem("Ball{2,$T}", T, ball2d(T), 1, 0, 1, 0, 1, 1, 1),
        SupportItem("Ball{3,$T}", T, ball3d(T), 1, 0, 0, 1, 1, 0, 1),
        SupportItem("Box{1,$T}", T, box1d(T), 1, 1, 0, 0, 1, 1, 1),
        SupportItem("Box{2,$T}", T, box2d(T), 1, 0, 1, 0, 1, 1, 1),
        SupportItem("Box{3,$T}", T, box3d(T), 1, 0, 0, 1, 1, 0, 1),
        SupportItem("Circle{$T}", T, circle(T), 1, 1, 0, 0, 1, 1, 1),
        SupportItem("Cylinder{$T}", T, cyl(T), 1, 0, 0, 1, 1, 0, 1),
        SupportItem("CylinderSurface{$T}", T, cylsurf(T), 1, 0, 1, 0, 1, 1, 1),
        SupportItem("Disk{$T}", T, disk(T), 1, 0, 1, 0, 1, 1, 1),
        # Frustum -- not yet supported
        SupportItem("ParaboloidSurface{$T}", T, parab(T), 1, 0, 1, 0, 1, 1, 1),
        # SimpleMesh
        SupportItem("Sphere{2,$T}", T, sphere2d(T), 1, 1, 0, 0, 1, 1, 1),
        SupportItem("Sphere{3,$T}", T, sphere3d(T), 1, 0, 1, 0, 1, 1, 1),
        SupportItem("Tetrahedron", T, tetra(T), 1, 0, 0, 1, 0, 1, 0),
        SupportItem("Triangle{$T}", T, triangle(T), 1, 0, 1, 0, 1, 1, 1),
        SupportItem("Torus{$T}", T, torus(T), 1, 0, 1, 0, 1, 1, 1)
    ]

    @testset "Float64 Geometries" begin
        map(autotest, SUPPORT_MATRIX(Float64))
    end
end
