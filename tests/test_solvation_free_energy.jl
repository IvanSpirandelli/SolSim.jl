function run_all_solvation_free_energy_tests()
    @testset verbose = true "Solvation Free Energy" begin
        @testset verbose = true "Integration" begin
            test_solvation_free_energy()
        end

        @testset verbose = true "Area" begin 
            test_area()
        end

        @testset verbose = true "Volume" begin
            test_volume()
        end

        @testset verbose = true "Mean Curvature" begin 
            test_mean()
        end

        @testset verbose = true "Gaussian Curvature" begin 
            test_gaussian()
        end
    end
end

function test_solvation_free_energy()
    c1 = @SVector[0.0, 0.0, 0.0]
    c2 = @SVector[2.0, 0.0, 0.0]
    r1 = 1.0
    r2 = 1.0

    centers = [c1, c2]
    solute_radii = [r1, r2]

    prefactors = @SVector[3.0/8.0,1.0/8.0,1.0/8.0,1.0/8.0]
    e = solvation_free_energy(centers, solute_radii, 0.0, prefactors)
    @test isapprox(e,4*pi)

    rs = 0.475
    η = 0.475
    radii = [r1, r2]
    centers = [@SVector[0.0, 0.0, 0.0], @SVector[2.0, 0.0, 0.0]]
    prefactors = get_wb_prefactors(rs, η)

    e1 = solvation_free_energy(centers, radii, rs, prefactors)
    e2 = solvation_free_energy(centers, radii, rs, η)
    @test isapprox(e1, e2)

    e4 = solvation_free_energy(nothing, Atoms{2}(centers, solute_radii), rs, η)
    @test isapprox(e2, e4)

    e5 = solvation_free_energy(1, Atoms{2}(centers, solute_radii), rs, η)
    @test isapprox(e4, e5)

    e6 = solvation_free_energy(2, Atoms{2}(centers, solute_radii), rs, η)
    @test isapprox(e5, e6)

    e7 = solvation_free_energy_for_hard_spheres(nothing, centers, solute_radii, rs, η)
    @test isapprox(e6, e7)

    e8 = solvation_free_energy_for_hard_spheres(1, centers, solute_radii, rs, η)
    @test isapprox(e7, e8)

    e9 = solvation_free_energy_for_hard_spheres(2, centers, solute_radii, rs, η)
    @test isapprox(e8, e9)
end

function test_volume()
    v = [(SA_F64[0, 0, 0]), (SA_F64[2, 0, 0])]
    r = [1.0,1.0]
    measures = get_measures(v,r,0.0)
    @test measures[1] == 2.0 * (4.0/3.0) * pi
end

function test_area()
    v = [(SA_F64[0, 0, 0]), (SA_F64[2, 0, 0])]
    r = [1.0,1.0]
    measures = get_measures(v,r, 0.0)
    @test measures[2] == 2.0 * 4.0 * pi
end

function test_mean()
    v = [(SA_F64[0, 0, 0]), (SA_F64[2, 0, 0])]
    r = [1.0,1.0]
    measures = get_measures(v,r, 0.0)
    @test measures[3] == 2.0 * 4.0 * pi
end

function test_gaussian()
    v = [(SA_F64[0, 0, 0]), (SA_F64[2, 0, 0])]
    r = [1.0,1.0]
    measures = get_measures(v,r, 0.0)
    @test measures[4] == 2.0 * 4.0 * pi
end