function run_all_perturbation_tests()
    @testset verbose = true "Perturbations" begin
        @testset verbose = true "Get Shift Vector" begin
            test_shifting()
        end
    end
end

function test_shifting()
    c1 = @SVector[0.52, 0.0, 0.0]
    c2 = @SVector[3.5, 0.0, 0.0]
    centers = [c1, c2]
    interaction_range = 2.0
    bounds = 4.0
    t = get_shift_vector(centers, interaction_range, bounds)
    @test t == [2.01, 0.0, 0.0]

    shift_centers!(centers, t, bounds)
    @test isapprox(centers[1][1], 2.51) 
    @test isapprox(centers[2][1], 1.49) 

    centers = [@SVector[0.1, 0.1, 0.1], @SVector[3.0, 0.1, 0.1]]
    interaction_range = 2.5
    bounds = 5.0
    t = get_shift_vector(centers, interaction_range, bounds)
    shift_centers!(centers, t, bounds)
    solute = Atoms{2}(
        centers,
        @SVector[1.0, 1.0]
    )
    @test !is_intersecting_set_of_atoms(solute)

    centers = [@SVector[0.1, 0.1, 0.1], @SVector[3.438748577248782, 0.11035852625725086, 0.17321492140895642]]
    interaction_range = 2.5
    bounds = 5.0
    t = get_shift_vector(centers, interaction_range, bounds)
    shift_centers!(centers, t, bounds)
    solute = Atoms{2}(
        centers,
        @SVector[1.0, 1.0]
    )
    @test is_intersecting_set_of_atoms(solute)

    centers = [
            @SVector[4.964062664063326, 4.924513050139993, 4.920997023332272],
            @SVector[3.0, 0.1, 0.1]
        ]
    interaction_range = 2.5
    bounds = 5.0
    t = get_shift_vector(centers, interaction_range, bounds)
    shift_centers!(centers, t, bounds)
    solute = Atoms{2}(
        centers,
        @SVector[1.0, 1.0]
    )
    @test is_intersecting_set_of_atoms(solute)
end