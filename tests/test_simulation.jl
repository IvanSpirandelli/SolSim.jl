function run_all_simulation_tests()
    @testset verbose = true "Simulation" begin
        @testset verbose = true "Initialization" begin
            test_initialization()
        end  
    end
end

function test_initialization()
    @test isnothing(initialize_atoms("", [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], 1.0))
    @test !isnothing(initialize_atoms("", [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], 6.0))
end