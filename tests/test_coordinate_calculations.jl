function run_all_coordinate_calculation_tests()
    @testset verbose = true "Coordinate Calculations" begin
        @testset verbose = true "Intersections" begin 
            test_are_balls_intersecting()
        end
    end
end

function test_are_balls_intersecting()
    center_one = (SA_F64[1.5, 0, 0])
    center_two = (SA_F64[0, 0, 0])
    r_one = Float64(1)
    r_two = Float64(1)
    @test are_balls_intersecting(center_one, center_two, r_one, r_two)

    center_one = (SA_F64[1, 0, 0])
    center_two = (SA_F64[-0.2, 0, 0])
    r_one = Float64(0.5)
    r_two = Float64(0.5)
    @test !are_balls_intersecting(center_one, center_two, r_one, r_two)

    center_one = (SA_F64[0.5, -0.5, 0.5])
    center_two = (SA_F64[-0.5, 0.5, -0.5])
    r_one = Float64(0.5)
    r_two = Float64(0.5)
    @test !are_balls_intersecting(center_one, center_two, r_one, r_two)

    center_one = (SA_F64[0.5, 1, 1])
    center_two = (SA_F64[-0.5, 1, 1])
    r_one = Float64(0.5)
    r_two = Float64(0.5)
    @test !are_balls_intersecting(center_one, center_two, r_one, r_two)
end
