function run_all_utilities_tests()
    @testset verbose = true "Utilities" begin
        @testset verbose = true "Vec to Ply" begin
            test_vec_to_ply()
        end
        @testset verbose = true "Ply to Vec" begin
            test_ply_to_vec()
        end
        @testset verbose = true "Poly to Vec" begin
            test_vec_to_poly()
        end
        @testset verbose = true "Vec to Poly" begin
            test_poly_to_vec()
        end
    end
end

function test_vec_to_ply()
    vec_to_ply([@SVector[0.0,0.0,0.0]], @SVector[0.5], "tests/test_files/single_created")

    f1 = open("tests/test_files/single.ply")
    f2 = open("tests/test_files/single_created.ply")

    @test readlines(f1) == readlines(f2)
end

function test_ply_to_vec()
    @test ply_to_vec("tests/test_files/single.ply") == ([@SVector[0.0,0.0,0.0]], @SVector[0.5])
end

function test_vec_to_poly()
    vec_to_poly([@SVector[0.0,0.0,0.0], @SVector[1.0,0.0,0.0], @SVector[0.0,1.0,0.0], @SVector[0.0,0.0,1.0]], @SVector[1.0, 1.0, 1.0, 1.0], "tests/test_files/cube_created")

    f1 = open("tests/test_files/cube.poly")
    f2 = open("tests/test_files/cube_created.poly")

    @test readlines(f1) == readlines(f2)
end

function test_poly_to_vec()
    @test poly_to_vec("tests/test_files/cube.poly") == ([@SVector[0.0,0.0,0.0], @SVector[1.0,0.0,0.0], @SVector[0.0,1.0,0.0], @SVector[0.0,0.0,1.0]], @SVector[1.0, 1.0, 1.0, 1.0])
end