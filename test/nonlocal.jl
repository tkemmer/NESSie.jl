# check eye!
m = zeros(Int, 3, 3)
eye!(m);    @test m == [1 0 0; 0 1 0; 0 0 1]
eye!(m, 2); @test m == [2 0 0; 0 2 0; 0 0 2]
m = zeros(Int, 2, 3)
eye!(m);    @test m == [1 0 0; 0 1 0]
eye!(m, 2); @test m == [2 0 0; 0 2 0]
m = zeros(Int, 3, 2)
eye!(m);    @test m == [1 0; 0 1; 0 0]
eye!(m, 2); @test m == [2 0; 0 2; 0 0]
