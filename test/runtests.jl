using TensorsLite,TensorsLiteGeometry
using Test

displacement = (1.3𝐢-1.3𝐣)
a = Vec() + displacement
b = 2.0𝐢 + displacement
c = 2.0𝐢 + 2.0𝐣 + displacement
d = 2.0𝐣 + displacement

@testset "Area calculation" begin
    @test area(a,b,c) ≈ 2.0
    @test area(a,b,c,d) ≈ 4.0
    @test area([a,b,c],[1,2,3]) ≈ 2.0
    @test area([a,b,c,d],[1,2,3,4]) ≈ 4.0
end

@testset "in_triangle function" begin
    @test in_triangle(a,a,b,c)[1]
    @test in_triangle(b,a,b,c)[1]
    @test in_triangle(c,a,b,c)[1]
    @test in_triangle((a+b+c)/3,a,b,c)[1]
    @test !(in_triangle(3𝐢+displacement,a,b,c)[1])
end

@testset "in_polygon function" begin
    @test in_polygon(a,[a,b,c,d],[1,2,3,4])
    @test in_polygon(b,[a,b,c,d],[1,2,3,4])
    @test in_polygon(c,[a,b,c,d],[1,2,3,4])
    @test in_polygon(d,[a,b,c,d],[1,2,3,4])
    @test !in_polygon(Vec(),[a,b,c,d],[1,2,3,4])
end