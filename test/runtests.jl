using TensorsLite,TensorsLiteGeometry
using Test

a = Vec()
b = 2.0𝐢
c = 2.0𝐢 + 2.0𝐣
d = 2.0𝐣

@testset "Area calculation" begin
    @test area(a,b,c) == 2.0
    @test area(a,b,c,d) == 4.0
end

@testset "in_triangle function" begin
    @test in_triangle(a,a,b,c)[1]
    @test in_triangle(b,a,b,c)[1]
    @test in_triangle(c,a,b,c)[1]
    @test in_triangle((a+b+c)/3,a,b,c)[1]
    @test !(in_triangle(3𝐢,a,b,c)[1])
end