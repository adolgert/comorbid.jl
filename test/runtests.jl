using comorbid
using Test

@testset "comorbid.jl" begin
    @test isapprox(exact_burden_term([.1], [.8], [true]), [.08])
end
