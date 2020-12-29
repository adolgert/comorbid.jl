using DataFrames, GLM

function construct_variables(cause_cnt, data_cnt)
    exact = zeros(Float64, cause_cnt * data_cnt)
    weight = zeros(Float64, cause_cnt * data_cnt)
    prevalence = zeros(Float64, cause_cnt * data_cnt)
    total = zeros(Float64, cause_cnt * data_cnt)
    bc = zeros(Float64, cause_cnt * data_cnt)
    bc2 = zeros(Float64, cause_cnt * data_cnt)
    w = rand(cause_cnt)
    for i in 1:data_cnt
        p = rand(cause_cnt)
        e = exact_burden(w, p)
        assign = ((i - 1) * cause_cnt + 1):(i * cause_cnt)
        exact[assign] = e
        weight[assign] = w
        prevalence[assign] = p
        bc[assign] = (w .* p) / sum(w .* p)
        bc2[assign] = (w .* p) / sum(e)
        total[assign] .= sum(e)
    end
    weight, prevalence, exact, total, bc, bc2
end


w, p, e, t, bc, bc2 = construct_variables(20, 20)
data = DataFrame(E = e, W = w, P = p, BC = bc, BC2 = bc2, T = t)
ols = lm(@formula(E ~ BC + W) , data)
scatter(e, predict(ols))


w1 = 0.3 * rand(10)
p1 = 0.6 * rand(10)
b1 = exact_burden(w1, p1)

w2 = 0.5 * rand(10)
p2 = 0.2 * rand(10)
b2 = exact_burden(w2, p2)

b12 = vcat(b1, b2)
b3 = exact_burden(vcat(w1, w2), vcat(p1, p2))
sum(b3)
1 - (1 - sum(b1)) * (1 - sum(b2))

# This is an oddly successful way to combine the two smaller exact
# calculations into a larger calculation. It's off by a constant.
b1p = b1 * (1 - sum(b2))
b2p = b2 * (1 - sum(b1))
scatter(b3[1:10], b1p, markercolor = :green)
scatter!(b3[11:20], b2p, markercolor = :blue)
xlabel!("Exact Burden")
ylabel!("Estimated Burden")

guess = vcat(b1p, b2p)
guess = guess * sum(b3) / sum(guess)
rel = (guess - b3) ./ b3
# Max relative error is 0.0363.
maximum(abs.(rel))
