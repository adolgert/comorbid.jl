
function exponential_burden_term(weights, prevalences, which)
    w = weights[which]
    yes = prevalences[which]
    no = prevalences[.!which]
    results = zeros(Float64, length(which))
    results[which] += w * ((1 - exp(-sum(w))) * prod(yes) * prod(1 .- no) / sum(w))
    results
end


function exponential_burden_term!(weights, prevalences, which, running_sum)
    w = weights[which]
    yes = prevalences[which]
    no = prevalences[.!which]
    running_sum[which] += w * ((1 - exp(-sum(w))) * prod(yes) * prod(1 .- no) / sum(w))
end


"""
    exponential_burden(weights, prevalences)

Calculate the burden using an alternative formula that has exponentials in it.
"""
function exponential_burden(weights, prevalences)
    n = length(prevalences)
    which = zeros(Bool, n)
    b = zeros(Float64, n)
    for i in 1:(1<<n - 1)
        to_flip = grey_next_flip(i)
        which[to_flip] = !which[to_flip]
        exponential_burden_term!(weights, prevalences, which, b)
    end
    b
end


"""
As the number of causes gets larger, the error of this method decreases.
It's a lot for 5 causes but around 1% for 20 causes. That's the maximum
relative error for any one burden.
"""
function exponential_burden_random(cnt)
    w = rand(cnt)
    p = rand(cnt)
    exact = exact_burden(w, p)
    exp_burden = exponential_burden(w, p)
    total = total_burden(w, p)
    trial = exp_burden * total / sum(exp_burden)
    [sum(exp_burden), total, sum(exact), maximum(abs.((exp_burden - exact) ./ exact))]
end
