"""
This file calculates exact burden and fake burden.
"""

"""
A Grey code is a map from a binary number to a different binary number.
The magic is that if you call this `n=2^i` times, you will see every binary
number between `0` and `n-1`, but each one will be only one digit off from
the last. This function, in particular, tells you which bit will flip next.
"""
function grey_next_flip(n)
    for i in 1:1:64
        if n & 1 == 1
            return i
        else
            n >>= 1
        end
    end
    65
end


"""
Calculates (b_j /sum(b_i))(1 - prod(1-b_i)) (prod(p_i) prod(1 - p_j)).
"""
function exact_burden_term(weights, prevalences, which)
    w = weights[which]
    yes = prevalences[which]
    no = prevalences[.!which]
    results = zeros(Float64, length(which))
    results[which] += w * ((1 - prod(1 .- w)) * prod(yes) * prod(1 .- no) / sum(w))
    results
end


"""
Calculates (b_j /sum(b_i))(1 - prod(1-b_i)) (prod(p_i) prod(1 - p_j))
and adds it to the given running sum.
"""
function exact_burden_term!(weights, prevalences, which, running_sum)
    w = weights[which]
    yes = prevalences[which]
    no = prevalences[.!which]
    running_sum[which] += w * ((1 - prod(1 .- w)) * prod(yes) * prod(1 .- no) / sum(w))
end


function exact_burden(weights, prevalences)
    n = length(prevalences)
    which = zeros(Bool, n)
    b = zeros(Float64, n)
    for i in 1:(1<<n - 1)
        to_flip = grey_next_flip(i)
        which[to_flip] = !which[to_flip]
        exact_burden_term!(weights, prevalences, which, b)
    end
    b
end


"""
Make random prevalences and burdens to calculate.
"""
function exact_burden_random(cnt)
    exact_burden(rand(cnt), rand(cnt))
end


"""
This is the super simple, exact total burden for a population. Given how simple it is,
it's surprising that we don't have a closed form for individual contributions to burden.
"""
function total_burden(weights, prevalences)
    1 - prod(1 .- weights .* prevalences)
end


"""
Maybe we can divide the total burden by some fraction in order to get
individual cause's burden. Nah.
"""
function fake_burden(weights, prevalences)
    total = total_burden(weights, prevalences)
    w = weights .* prevalences
    w * total / sum(w)
end


function total_burden_random(cnt)
    w = rand(cnt)
    p = rand(cnt)
    exact = exact_burden(w, p)
    [total_burden(w, p), sum(exact)]
end


"""
Error in the fake burden is above 10%, even for 20 causes.
"""
function fake_burden_random(cnt)
    w = rand(cnt)
    p = rand(cnt)
    exact = exact_burden(w, p)
    fake = fake_burden(w, p)
    [sum(fake), sum(exact), maximum(abs.((fake .- exact) ./ exact))]
end
