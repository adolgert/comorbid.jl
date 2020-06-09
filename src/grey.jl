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


function exact_burden_term(weights, prevalences, which)
    w = weights[which]
    yes = prevalences[which]
    no = prevalences[.!which]
    results = zeros(Float64, length(which))
    results[which] += w * ((1 - prod(1 .- w)) * prod(yes) * prod(1 .- no) / sum(w))
    results
end


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


function exact_burden_random(cnt)
    exact_burden(rand(cnt), rand(cnt))
end


function total_burden(weights, prevalences)
    1 - prod(1 .- weights .* prevalences)
end



function total_burden_random(cnt)
    w = rand(cnt)
    p = rand(cnt)
    exact = exact_burden(w, p)
    [total_burden(w, p), sum(exact)]
end
