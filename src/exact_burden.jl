using Plots
using ColorSchemes
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
    weights = rand(cnt)
    prevalences = rand(cnt)
    burden = exact_burden(weights, prevalences)
    weights, prevalences, burden
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


relerr(observed, expected) = (observed .- expected) ./ observed


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


function plot_compare(weight, prevalence, burden)
    total = sum(burden)
    estimated = (weight .* prevalence) * total / sum(weight .* prevalence)
    rel = relerr(estimated, burden)
    normalized = (rel .- minimum(rel)) / (maximum(rel) - minimum(rel))
    colors = get(ColorSchemes.buda, normalized)
    scatter(weight, prevalence, markercolor = colors)
end


"""
Error in the fake burden is above 10%, even for 20 causes.
"""
function fake_burden_basic_plot(cnt)
    plotter = scatter
    plot_cnt = 20
    for i in 1:plot_cnt
        w = rand(cnt)
        p = rand(cnt)
        exact = exact_burden(w, p)
        fake = fake_burden(w, p)
        rplot = plotter(
                fake,
                exact,
                markercolor = get(ColorSchemes.tab20b, (i-1)/(plot_cnt - 1)),
                legend = false
                )
        plotter = scatter!
        if i == plot_cnt
            display(rplot)
        end
    end
end


function generate_prevalences(cause_cnt, trial_cnt)
    weight = rand(cause_cnt)
    relerr_and_prevalence = zeros(Float64, cause_cnt, 2, trial_cnt)
    for t in 1:trial_cnt
        prevalence = rand(cause_cnt)
        burden = exact_burden(weight, prevalence)
        total = sum(burden)
        estimated = (weight .* prevalence) * total / sum(weight .* prevalence)
        rel = relerr(estimated, burden)
        relerr_and_prevalence[:, 1, t] = rel
        relerr_and_prevalence[:, 2, t] = prevalence
    end
    weight, relerr_and_prevalence
end


"""
This plot asks how the error of a single cause depends on its
prevalence. It plots relative error by prevalence. It shows that
larger prevalence is estimated below what it should be. There is scatter in
this plot. I'll bet that scatter is a function of the total prevalence for
all causes for each trial. Either as a sum or as it appears in the total burden
calculation.
"""
function plot_by_cause(weight, relerr_and_prevalence)
    cause_cnt = length(weight)
    trial_cnt = size(relerr_and_prevalence, 3)
    callit = scatter
    for cidx in 1:5:cause_cnt
        color = get(ColorSchemes.tab20b, (cidx - 1) / (cause_cnt - 1))
        rel = reshape(relerr_and_prevalence[cidx, 1, :], trial_cnt)
        prev = reshape(relerr_and_prevalence[cidx, 2, :], trial_cnt)
        display(callit(prev, rel, markercolor = color))
        callit = scatter!
    end
end


"""
This plots prevalence by weight. It colors by relative error. We see that
larger weights tend to skew towards positive relative error.
"""
function plot_trials(weight, relerr_and_prevalence)
    cause_cnt = length(weight)
    trial_cnt = size(relerr_and_prevalence, 3)
    w = repeat(weight, trial_cnt)
    rel = reshape(relerr_and_prevalence[:, 1, :], cause_cnt * trial_cnt)
    normrel = (rel .- minimum(rel)) / (maximum(rel) - minimum(rel))
    p = reshape(relerr_and_prevalence[:, 2, :], cause_cnt * trial_cnt)
    colors = get(ColorSchemes.viridis, normrel)

    @show length(w)
    @show length(rel)
    @show length(p)
    display(scatter(w, p, markercolor = colors))
end
