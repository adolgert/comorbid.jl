"""
These are two ways to calculate the total burden of a population, attributed to
each cause, by short-circuiting the complete sum. The complete sum has 2^n terms,
where n is the number of causes. These examples try to cut off after burden
terms get small or the number of comorbidities gets high. Neither works super well.
"""

function nextburden(state, leaf = false)
    which = copy(state)
    cstate = cumsum(which)
    depth = cstate[end]
    if depth == 0
        which[1] = true
        depth = 1
    else
        idx = findfirst(cstate .== depth)
        if idx == length(which)
            if depth == 1
                fill!(which, false)
            else
                depth -= 1
                which[idx] = false
                up_one = findfirst(cstate .== depth)
                which[up_one] = false
                which[up_one + 1] = true
            end
        elseif leaf
            which[idx] = false
            which[idx + 1] = true
        else
            which[idx + 1] = true
        end
    end
    which
end


function show_iterate_burdens(n)
    state = zeros(Bool, n)
    for i in 1:2^n
        state = nextburden(state)
        println(state)
    end
end

"""
  treeburden(weights, prevalences)

Calculate the burden by searching term-by-term until terms get too small.
"""
function treeburden(weights, prevalences, epsilon)
    n = length(prevalences)
    b = zeros(Float64, n)
    which = nextburden(zeros(Bool, n))
    iter_count = 0
    while any(which)
        term = exact_burden_term(weights, prevalences, which)
        b += term
        if sum(term) >= epsilon
            which = nextburden(which)
        else
            which = nextburden(which, true)
        end
        iter_count += 1
    end
    iterfrac = iter_count / (2^n - 1)
    println("treeburden $(iterfrac)")
    b
end


function depthburden(weights, prevalences, maxdepth)
    n = length(prevalences)
    b = zeros(Float64, n)
    which = nextburden(zeros(Bool, n))
    iter_count = 0
    while any(which)
        term = exact_burden_term(weights, prevalences, which)
        b += term
        if sum(which) < maxdepth
            which = nextburden(which)
        else
            which = nextburden(which, true)
        end
        iter_count += 1
    end
    iterfrac = iter_count / (2^n - 1)
    println("treeburden $(iterfrac)")
    b
end


function tree_burden_random(cnt, epsilon)
    w = rand(cnt)
    p = rand(cnt)
    exact = exact_burden(w, p)
    tree = treeburden(w, p, epsilon)
    [sum(tree), total_burden(w, p), sum(exact), maximum(abs.(tree .- exact))]
end



"""
Using an arbitrary depth cutoff performs pretty poorly too. Even evaluating
50% of the total possible terms results in 30% error.
"""
function depth_burden_random(cnt, depth)
    w = rand(cnt)
    p = rand(cnt)
    exact = exact_burden(w, p)
    tree = depthburden(w, p, depth)
    total = sum(exact)
    trial = tree * total / sum(tree)
    [sum(tree), total_burden(w, p), sum(exact), maximum(abs.((trial .- exact) ./ exact))]
end
