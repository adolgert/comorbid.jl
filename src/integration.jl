using DifferentialEquations

function rescale_jacobian!(dbi_dpj_n, db_dpj, dbi_dpi)
    for col in 1:size(dbi_dpj_n, 2)
        true_sum_of_non_diagonals = db_dpj[col] - dbi_dpi[col]
        actual_sum_of_non_diagonals = 0.0
        for row in 1:size(dbi_dpj_n, 1)
            if row != col
                actual_sum_of_non_diagonals += dbi_dpj_n[row, col]
            end
        end
        scale = true_sum_of_non_diagonals / actual_sum_of_non_diagonals
        for row in 1:size(dbi_dpj_n, 1)
            if row == col
                dbi_dpj_n[row, col] = dbi_dpi[col]
            else
                dbi_dpj_n[row, col] *= scale
            end
        end
    end
end

function integrated_explanation()
    cnt = 5
    d = rand(cnt)
    p = rand(cnt)

    b = total_burden(d, p)
    b_i = exact_burden(d, p)

    # Approximate values of the Jacobian.
    dbi_dpj_n = order_n_partial(d, p, cnt)
    # Approximate column sum of the Jacobian.
    db_dpj_n = sum(dbi_dpj_n, dims = 1)

    # The true column sum of the Jacobian.
    db_dpj = similar(p)
    partial_wrt_prevalence!(db_dpj, d, p)

    # The true diagonal of the Jacobian.
    dbi_dpi = b_i ./ p
    dbi_dpj_n
    rescale_jacobian!(dbi_dpj_n, db_dpj, dbi_dpi)
end


"""
t from 0 to 1.
p(t) = prevalence * t
dp = prevalence * dt
"""
struct BurdenPath
    d::Vector{Float64}
    prevalence::Vector{Float64}
    order::Int
end


"""
    (bp::BurdenPath)(db, b, params, t)

This is the derivative along the path from prevalence=0
at t=0 to full prevalence at t=1.
The arguments are the output derivatives, `db`, the input
total burden, `b` at time t, and parameters, which are unused.
"""
function (bp::BurdenPath)(db, b, params, t)
    if t > 0
        # Calulate an approximate Jacobian and then dot it with
        # the line integral direction, from (0,0,0...) to (p_1, p_2, p_3...).
        p = bp.prevalence * t
        # Approximate values of the Jacobian for first `db.order` terms.
        dbi_dpj_n = order_n_partial(bp.d, p, bp.order)

        # The true column sum of the Jacobian is easy-enough to calculate.
        db_dpj = similar(p)
        partial_wrt_prevalence!(db_dpj, bp.d, p)

        # The true diagonal of the Jacobian is very simple.
        dbi_dpi = b ./ p

        # We use the true values to fix the approximate values.
        rescale_jacobian!(dbi_dpj_n, db_dpj, dbi_dpi)
        db .= dbi_dpj_n * bp.prevalence
    else
        db .= bp.d .* bp.prevalence
    end
end


"""
    integrated_burden(weights, prevalences, order)

Calculates the burden for given weights and prevalences using integration.
The total sum of burden, by cause, is difficult to calculate because
individuals with many comorbidities contribute significantly to the total
burden. This means many of the combinatorial terms have to be included.
This function works around that problem by calculating the exact burden at
low prevalence and calculating the change in burden as prevalences increase.
In this way, the cross-terms, coming from comorbidities, contribute gradually
to the total burden, even if we don't explicitly calculate comorbidities for
individuals with many health outcomes.
"""
function integrated_burden(weights, prevalences, order)
    bp = BurdenPath(weights, prevalences, order)
    b0 = zeros(length(weights))
    tspan = (0.0, 1.0)
    prob = ODEProblem(bp, b0, tspan)
    # Tsit5(), BS3()
    sol = solve(prob, Tsit5(), save_everystep=false)
    sol[end]
end
