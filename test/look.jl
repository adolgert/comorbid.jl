
fake_burden_basic_plot((d, p) -> integrated_burden(d, p, 6), 20)
fake_burden_basic_plot(fake_gen, 20)
weight, relerr_and_prevalence = generate_prevalences((d, p) -> integrated_burden(d, p, 2), 10, 10)
plot_by_cause(weight, relerr_and_prevalence)
