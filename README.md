# comorbid.jl

[![Build Status](https://travis-ci.com/adolgert/comorbid.jl.svg?branch=master)](https://travis-ci.com/adolgert/comorbid.jl)

Calculates the contribution of multiple independent diseases to comorbidity.

It's not possible to calculate the burden of disease on a population exactly, due
to the exponential size of the number of terms we need to calculate. This repository
looks at several ways to approximate that calculation. It isn't the code we use to do
the calculation. That code spends most of its time streaming population prevalences from
databases and streaming results back, and it works with health states, not causes as
such. A cause has health states; a health state has burden.
Those details don't matter for finding faster math, though.

When we assess the cost of disease for a population, we combine two factors, the
prevalence of each disease and the burden associated with each disease. The prevalence
is a fraction of people who have it, and we consider separate diseases as independently
likely for this calculation, even though we know that's not the case. The burden is
a measure of how much that disease affects a person. There is an international standard
for this measure, even though different countries may value health in different ways.
That measure is a number between zero and one for each cause of disease. For people
who are comorbid, meaning they have more than one cause of disease at the same time,
we need to add the burden numbers in such a way that they always add up to less than one.

If a single person has multiple causes of disease, we say their total burden
is `1 - product(1 - b_i)`, where `b_i` is the burden of each of the causes this person
has. That keeps the total under one, where one represents the burden of death.
We take that total and reassign it back to the causes proportionally. So the cause
from burden j is `(b_j / sum(b_i)) (1 - product(1 - b_i))`. That's the burden for a person
who has the set of causes given by the `i`.

To calculate the probability a person has the set of causes given by `i`, we compare it
to the classic problem of how many condiments you can put on a hamburger. You can have
ketchup or no ketchup, mustard or no mustard, pickles or no pickles. If there are `n`
condiments, there are `2^n` ways to make a burger. In this case, there are `n` causes
of disease, so there are `2^n` ways to be sick. The probability of any one comorbidity
combination is `prod(p_i)prod(1-p_j)` where `p_i` is the probability to be sick with
a cause and `p_j` is the probability to be sick with a different cause.

This repository finds ways to calculate the sum of the terms
`(b_j / sum(b_i)) (1 - product(1 - b_i))prod(p_i)prod(1-p_j)`
over all combinations of causes.

* Exact calculation. This takes about 40s for 25 causes. We need it for over 500 causes.

* Fake calculation. We have an exact expression for total comorbidity. This reallocates it.
  It has error above ten percent in some sample cases.

* Exponential version. This uses a modified version of the above equation that's close
  in value and can be calculated exactly. This has large error for some cases but seems
  maybe OK when there are more causes. It has potential.

* Sum with cutoffs. This looks at two different ways to short-circuit the exact
  calculation. One throws out cases with too many comorbidities. The other stops
  adding terms when the burden estimates are under a cutoff. Both perform poorly.
