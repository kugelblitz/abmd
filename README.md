# The ABMD Integrator

## About

This is an implementation of Adams–Bashforth–Moulton multistep integration
method with support for delay differential equations and advanced differential
equations.

The algorithm is described in the paper
[“On the Extension of Adams–Bashforth–Moulton Methods for Numerical Integration
of Delay Differential Equations and Application to the Moon's Orbit”](
https://arxiv.org/pdf/1903.02098)
by Dan Aksim and Dmitry Pavlov.


## Example

A usage example is given in [test/main.c](test/main.c).
The equation being integrated is a simplified 2D delay equation of
lunar motion in Earth's gravity field:
```math
\ddot{\mathbf{r}} = -\frac{(\mu + \mu_\oplus) \mathbf r}{r^3}
                    - \frac{3k_2\mu}{r^3} \left( 1 + \frac{\mu}{\mu_\oplus} \right)
                      \left( \frac{a_e}{r} \right)^5 \mathbf r(t-\tau),
```
which is inspired by eq. (5.213-1) from
[“Explanatory Supplement to the Astronomical Almanac”](
https://books.google.ru/books/about/?id=WBiqdNy_2KIC).

First, the equation is integrated from $`t_0=0`$ to $`t_1=365~\mathrm{d}`$.
Then, the final state at $`t_1`$ is taken as the new initial condition
and the equation is integrated backward in time from $`t_1`$ to $`t_0`$,
which is effectively the same as integrating an advanced equation
(with $`\mathbf r(t+\tau)`$) forward in time.
The difference in Moon positions between forward and backward runs
is plotted to `diff.png`.
