# Mini Circuit Simulation Tutorial

Hi! Here's a little tutorial for those curious what SPICE does to model nonlinear parts of your circuits.

## Prereq
- Some calculus skills
- Some time
- Comfortable with matrices maybe

## Example
Let's take this circuit as our first example:

<img src="https://github.com/meowstr/circuitsim_tutorial/blob/main/schem.png" width="300">

The two diodes are what make this circuit worth modeling. They are nonlinear (exponential in fact), with most of their current range sitting around 0.6v. This of course is why we can usually assume them to be 0.6v drops instead of doing any math with exponentials since we often only need a rough approximation. But in the situation that we need a more exact answer, we can't get away with this. A fairly nice model for a diode is:

```I = I_s e^(V/V_T)```

Where 
```
I_s = 10e-13 A
V_T = 25 mV
```

## Mini-example
It's very simple at a glance. You'll find however that even just one diode and one resistor in series is not intuitive to solve for. Say theres a 1k resistor and diode connected in series, with 5v across them. What is the voltage across the diode? Let x be that voltage. So `(5 - x) / 1k = I_s e^(x/V_T)`. We just need to solve for x here, but it's difficult without using something like a lambert function (which you often won't find in your calculator or programming language). Simplest approach would be to solve for x numerically, using a root finding method like [bisection method](https://en.wikipedia.org/wiki/Bisection_method). If you wanna test yourself here, you should get 612.615 mV for that example.

## Solving the main example
Back to our original example. We have 3 voltages, v1, v2, and v3, that we want to find. At each of those points we know that KCL applies so we have 3 equations (labelling them F_1, F_2, F_3):

```
F_1 = i1 - i2 = 0
F_2 = i3 - i4 = 0
F_3 = i2 + i4 - i5 = 0

F = [ F_1  F_2  F_3 ]
```

where

```
i1 = (5 - v1) / 1k
i2 = (v1 - v3) / 1k
i3 = (9 - v2) / 1k
i4 = I_s e^((v2 - v3) / V_T)
i5 = I_s e^(v3 / V_T)
```

So we now have a system of equations that each equal 0. In other words, we want to find the root/zeroes of this big multivariate function called F. We can't just use bisection method this time because that only works for one variable (though there is something similar that extends to n-dimensions, I'm not sure I would recommend it). Instead we can use something called [Newton-Raphson method](https://en.wikipedia.org/wiki/Newton%27s_method]) which works for any number of variables. For good intuition, I think it helps to view the method with one variable first. Essentially it works by starting with a guess, then taking the slope on the function at that guess. Using that slope we construct a line, which acts as an approximation of the function at our guess. We can solve for the zero of our approximated function much easier since its just a line, of which will be our next guess. You do this over and over until you get the precision that you want, taking the current guess, making a line approximation of the function at the guess, finding the zero of that and making it the next guess.

Extending this to n-variables means extending the concept of slope and line approximation to generalized derivatives and linear approximation. A [Jacobian](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant) is a generalized derivative (or slope). Without describing too much theory (just read the wikipedia page for that), the Jacobian is a matrix of derivatives/slopes of our big multi-variable, multi-output function. because there are multiple variables, there are multiple slopes to think about. We want to know how the function changes as we change each variable, i.e. we want the rate of change with respect to each variable (also known as a partial derivative). The Jacobian in our case is a 3x3 matrix:

```
J(v1, v2, v3) = [ grad( F_1 ) ]
                [ grad( F_2 ) ]
                [ grad( F_3 ) ]
```
where
```
grad (f) = [ derivative of f w.r.t. v1,  derivative of f w.r.t. v2,  derivative of f w.r.t. v3 ]
```

The linear approximation of our big F function about v1, v2, v3 (which is our system of equations we're looking to solve) is:

```
F_approx(V1 V2 V3) = J(v1, v2, v3) * ( [V1 V2 V3] - [v1 v2 v3] ) + F(v1 v2 v3)
```

Solving for the zero of F_approx is simply a linear matrix solve:

```
F_approx = 0
J * [V1 V2 V3] + F(v1 v2 v3) = 0
J * [V1 V2 V3] = -F(v1 v2 v3)
```

There's many linear algebra libraries that can solve matrix equations, so it's up to you of course. But I used LAPACK's `dgesv` function for this. Was simple enough after struggling to understand how the leading dimension stuff worked in it. After 5 or 6 Newton steps, it's likely you can get the solution to our original equations, and therefore the voltages of our example circuit. However there is some nuance.

It's possible that using too small of initial guesses will cause the Newton steps to overshoot the solution, then be too far up one of our exponentials to ever return. A bandaid to this is to limit how much each step is allowed to move up or down the function. If the step causes a change in one the values of our function above a fixed threshold, then cut the step in half and test again. This worked for me, though I'm sure there's a more sophisticated way to decide if a change is significant.

Anyway, that's the gist of it I think. At least that should cover what is written in the code and provide you a good starting point for circuit simulation. SPICE uses Newton-Raphson much like this as far as I know. To check your work, my solution to the example circuit is

```
v1, v2, v3 = 2.816510 V  1.259826 V  0.633019 V
i1, i2, i3, i4, i5 = 2.183490 mA  2.183490 mA  7.740174 mA  7.740174 mA  9.923664 mA
```

Hope the rest of your day goes well, and thanks for reading!
