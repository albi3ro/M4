# Ordinary Differential Equation Solvers
## Christina Lee

So what's an <i>Ordinary Differential Equation</i>?

Differential Equation means we have some equation (or equations) that have derivatives in them.  

The <i>ordinary</i> part differentiates them from <i>partial</i> differential equations (the ones with curly $\partial$ derivatives).  Here, we only have one <b>independent</b> variable (let's call it $t$), and one or more <b>dependent</b>  variables (let's call them $x_1, x_2, ...$).  In partial differential equations, we can have more than one independent variable.

This ODE can either be written as a system of the form
$$
\frac{d x_1}{dt}=f_1(t,x_1,x_2,...,x_n)
$$
$$
\frac{d x_2}{dt}=f_2(t,x_1,x_2,...,x_n)
$$
...
$$
\frac{d x_n}{dt}=f_n(t,x_1,x_n,...,x_n)
$$
or a single n'th order ODE of the form
$$
f_n(t,x) \frac{d^n x}{dt^n}+...+f_1(t,x) \frac{dx}{dt}+f_0(t,x)=0
$$
that can be rewritten in terms of a system of first order equations by performing variable substitutions such as
$$
\frac{d^i x}{dt^i}=x_i
$$

Though STEM students such as I have probably spent thousands of hours pouring of ways to analytically solve both ordinary and partial differential equations, unfortunately the real world is rarely so kind as to provide us with an exactly solvable differential equation.  At least for interesting problems.  

We can sometimes approximate the real world as an exactly solvable situation, but for the situation we are really interested in, we have to turn to numerics.  This isn't saying those thousand different analytic methods are for nothing.  We need an idea ahead of time of what the differential equation should be doing, in order to tell if it's misbehaving or not.  We can't just blindly plug and chug.  

Today will be about introducing four different methods based on Taylor expansion to a specific order, also known as Runge-Kutta Methods.  We can improve these methods with adaptive stepsize control, but that is a topic for another time.  Just like the other modern types of solvers, Richardson extrapolation and predictor-corrector.  

Nonetheless, in order to work with ANY computational differential equation solver, you need to understand the fundementals of routines like Euler and Runge Kutta, their error propogation, and where they can go wrong. Otherwise, you might misinterpret the results of more advanced routines. 

<b>WARNING:</b> If you actually need to solve a troublesome differential equation for a research problem, use a package, like [DifferentialEquations](https://github.com/JuliaDiffEq/DifferentialEquations.jl).  These packages have much better error handling and optimization.

Lets first add our plotting package and colors.