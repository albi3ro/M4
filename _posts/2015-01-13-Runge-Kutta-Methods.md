---
title: Runge-Kutta Methods
layout: post
comments: True
category: Numerics
author: Christina C. Lee
image: "https://albi3ro.github.io/M4/Images/ODE/graphic.png"
tags: [ODE]
description:  An introduction to the Runge-Kutta class of Ordinary Differential Equation solvers, specifically covering Euler, Implicit Euler, RK2, and RK4.  
---

# Ordinary Differential Equation Solvers: Runge-Kutta Methods

### Christina Lee

So what's an <i>Ordinary Differential Equation</i>?

Differential Equation means we have some equation (or equations) that have derivatives in them.  

The <i>ordinary</i> part differentiates them from <i>partial</i> differential equations, the ones with curly $\partial$ derivatives.  Here, we only have one <b>independent</b> variable, let's call it $t$, and one or more <b>dependent</b>  variables, let's call them $x_1, x_2, ...$.  In partial differential equations, we can have more than one independent variable.

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

Though STEM students such as I have probably spent thousands of hours pouring of ways to analytically solve both ordinary and partial differential equations, unfortunately, the real world is rarely so kind as to provide us with an exactly solvable differential equation.  At least for interesting problems.  

We can sometimes approximate the real world as an exactly solvable situation, but for the situation we are interested in, we have to turn to numerics.  I'm not saying those thousand different analytic methods are for nothing.  We need an idea ahead of time of what the differential equation should be doing, to tell if it's misbehaving or not.  We can't just blindly plug and chug.  

Today will be about introducing four different methods based on Taylor expansion to a specific order, also known as Runge-Kutta Methods.  We can improve these methods with adaptive stepsize control, but that is a topic for another time, just like the other modern types of solvers such as Richardson extrapolation and predictor-corrector.  

Nonetheless, to work with ANY computational differential equation solver, you need to understand the fundamentals of routines like Euler and Runge-Kutta, their error propagation, and where they can go wrong. Otherwise, you might misinterpret the results of more advanced methods.

<b>WARNING:</b> If you need to solve a troublesome differential equation for a research problem, use a package, like [DifferentialEquations](https://github.com/JuliaDiffEq/DifferentialEquations.jl).  These packages have much better error handling and optimization.

Let's first add our plotting package and colors.

```julia
Pkg.update()
Pkg.add("Gadfly")
Pkg.add("Colors")
Pkg.add("Lazy")
using Gadfly
# Solarized Colors that I like working with
# http://ethanschoonover.com/solarized
using Colors
base03=parse(Colorant,"#002b36");
base02=parse(Colorant,"#073642");
base01=parse(Colorant,"#586e75");
base00=parse(Colorant,"#657b83");
base0=parse(Colorant,"#839496");
base1=parse(Colorant,"#839496");
base2=parse(Colorant,"#eee8d5");
base3=parse(Colorant,"#fdf6e3");

yellow=parse(Colorant,"#b58900");
orange=parse(Colorant,"#cb4b16");
red=parse(Colorant,"#dc322f");
magenta=parse(Colorant,"#d33682");
violet=parse(Colorant,"#6c71c4");
blue=parse(Colorant,"#268bd2");
cyan=parse(Colorant,"#3aa198");
green=parse(Colorant,"#859900");
```

We will be benchmarking our solvers on one of the simplest and most common ODE's,

$$
\frac{d}{d t}x=x \;\;\;\;\;\;\; x(t)=C e^t
$$

Though this only has one dependent variable, we want to structure our code so that we can accommodate a series of dependent variables, $y_1,y_2,...,y_n$, and their associated derivative functions.  Therefore, we create a function for each dependent variable, and then `push` it into an array declared as type `Function`.

```julia
function f1(t::Float64,x::Array{Float64,1})
    return x[1]
end
f=Function[]
push!(f,f1)
```

### Euler's Method
{% include image.html img="M4/Images/ODE/graphic.png" title="Stepping" caption="Approximating the slope each step." %}
First published in Euler's <i>Instutionum calculi integralis</i> in 1768, this method gets a lot of milage, and if you are to understand anything, this method is it.  

We march along with step size $h$, and at each new point, calculate the slope.  The slope gives us our new direction to travel for the next $h$.

We can determine the error from the Taylor expansion of the function

$$
x_{n+1}=x_n+h f(x_n,t) + \mathcal{O}(h^2).
$$

In case you haven't seen it before, the notation $\mathcal{O}(x)$ stands for "errors of the order x".
Summing over the entire interval, we accumuluate error according to

$$N\mathcal{O}(h^2)= \frac{x_f-x_0}{h}\mathcal{O}(h^2)=h $$,

making this a <b>first order</b> method.  Generally, if a technique is $n$th order in the Taylor expansion for one step, its $(n-1)$th order over the interval.

```julia
function Euler(f::Array{Function,1},t0::Float64,x::Array{Float64,1},h::Float64)
    d=length(f)
    xp=copy(x)
    for ii in 1:d
        xp[ii]+=h*f[ii](t0,x)
    end

    return t0+h,xp
end
```

## Implicit Method or Backward Euler


If $f(t,x)$ has a form that is invertible, we can form a specific expression for the next step.  For example, if we use our exponential,

$$
x_{n+1}=x_n+ h f(t_{n+1},x_{n+1})
$$

$$
x_{n+1}-h x_{n+1}=x_n
$$

$$
x_{n+1}=\frac{x_n}{1-h}
$$

This expression varies for each differential equation and only exists if the function is invertible.

```julia
function Implicit(f::Array{Function,1},t0::Float64,x::Array{Float64,1},h::Float64)
    return t0+h,x[1]/(1-h)
end
```

## 2nd Order Runge-Kutta

So in the Euler Method, we could just make more, tinier steps to achieve more precise results. Here, we make <i>bettter</i> steps.  Each step itself takes more work than a step in the first order methods, but we win by having to perform fewer steps.

This time, we are going to work with the Taylor expansion up to second order,

$$
x_{n+1}=x_n+h f(t_n,x_n) + \frac{h^2}{2} f^{\prime}(t_n,x_n)+ \mathcal{O} (h^3).
$$

Define

$$
k_1=f(t_n,x_n),
$$

so that we can write down the derivative of our $f$ function, and the second derivative (curvature), of our solution,

$$
f^{\prime}(t_n,x_n)=\frac{f(t_n+h/2,x_n+h k_1/2)-k_1}{h/2}+\mathcal{O}(h^2).
$$

Plugging this expression back into our Taylor expansion, we get a new expression for $x_{n+1}$

$$
x_{n+1}=x_n+hf(t_n+h/2,x_n+h k_1/2)+\mathcal{O}(h^3)
$$

We can also interpret this technique as using the slope at the center of the interval, instead of the start.

```julia
function RK2(f::Array{Function,1},t0::Float64,x::Array{Float64,1},h::Float64)
    d=length(f)

    xp=copy(x)
    xk1=copy(x)

    for ii in 1:d
        xk1[ii]+=f[ii](t0,x)*h/2
    end
    for ii in 1:d
        xp[ii]+=f[ii](t0+h/2,xk1)*h
    end

    return t0+h,xp
end
```

## 4th Order Runge-Kutta
Wait! Where's 3rd order? There exists a 3rd order method, but I only just heard about it while fact-checking for this post.  RK4 is your dependable, multi-purpose workhorse, so we are going to skip right to it.  

$$
k_1= f(t_n,x_n)
$$

$$
k_2= f(t_n+h/2,x_n+h k_1/2)
$$

$$
k_3 = f(t_n+h/2, x_n+h k_2/2)
$$

$$
k_4 = f(t_n+h,x_n+h k_3)
$$

$$
x_{n+1}=x_n+\frac{h}{6}\left(k_1+2 k_2+ 2k_3 +k_4 \right)
$$

I'm not going to prove here that the method is fourth order, but we will see numerically that it is.

<i>Note:</i> I premultiply the $h$ in my code to reduce the number of times I have to multiply $h$.

```julia
function RK4(f::Array{Function,1},t0::Float64,x::Array{Float64,1},h::Float64)
    d=length(f)

    hk1=zeros(x)
    hk2=zeros(x)
    hk3=zeros(x)
    hk4=zeros(x)

    for ii in 1:d
        hk1[ii]=h*f[ii](t0,x)
    end
    for ii in 1:d
        hk2[ii]=h*f[ii](t0+h/2,x+hk1/2)
    end
    for ii in 1:d
        hk3[ii]=h*f[ii](t0+h/2,x+hk2/2)
    end
    for ii in 1:d
        hk4[ii]=h*f[ii](t0+h,x+hk3)
    end

    return t0+h,x+(hk1+2*hk2+2*hk3+hk4)/6
end
```

This next function merely iterates over a certain number of steps for a given method.  

```julia
function Solver(f::Array{Function,1},Method::Function,t0::Float64,
        x0::Array{Float64,1},h::Float64,N::Int64)
    d=length(f)
    ts=zeros(Float64,N+1)
    xs=zeros(Float64,d,N+1)

    ts[1]=t0
    xs[:,1]=x0

    for i in 2:(N+1)
        ts[i],xs[:,i]=Method(f,ts[i-1],xs[:,i-1],h)
    end

    return ts,xs
end
```

```julia
N=1000
xf=10
t0=0.
x0=[1.]
dt=(xf-t0)/N

tEU,xEU=Solver(f,Euler,t0,x0,dt,N);
tIm,xIm=Solver(f,Implicit,t0,x0,dt,N);
tRK2,xRK2=Solver(f,RK2,t0,x0,dt,N);
tRK4,xRK4=Solver(f,RK4,t0,x0,dt,N);

xi=tEU
yi=exp(xi);

errEU=reshape(xEU[1,:],N+1)-yi
errIm=reshape(xIm[1,:],N+1)-yi
errRK2=reshape(xRK2[1,:],N+1)-yi;
errRK4=reshape(xRK4[1,:],N+1)-yi;
```

```julia
plot(x=tEU,y=xEU[1,:],Geom.point,
Theme(highlight_width=0pt,default_color=green,
default_point_size=3pt))
```
<iframe src="/M4/Images/ODE/exp.js.svg"  style="border:none; background: #ffffff"
width="600px" height="600px"></iframe>

```julia
lEU=layer(x=tEU,y=xEU[1,:],Geom.point,
Theme(highlight_width=0pt,default_color=green,
default_point_size=3pt))

lIm=layer(x=tIm,y=xIm[1,:],Geom.point,
Theme(highlight_width=0pt,default_color=yellow,
default_point_size=3pt))

lRK2=layer(x=tRK2,y=xRK2[1,:],Geom.point,
Theme(highlight_width=0pt,default_color=cyan,
default_point_size=2pt))

lRK4=layer(x=tRK4,y=xRK4[1,:],Geom.point,
Theme(highlight_width=0pt,default_color=violet,
default_point_size=4pt))

lp=layer(x->e^x,-.1,10,Geom.line,Theme(default_color=red))


plot(lp,lEU,lIm,lRK2,lRK4,
Guide.manual_color_key("Legend",["Euler","Implicit","RK2","RK4","Exact"],
[green,yellow,cyan,violet,red]),
Coord.cartesian(xmin=9.5,xmax=10.1))
```
<iframe src="/M4/Images/ODE/comp_exp.js.svg"  style="border:none; background: #ffffff"
width="600px" height="600px"></iframe>

```julia
lEU=layer(x=xi,y=errEU,Geom.point,
Theme(highlight_width=0pt,default_color=green,
default_point_size=1pt))

lIm=layer(x=xi,y=errIm,Geom.point,
Theme(highlight_width=0pt,default_color=yellow,
default_point_size=1pt))

lRK2=layer(x=xi,y=errRK2,Geom.point,
Theme(highlight_width=0pt,default_color=cyan,
default_point_size=1pt))

lRK4=layer(x=xi,y=errRK4,Geom.point,
Theme(highlight_width=0pt,default_color=violet,
default_point_size=1pt))

plot(lEU,lIm,lRK2,lRK4,Scale.y_asinh,
Guide.manual_color_key("Legend",["Euler","Implicit","RK2","RK4"],
[green,yellow,cyan,violet]))
```
<iframe src="/M4/Images/ODE/comp_err.js.svg"  style="border:none; background: #ffffff"
width="600px" height="600px"></iframe>

## Scaling of the Error

I talked above about the error scaling either as $h,h^2$, or $h^4$.  I won't just talk but here will numerically demonstrate the relationship as well.  For a variety of different step sizes, the below code calculates the final error for each method.  Then we will plot the error and see how it scales.

```julia
t0=0.
tf=1.
dx=tf-t0
x0=[1.]

dt=collect(.001:.0001:.01)

correctans=exp(tf)
errfEU=zeros(dt)
errfIm=zeros(dt)
errfRK2=zeros(dt)
errfRK4=zeros(dt)



for ii in 1:length(dt)
    N=round(Int,dx/dt[ii])
    dt[ii]=dx/N

    tEU,xEU=Solver(f,Euler,t0,x0,dt[ii],N);
    tIm,xIm=Solver(f,Implicit,t0,x0,dt[ii],N);
    tRK2,xRK2=Solver(f,RK2,t0,x0,dt[ii],N);
    tRK4,xRK4=Solver(f,RK4,t0,x0,dt[ii],N);

    errfEU[ii]=xEU[1,end]-correctans
    errfIm[ii]=xIm[1,end]-correctans
    errfRK2[ii]=xRK2[1,end]-correctans
    errfRK4[ii]=xRK4[1,end]-correctans
end
```

```julia
lEU=layer(x=dt,y=errfEU,Geom.point,
Theme(highlight_width=0pt,default_color=green,
default_point_size=1pt))

lIm=layer(x=dt,y=errfIm,Geom.point,
Theme(highlight_width=0pt,default_color=yellow,
default_point_size=1pt))

lRK2=layer(x=dt,y=errfRK2*10^5,Geom.point,
Theme(highlight_width=0pt,default_color=cyan,
default_point_size=1pt))

lRK4=layer(x=dt,y=errfRK4*10^10,Geom.point,
Theme(highlight_width=0pt,default_color=violet,
default_point_size=1pt))

tempEU(x)=(errfEU[end]*x/.01)
tempIm(x)=(errfIm[end]*x/.01)
#le=layer([tempEU,tempIm],0,.01,Geom.line,Theme(default_color=base01))
le=layer(tempEU,0,.01,Geom.line,Theme(default_color=base01))
lei=layer(tempIm,0,.01,Geom.line,Theme(default_color=base01))


temp2(x)=(errfRK2[end]*(x/.01)^2*10^5)
l2=layer(temp2,0,.01,Geom.line,Theme(default_color=base00))

temp4(x)=(errfRK4[end]*(x/.01)^4*10^10)
l4=layer(temp4,0,.01,Geom.line,Theme(default_color=base00))

xl=Guide.xlabel("h")
ylrk2=Guide.ylabel("Error e-5")
ylrk4=Guide.ylabel("Error e-10")
yl=Guide.ylabel("Error")

pEU=plot(lEU,lIm,le,lei,xl,yl,Guide.title("Euler and Implicit, linear error"))
p2=plot(lRK2,l2,xl,ylrk2,Guide.title("RK2, error h^2"))
p4=plot(lRK4,l4,xl,ylrk4,Guide.title("RK4, error h^4"))
vstack(hstack(p2,p4),pEU)
```

<iframe src="/M4/Images/ODE/err_scale.js.svg"  style="border:none; background: #ffffff"
width="600px" height="600px"></iframe>

## Arbitrary Order
While I have presented 4 concrete examples, many more exist. For any choice of variables $a_i, \beta_{i,j},a_i$ that fulfill

$$
x_{n+1}=x_n+h\left(\sum_{i=1}^s a_i k_i \right)+ \mathcal{O}(h^p)
$$

with

$$
k_i=f\left(t_n+\alpha_i h,x_n+h\left(\sum_{j=1}^s \beta_{ij} k_j \right) \right)
$$

we have a Runge-Kutta method of order $p$, where $p\geq s$.  The Butcher tableau provides a set of consistent coefficients.


These routines are also present in the M4 folder in a module named `diffeq.jl`. For later work, you may simply import the module.

Stay tuned for when we tuned these routines to the stiff van der Pol equations!
