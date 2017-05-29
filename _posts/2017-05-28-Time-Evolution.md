---
title: Time Evolution Split Operator Method
layout: post
comments: True
category: Prerequisites
tags: [Quantum]
author: Christina C. Lee
description: Quantum mechanical time evolution through the split operator method looking at Rabi Oscillations between two harmonic oscillator wells.
image: "https://albi3ro.github.io/M4/Images/TimeEvolution/shifting.png"
---

When cleaning my apartment, sometimes I just grab the nearest dirty thing to me and try to do something to it.  But that is not the most efficient way to get things done.  If I'm planning, I'll first dedicate my attention to one problem, like putting clothing away, then rotate my attention to something else, like dirty dishes.  I can keep focused on just one task and do it well.  Each problem I solve optimally in shorter intervals instead of tackling everything at once.  

That same principle applies to solving partial differential equations.  [<u>Splitting Methods in Communication, Imaging, Science, and Engineering</u>](#citation1) called this principle one of the big ideas of numerical computation.  In numerics, we call it <b>Strang splitting</b>.  

We will be applying Strang splitting to solve the Schrondinger equation, but people use the same idea to a variety of problems, like ones with different timescales, length scales, or physical processes.  We will be using it to separate out terms diagonal in position space from terms diagonal in momentum space.  

We can describe a class of general partial differential equations by
<img src="/M4/Images/TimeEvolution/eq5.png" class="eqim">

Over a small time step, the following approximation holds
<img src="/M4/Images/TimeEvolution/eq6.png" class="eqim2">

For <b>Strang splitting</b>, instead of applying both operators together, we break them up into two.  I'll discuss non-commutivity later.
<img src="/M4/Images/TimeEvolution/eq7.png" class="eqim2">
$U_1$ and $U_2$ are evolution operators. We can define

$$
\tilde{y}(0)= U_2 y(0)
$$

so that

$$
y (\delta t) = U_1 \tilde{y}(0)
$$

### Applying to Quantum Mechanics
Now let's take a look at the Schrodinger Equation:
<img src="/M4/Images/TimeEvolution/eq8.png" class="eqim">
The Hamiltonian gets seperated into position terms and momentum terms.  For ease, let's define our unitary evolution operators,
<img src="/M4/Images/TimeEvolution/eq9.png" class="eqim2">

I mentioned earlier that I would discuss non-communitvity.  We need to do that now.  We can't simply seperate the evolution operator for the full Hamiltonian into two parts, because we would introduce terms proportional to the commutator.  

<img src="/M4/Images/TimeEvolution/eq10.png" class="eqim2">

$e^{A+B}$ expanded has terms that look like $AB$ <b>and</b> $BA$, whereas $e^{A}e^{B}$ only has terms that look like $AB$.  We lose the symmetry of the expression.  We can gain back an order of accuracy by symmetrizing our formula, calculating a time step by

<img src="/M4/Images/TimeEvolution/eq11.png" class="eqim2">

But the next step will the start with $U_x (\delta t/2)$ !

<img src="/M4/Images/TimeEvolution/eq12.png" class="eqim">
<img src="/M4/Images/TimeEvolution/eq13.png" class="eqim">

All we need to do to add an order of accuracy is start the simulation with $U_x(\delta t/2)$ and end it with $U_x(\delta t/2)$, leaving everything else the same.  Pretty remarkable you can get that much of an improvement for that little work.  Once we apply this to a bunch of time steps, we get

<img src="/M4/Images/TimeEvolution/eq14.png" class="eqim">

We have to apply a few operators before starting the loop. Between the loop and a measurement, we have to apply an additional operator.

In the spatial domain, the momentum operator involves derivatives and is rather icky.  But in the momentum domain, we only have to multiply by $k^2/2m$.  Thanks to some nicely optimized libraries, we can just transform into the momentum domain with `fft`, solve the momentum problem there, and transform back with `ifft`.

## Rabi Oscillations

To demonstrate time evolution in a simple system with interesting physics, I chose to apply the split operator to Rabi Oscillations between two harmonic oscillators.  

To get an idea of what will happen, we will use a qualitative model of two states weakly coupled to each other by a parameter $\epsilon$.  If we have the two minima sufficifiently seperated from each other, tunneling will happen slowly and will not significantly affect the shape of the eigenfunctions and their energies $E_0$.  Instead of of solving for the shape of the wavefunction, we solve a two-state Hamiltonian that looks like this,

<img src="/M4/Images/TimeEvolution/eq1.png" class="eqim">

The eigenvalues and corresponding eigenvectors of the matrix are,

<img src="/M4/Images/TimeEvolution/eq15.png" class="eqim2">

<img src="/M4/Images/TimeEvolution/eq2.png" class="eqim">


If a wavefunction starts purely in the right state, we want to choose a combination of our eigenvectors that sets the left state to zero at $t=0$.  The resulting wavefunction will evolve as,

<img src="/M4/Images/TimeEvolution/eq3.png" class="eqim">

<img src="/M4/Images/TimeEvolution/eq4.png" class="eqim">


Thus this simple phenomological model shows us how we can expect the wavefunction to move back and forth between the two wells in a cyclical manner, with a period proportional to the rate of tunneling.

#### Packages
For the first time, I'm using `Plots.jl`, with `PlotlyJS` as the backend for `Plots`.  I had to coax my computer a bit to get the packages to work, but the errors were specific to my computer.  I choose to switch as I believe the `PlotlyJS` output will provide a better expierence for those viewing the GitHub-Pages static site, though it might cause more trouble for anyone using the jupyter notebooks.  If you are having trouble plotting yourself, I recommend just switching back to whatever package is easiest for you.

I used `Plots` to generate a gif for me directly, but I found my old method of generating a file of .png files and then using `ffmpeg` from the command line much fast and easier.


```julia
using Plots
plotlyjs()
#using Rsvg
#using PyPlot # if you easiest for you
```

## Input Parameters


```julia
# Set Time Parameters
t0=0
tf=40000
dt=.1

# Set Space Grid Parameters
dx=.1
xmax=8
# xmin will be -xmax.  Making the situation symmetric

# How far seperated are the potential minima
seperation=6;
# minima at seperation/2 and -seperation/2

# How often we measure occupation and view the state
nmeasure=1000;
```

## Automatic Evaluation Parameters

Given the parameters above, we can calculate the following variables that we will use in the code.

Note: `k` Gave me a bit of a headache. The algorithm depends quite a bit on the conventions `fft` decides to use ;(
Currently, I'm using odd `N`.  You'll have to change the formula if you use even `N`.


```julia
t=collect(t0:dt:tf)
x=collect(-xmax:dx:xmax)

nt=length(t)
N=length(x)

k = [ collect(0:((N-1)/2)) ; collect(-(N-1)/2:-1) ] *2*π/(N*dx);

occupation=zeros(Complex{Float64},floor(Int,nt/nmeasure),2);
```

## The Potentials and Evolution Operators


```julia
Vx=.5*(abs(x)-seperation/2).^2;
Vk=k.^2/2

Uxh=exp(-im*Vx*dt/2);
Ux=exp(-im*Vx*dt);
Uf=exp(-im*Vk*dt);
"potentials and evolvers defined"
```


```julia
plot(x, Vx)

plot!(xlabel="x", ylabel="V",
plot_title="Double Well Potential")
```



<iframe src="/M4/Images/TimeEvolution/potential.html"  style="border:none; background: #ffffff"
width="700px" height="450px"></iframe>

## The Unperturbed Wavefunctions

The ground state for a harmonic oscillator is a Gaussian
\begin{equation}
\langle x | \phi \rangle= \phi (x) = \frac{1}{\pi^{1/4}} e^{-\frac{x^2}{2}}
\end{equation}
We assume $\omega = \hbar = m = 1$ for sake of convenience.



```julia
ϕ(x)=π^(-.25)*exp(-x.^2/2)
```


```julia
ϕl=ϕ(x+seperation/2);
ϕr=ϕ(x-seperation/2);
Ψ0=convert(Array{Complex{Float64},1},ϕl);
```


```julia
plot(x,ϕl,label="ϕl")
plot!(x,ϕr,label="ϕr")

plot!(xlabel="x", ylabel="ϕ",
plot_title="Left and Right Wavefunctions")
```
<iframe src="/M4/Images/TimeEvolution/lrpsi.html"  style="border:none; background: #ffffff"
width="700px" height="450px"></iframe>

## FFT's

This algorithm runs a large number of Fast Fourier Transforms and Inverse Fast Fourier Transforms.  To speed the process, we can tell the computer to spend some time, in the beginning, allocating the right amount of space and optimizing the routine for the particular size and type of array we are passing it.  

The next cell does this, by using `plan_fft` and `plan_ifft` to generate objects that can act on our arrays as operators.


```julia
ft=plan_fft(Ψ0);
Ψf=ft*Ψ0;
ift=plan_ifft(Ψf);
```

## Occupancy of each state

To measure the occupancy of the total wavefunction in either the left or right well groundstate, I compute the value
\begin{equation}
c_r=\langle \Psi | \phi_r \rangle  = \int \Psi^* (x) \phi_r(x) dx \;\;\;\;\;\;\; p_r=c_r^2
\end{equation}
\begin{equation}
c_l = \langle \Psi | \phi_l \rangle = \int \Psi^* (x) \phi_l (x) dx\;\;\;\;\;\;\; p_l = c_l^2
\end{equation}
The probability of being in the state is the coefficient squared.

Though in theory these values will always be real, numerical errors introduce some errors, and Julia will assume that the answer is complex.  Therefore, we need to apply `abs` to make the numbers `Float64` instead of `Complex{Float64}`.


```julia
nmeas=1000
c=zeros(Float64,floor(Int,nt/nmeas),2);
```


```julia
# Uncomment the # lines to generate a gif.  Note: It takes a long time

Ψ=Ψ0;
jj=1

# The operators we have to start off with
Ψ=Ψ.*Uxh

Psif=ft*Ψ
Psif=Psif.*Uf    
Ψ=ift*Psif

#@gif for ii in 1:nt
for ii in 1:nt

    Ψ=Ψ.*Ux

    Psif=ft*Ψ
    Psif=Psif.*Uf    
    Ψ=ift*Psif   

    if ii%nmeas == 0

        # Every time we measure, we have to finish with a U_x half time step
        Ψt=Ψ.*Uxh

        c[jj,1]=abs(sum( conj(Ψt).*ϕl )) *dx
        c[jj,2]=abs(sum( conj(Ψt).*ϕr )) *dx


        jj+=1
    end

    #plot(x[21:141],Vx[21:141]/6, label="Vx scaled")
    #plot!(x,abs(conj(Ψt).*Ψt), label="Wavefunction")
    #plot!(xlabel="x", ylabel="Ψ",
    #    plot_title="Wavefunction evolution")

end    
#end every nmeas

Ψ=Ψ.*Uxh;
```
{% include image.html img="M4/Images/TimeEvolution/evolution_plots.gif" title="Time Evolution" caption="Time Evolution in a double well potential.  This gif was generated by `Plots.jl`"%}

{% include image.html img="M4/Images/TimeEvolution/doublewell.gif" title="Time Evolution" caption="Time Evolution in a double well potential.  This gif was generated splicing together PyPlot png's."%}

```julia
plot(c[:,1].^2,label="Left Prob")
plot!(c[:,2].^2, label="Right Prob")
plot!(xlabel="x",ylabel="Probability",
    plot_title="Rabi Oscillations for a Double Harmonic Oscillator")
```

<iframe src="/M4/Images/TimeEvolution/rabioscil.html"  style="border:none; background: #ffffff"
width="700px" height="450px"></iframe>

The faster one can perform Fourier Transforms, the faster one can carry out this algorithm.  Therefore, scientists, such as [2] and [3], will use multiple cores or GPU's.  

In addition to real-time evolution, algorithms like this can determine the ground state of an arbitrary system by imaginary time evolution.  Soon, I will take the work covered here and look at this aspect of the algorithm.  

I'd particularly like to thank Dr. Lee J. O'Riordan [3] for sharing his experience on this topic.

<a name="citation1">[1]</a> Glowinski, Roland, Stanley J. Osher, and Wotao Yin, eds. Splitting Methods in Communication, Imaging, Science, and Engineering. Springer, 2017.

[2] Heiko Bauke and Christoph H. Keitel. Accelerating the Fourier split operator method via graphics processing unit. Computer Physics Communications, 182(12):2454–2463 (2011)

[3] Lee James O'Riordan et al., GPUE: Phasegineering release, Zenodo. (2016) [https://github.com/mlxd/GPUE](https://github.com/mlxd/GPUE) DOI:10.5281/zenodo.57968
