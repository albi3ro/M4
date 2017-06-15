---
title: Phase Transitions
layout: post
comments: True
category: Prerequisites Required
author: Christina C. Lee
image: "https://albi3ro.github.io/M4/Images/PhaseTransitions/output_23_0.png"
tags: [Magnet]
description:  A look at phase transitions in magnets using a Metropolis-Hastings Algorithm, specifically at the square lattice Ising model.
---

Post Prerequisites: [Monte Carlo Calculation of pi]({{base.url}}/M4/numerics/Monte-Carlo.html), [Monte Carlo Markov Chain]({{base.url}}/M4/numerics/MCMC.html) [Monte Carlo Ferromagnet]({{base.url}}/M4/prerequisites/Monte-Carlo-Ferromagnet.html)

Prerequisites: Statistical Mechanics


In Monte Carlo Ferromagnets, we only looked at the state of a magnet at one temperature.  While I displayed how it looked at both a high temperature and a low-temperature situation, that had to be done manually.  Today we are going to automate our study of temperature's effect.

First off, what's the theory?

The standard apparatus for phase transitions is <b>Landau Theory</b>.  In Landau theory, we assume that we have an <b>order parameter</b>.

Order parameters have
* local physical meaning
* global observational consequences

For example, in magnets we use the local magnetization $m$; in freezing fluids, the difference in density $\rho$.  The value of an order parameter should undergo a jump across a phase transition, usually from a zero to a non-zero value.

To construct Landau theory, we postulate that we can create a <b>Free Energy</b> for the system that is
* continuous in the order parameter
* respects the symmetries of the system

Given those constraints, we get
\begin{equation}
F=f_0(T) + \alpha_0 (T-T_c) m^2+ \frac{1}{2} \beta m^4 - h m + ...
\end{equation}
From this simple equation, we can derive many of the qualitative properties of phase transitions.  I won't repeat that here, but you can look it up in a Statistical Mechanics textbook.

Non-zero order parameters also represent symmetry breaking.  In an unmagnetized magnet, all directions are equivalent, but once it magnetizes, it chooses a preferred direction.  That direction is arbitrary up until the magnet spontaneously chooses it.

Back in the early days of cars, streets didn't have lanes and a specified side for cars to drive on.  Driving was symmetric.   When two cars approached each other on the same road, the drivers had to figure out how to avoid a collision, either by everyone going to their left or everyone going to their right.  Eventually, all the cars in a specific area picked up one convention. In the United States, we choose the right side of the road.  In Japan, they chose left.  Neither is wrong; it was pure chance.  The solution breaks the symmetry of the problem.

{% include image.html img="M4/Images/PhaseTransitions/symmetrybreaking.svg" title="Symmetry Breaking" caption="Though both potentials are symmetric, the low energy solution on the right breaks the symmetry." %}

### Specific Heat and Magnetic Susceptibility

In addition to magnetization and energy, today I introduce the physical observables specific heat and magnetic susceptibility.  These numbers describe how the old quantities change with temperature.

Let's take a look at the first one, <b>Specific Heat</b>:

$$
C_h=\frac{d \langle E \rangle}{dT}
$$

While the above is the formal definition, we can calculate a more convenient expression by plugging in the Maxwell-Boltzmann Distribution:

$$
\frac{ d
\frac{\sum_{i} E_{i} e^{-\beta E_i}}{ \sum_j e^{-\beta E_j} }
}{dT}=\beta^2 \left( \langle E^2 \rangle - \langle E \rangle ^2 \right)
$$

This equation we use to calculate the specific heat, as it uses the values calculated at only one temperature, eliminating the need for a numerical derivative.  For each temperature, we will need to calculate both the expectation of the Energy and the expectation of the Energy Squared.

While in general, a material may be difficult to describe, near phase transitions certain properties diverge according to corresponding <b>critical exponents</b>.  These fall out of the phenomenology of Landau theory, which you can find more about in a textbook.

These critical exponents obey <b>universality</b>.  Phase transitions in the same class will follow the same critical exponents.  Determining the exponent for the 2D square lattice through our simulations also determines the exponent for situations like the mixing of two materials or the percolation problem.

We define our specific heat critical exponent this way,

$$
C_h \propto \left( T- T_c \right)^{-\alpha_{\pm}}.
$$

The divergent scaling of most properties around a critical temperature depends on the divergence of the <b>correlation length</b> $\xi$, the size of fluctuations in the system.  At the critical temperature, fluctuations exist at <i>every length scale</i>.   Once we get close enough to the transition temperature for the fluctuations to become greater than our finite simulation lattice, our computationed numbers no longer agree with infinite-lattice-assumed theory.  The divergence gets modified by some function of the ratio between the correlation length and our system size,

$$
C_h \propto \left( T - T_c \right)^{-\alpha_{\pm}} g\left( \frac{L}{\xi (T)} \right).
$$

Magnetic Susceptibility fulfills the same general role as specific heat, but for magnetization instead of energy.

$$
\chi =\frac{d \langle M \rangle}{dT}=\beta^2 \left( \langle M^2 \rangle - \langle M \rangle ^2 \right)
$$

Again, this property also has a corresponding critical exponent.

$$
\chi \propto \left( T - T_c \right)^{-\gamma_{\pm}}
$$

## Exact Solution

The 2D Square Ising model possesses an exact solution, first discovered by Lars Onsager in 1944, but since rewritten in a variety of different ways.  I won't go over how to achieve the solution, but I will show the results for comparison.

The critical temperature:

$$
J \beta_c= \frac{\ln \left( \sqrt{2}+1 \right)}{2} \approx .4406
$$

$$
T_c \approx 2.269 J
$$

The magnetization curve when $T<T_c$ also $\beta > \beta_c$

$$
M=\left(1-\sinh^{-4} \left(2 \beta J \right) \right)^{1/8}
$$

And finally, the much more complicated expression for internal energy per site,

$$
k=\sinh (2 \beta J)^{-2}
$$

$$
U=-J \coth (2 \beta J) \Bigg(
1 + \frac{2}{\pi} \left( 2 \tanh^2 (2 \beta J) - 1 \right)
$$

$$
\;\;\;\;\;\; \int_0^{\pi/2} \frac{1}{\sqrt{1-4k (1+k)^{-2} \sin^2 (\theta)}}
d \theta
\Bigg)
$$

If we can achieve these results to perfect accuracy on pen and paper, why do we even bother simulating them in a computer?  If you try some other lattices and coupling constants, you might soon find out why.  We want to make <i>sure</i> our code is running well on something where we know the results before venturing into unknown territory.

To start off our coding, let's just type those equations up.


```julia
# The critical inverse temperature
betac=log(1+sqrt(2))/2

# The exact magnetization
# Note: Does not take arrays
function Mexact(beta::Float64)
    if beta>betac
        return (1-sinh(2*beta*J)^(-4.))^(1./8.)
    else
        return 0
    end
end

# The exact energy
# Note: Does not take arrays
function Eexact(beta::Float64)
    k=sinh(2*beta*J)^(-2.)

    # the inside of the integral.  Define a function so we can use quadgk
    insides(θ)=(1-4*k*(1+k)^(-2)*sin(θ)^2)^(-.5)
    #we run a numerical integration
    integ=quadgk(insides,0,π/2)[1]

    return -J*coth(2*beta*J)*(1+2/π*(2*tanh(2*beta*J)^2-1)*integ )
end
```



Load our packages.

`Lattices.jl` is the same class I used in the previous post.  Well, I got around to updating the Julia syntax.  So mostly the same.


```julia
push!(LOAD_PATH,"../Packages")
using Lattices;
using PyPlot
```

Here we set up our lattice.

If looking for a mind bender, look at an `"Checkerboard"`, with an anti-ferromagnetic coupling constant, like `J=-1`.  Otherwise, just stick with a square.


```julia
## Define l here
l=20;

lt=MakeLattice("Square",l);
S=ones(Int8,l,l);  #Our spins
dt=1/(lt.N);
```

Same as last time.

These functions calculate properties of our lattice.


```julia
# The energy contribution of just one site
function dE(i::Int)
    Eii=0;
    for j in 1:lt.nnei
        Eii+=S[lt.neigh[i,j]];
    end
    Eii*=-J*S[i];  # we are computing J sz_i sz_j for one i
    return Eii;
end
# The energy of the entire lattice
function E()
    Evar=0;
    for k in 1:lt.N
        Evar+=.5*dE(k);
    end
    return Evar;
end
# The magnetization of the entire lattice
function M()
    Mvar=0;
    for k in 1:lt.N
        Mvar+=S[k];
    end
    return Mvar;
end
"defined functions"
```



### Adjustable Parameters


```julia
beta0=.1;
betaf=1;
dbeta=.1;

J=1;
t=10000;
nskip=10;   # don't measure every sweep= better decorrelation
"Parameters set"
```


```julia
nmeas=Int64(t/nskip); # how many times we will measure
betas=collect(beta0:dbeta:betaf)
"done"
```

I took what we ran last time and wrapped it into a function, `MCMCMagnet`.  We aren't looking at spin configurations this time; we only want the final measureables, which get returned at the end.


```julia
function MCMCMagnet(beta::Float64)
    tm=1; #Our measurement time step

    Ma=Array{Int32}(nmeas); # our magnetization measurements
    Ea=Array{Int32}(nmeas); # our energy measurements
    Ma2=Array{Int32}(nmeas); # magnetization squared
    Ea2=Array{Int32}(nmeas); # energy squared

    for ti in 1:t
        for j in 1:lt.N
            i = rand(1:lt.N); #Choosing a random site
            de=dE(i);
            if(de>0 || rand()<exp(2*beta*de) )
                S[i]=-S[i]; #Switch the sign
            end
        end
        if isapprox(mod(ti,nskip),0)
            Ma[tm]=M();
            Ma2[tm]=Ma[tm]^2;

            Ea[tm]=E();
            Ea2[tm]=Ea[tm]^2;

            tm+=1;

        end
    end
    Mave=mean(Ma);
    Mstd=std(Ma)/lt.N;
    Eave=mean(Ea);
    Estd=std(Ea)/lt.N;

    E2ave=mean(Ea2);
    M2ave=mean(Ma2);

    Ch=beta^2*(E2ave-Eave^2)/lt.N;
    χ=beta*(M2ave-Mave^2)/lt.N;

    return Mave/lt.N,Mstd,Eave/lt.N,Estd,Ch,χ
end
```



# The Temperature Loop

This cell looks pretty simple: initialization and one for loop.  But this is where everything ties together.


```julia
Mm=zeros(betas)
Mstd=zeros(betas)
Ee=zeros(betas)
Estd=zeros(betas)
Ch=zeros(betas)
χ=zeros(betas)
for ii in 1:length(betas)
    Mm[ii],Mstd[ii],Ee[ii],Estd[ii],Ch[ii],χ[ii]=MCMCMagnet(betas[ii])
    println("beta: ",betas[ii],"\tM: ",Mm[ii],"\tE: ",Ee[ii])
end
```

    beta: 0.1	M: 0.00039	E: -0.20528
    beta: 0.2	M: 0.003225	E: -0.42512
    beta: 0.3	M: 0.001325	E: -0.6998099999999999
    beta: 0.4	M: 0.003525	E: -1.11878
    beta: 0.5	M: 0.8965900000000001	E: -1.73719
    beta: 0.6	M: 0.973875	E: -1.9091300000000002
    beta: 0.7	M: 0.989965	E: -1.96268
    beta: 0.8	M: 0.996135	E: -1.98525
    beta: 0.9	M: 0.99832	E: -1.99347
    beta: 1.0	M: 0.99934	E: -1.99739


That was the temperature loop for the MCMC method.

This is the temperature loop to calculate the exact numbers.

```julia
Mex=zeros(betas)
Eex=zeros(betas)
for ii in 1:length(betas)
    Mex[ii]=Mexact(betas[ii])
    Eex[ii]=Eexact(betas[ii])
end
```


### Storing and Restoring Data

Since running these simulations take slightly longer than some of my previous posts, I've outputed data files of the results of three different sizes lattices (l=10, l=20, l=50), for runs of length 50,000.  All have the same `betas` variable. Feel free to take a look at them, or make your own data :)


```julia
#writedlm("PhaseTransitionsL20T50000.dat",[betas Mm Mstd Ee Estd Ch χ])
```


```julia
temp=readdlm("PhaseTransitionsL50T50000.dat")
betas=temp[:,1];
Mm50=temp[:,2];
Mstd50=temp[:,3];
Ee50=temp[:,4];
Estd50=temp[:,5];
Ch50=temp[:,6];
χ50=temp[:,7];

temp=readdlm("PhaseTransitionsL10T50000.dat")
Mm10=temp[:,2];
Mstd10=temp[:,3];
Ee10=temp[:,4];
Estd10=temp[:,5];
Ch10=temp[:,6];
χ10=temp[:,7];

temp=readdlm("PhaseTransitionsL20T50000.dat")
Mm20=temp[:,2];
Mstd20=temp[:,3];
Ee20=temp[:,4];
Estd20=temp[:,5];
Ch20=temp[:,6];
χ20=temp[:,7];
```




# Plotting Section


```julia
plot(betas,Mex,label="Exact",linewidth="5")
plot(betas,abs(Mm10),label="l=10")
plot(betas,abs(Mm20),label="l=20")
plot(betas,abs(Mm50),label="l=50")

xlabel("Inverse Temperature β")
ylabel("Magnetization")
legend()
title("2D Square Ferromagnet Phase Transition")
```


{% include image.html img="M4/Images/PhaseTransitions/output_23_0.png" title="Magnetization" caption="The order parameter magnetization jumps at the critical temperature $\beta_c\approx.44$." %}




```julia
plot(betas,Eex,label="Exact",linewidth=5)
plot(betas,Ee10,label="l=10")
plot(betas,Ee20,label="l=20")
plot(betas,Ee50,label="l=50")
xlabel("Inverse Temperature β")
ylabel("Magnetization")
legend()
title("2D Square Ferromagnet Phase Transition")
```

{% include image.html img="M4/Images/PhaseTransitions/output_24_0.png" title="Energy" caption="The Energy curve for a variety of sizes and the exact infinite lattice." %}


```julia
plot(betas,Ch10,label="l=10")
plot(betas,Ch20,label="l=20")
plot(betas,Ch50,label="l=50")
xlabel("Inverse Temperature β")
ylabel("Specific Heat")
legend()
title("Specific Heat Divergence with Finite Size Scaling")
```

{% include image.html img="M4/Images/PhaseTransitions/output_25_0.png" title="Specific Heat" caption="The specific heat diverges around the critical temperature, but a finite simulation lattice prevents a true divergence from occurring.  The larger the lattice, the more the specific heat diverges." %}



I presented two different equivalent formulas for the specific heat.  Let's check to make sure their equivalent.

For the formal defination, we'll take a numeric derivative with respect to $\beta$,
\begin{equation}
\frac{d \langle E \rangle}{dT} \approx \frac{ \Delta \langle E \rangle}{\Delta \beta}\frac{d \beta}{d T}
\end{equation}

Note: I corrected the scaling by hand.


```julia
Ch2=(Ee50[1:end-1]-Ee50[2:end]).*betas[2:end].^2

plot(betas[2:end],Ch2*35,label="dE/dT")
plot(betas,Ch50,label="beta2(<E2>-<E>2)")
xlabel("Temperature")
ylabel("Specific Heat")
legend()
title("Comparing Specific Heat Calculations")

```


{% include image.html img="M4/Images/PhaseTransitions/output_27_0.png" title="Specific Heat Calculations" caption="Both the derivative formula and the expectation value formula for the specific heat give the same result, after correcting for a scaling factor." %}




```julia
plot(betas,χ10,label="l=10")
plot(betas,χ20,label="l=20")
plot(betas,χ50,label="l=50")
xlabel("Inverse Temperature β")
ylabel("Susceptibility")
legend()
title("Magnetic Susceptibility with Finite Size Scaling")
```


{% include image.html img="M4/Images/PhaseTransitions/output_28_0.png" title="Magnetic Susceptibility" caption="The magnetic susceptibility, unlike specific heat, is not symmetric around the critical temperature.  The apparent value for critical temperature also changes with system size, approaching the infinite lattice value." %}




```julia
M2=(Mm50[1:end-1]-Mm50[2:end]).*(betas[1:end-1]).^2
plot(betas[1:end-1],M2*750,label="dM/dT")
plot(betas,χ50,label="(<M2>-<M>2)beta2")
legend()
xlabel("Inverse Temperature β")
ylabel("Specific Heat")
title("Comparing Methods of Computing Susceptibility")
```


{% include image.html img="M4/Images/PhaseTransitions/output_29_0.png" title="Magnetic Susceptibility Methods" caption="Both the derivative method and the expectation value method give approximately the value for the susceptibility." %}



That's all for now.  Plenty more to talk about on this subject, but hopefully that can keep you occupied.
