---
title: Jacobi Transformation of a Symmetric Matrix
layout: post
comments: True
category: Numerics
tags: ED
author: Christina C. Lee
description:
---

<i>Based on Numerical Recipes in C++, Sec 11.1</i>



So you want to diagonalize a matrix, do you?
Well, if you have a tiny symmetric matrix, you REALLY want to write up the algorithm by hand, and don't want to spend much time trying to understand the algorithm, then you have come to the right place.

Otherwise, use LAPACK/BLAS to call a highly optimized routine that can work extremely quickly on large matrices. Julia has those libraries built in already. Even if you do call those matrices, you can make them work better by understanding what's going on underneath the hood, which is why we are going through this now.

Start with a base <i>rotation matrix</i> of the Form
$$
    P_{pq} =
    \begin{pmatrix}
           1& &  &  & & & & 0 \\
           & \ddots & & & &  &  \\
            & & c & \cdots & s & & \\
                & &\vdots& 1 & \vdots & & \\
               & & -s & \cdots & c & &  \\
               & & & & & \ddots &  \\
               0 & & & &  & & 1\\
    \end{pmatrix}
$$

From our starting arbitrary symmetric $$A$$,

$$
A^T = A
$$


we will run a series of transformations,

$$
A^{\prime}= P^{T}\_{pq} \cdot A \cdot P\_{pq}
$$

where each iteration brings $$A$$ closer to diagonal form.  Thus in our implementing our algorithm, we need to determine two things

* The values of $$c$$ and $$s$$
* The pattern of sweeping $$p$$ and $$q$$

And in the end we will need to finally determine if this actually converges, and if has any sort of efficiency.

So lets expand one transformation, and we if we can solve for $c$ and $s$.

$$
a^{\prime}_{rp}  = c a_{rp} - s a_{rq}
$$

$$
a^{\prime}_{rq}  = c a_{rq} + s a_{rp}
$$

$$
a^{\prime}_{pp}  = c^2 a_{pp} + s^2 a_{qq} -2 sc a_{pq}
$$

$$
a^{\prime}_{qq}  = s^2 a_{qq} + c^2 a_{qq} + 2sc a_{pq}
$$

$$
a^{\prime}_{pq}  = \left( c^2-s^2 \right) a_{pq} + sc \left(a_{pq} - a_{qq} \right)
$$

## Determining $s$ and $c$
Given we specifically want $a^{\prime}_{pq}$ to be zero, we re-arrange the last equation,

$$
        \frac{c^2-s^2}{2 sc} = \frac{a_{pq}-a_{qq}}{2 a_{pq}} = \theta
$$

At first glance, this equation might not look easier to solve for $s$ or $c$.  At second glance either. We define a new parameter $t = s/c$, which now makes the equation,

$$
\frac{1-t^2}{2 t} = \theta \;\;\;\; \implies \;\;\; t^2 -2 \theta t -1=0,
$$

now quite easily solvable by our friendly quadratic formula.  Though the book does recommend using a form that pulls out smaller roots through

$$
t=\frac{\text{sgn}( \theta )}{| \theta | + \sqrt{ \theta ^2 + 1} }.
$$

Then reverse solve back to

$$
c=\frac{1}{\sqrt{t^2+1}} \;\;\; s=tc
$$

Though we could use the expressions above, if we simplify them with our new expressions for $c$ and $s$ analytically, we reduce computational load and round off error. These new expressions are

$$
a^{\prime}_{pq}  = 0
$$

$$
a^{\prime}_{qq}  = a_{qq} + t a_{qp}
$$

$$
a^{\prime}_{pp} = a_{pp} - t a_{pq}
$$

$$
a^{\prime}_{rp} = a_{rp} - s \left( a_{rq} +\tau a_{rp} \right)
$$

$$
a^{\prime}_{rq} = a_{rq} + s \left( a_{rp} -\tau a_{rq} \right)
$$

with the new variable

$$
\tau = \frac{s}{1+c}
$$

## Convergence

The sum of the squares of the off diagonal elements, choosen in either upper or lower triangles arbitrarily,

$$
S=\sum\limits_{r < s} |a_{rs}|^2
$$

## Eigenvectors

By forming a product of every rotation matrix, we also come to approximate the matrix $V$ where

$$
D = V^{T} \cdot A \cdot V
$$

and $D$ is the diagonal form of $A$.  $V$ is computed through iterative computation

$$
V^{\prime} = V \cdot P_i
$$

$$
v^{\prime}_{rs} = v_{rs}
$$

$$
v^{\prime}_{rp} = c v_{rp} - s v_{rq}
$$

$$
v^{\prime}_{rq} = s v_{rp} + c v_{rq}
$$

### Enough with the talking! LETS COMPUTE STUFF

{% highlight julia %}
using LinearAlgebra
{% endhighlight %}


{% highlight julia %}
# First, Lets make our nice, helpful functions

## A function to look at the convergence
function convergence(A::Array)
    num=0.0
    l=size(A)[1]
    for ii in 1:(l-1)
        for jj in (ii+1):l ## just looking at the lower triangle
            num+=A[ii,jj]^2
            #println(ii,' ',jj,' ',num,' ',A[ii,jj])
        end
    end
    return num
end
{% endhighlight %}


This makes a matrix easier to look at than when its filled
with 1.043848974e-12 everywhere

{% highlight julia %}
# This makes a matrix easier to look at when its filled 
# with 1.043848974e-12 everywhere
function roundmatrix(A::Array,rtol::Real)
    Ap=copy(A)
    for ii in 1:length(A)
        if abs(Ap[ii])<rtol 
            Ap[ii]=0
        end
    end
    return Ap;
end
{% endhighlight %}




{% highlight julia %}
## Here we create a random symmetric matrix
function makeA(n)
    A=randn(n,n);
    for ii in 1:n
        A[ii,1:ii]=transpose(A[1:ii,ii]) 
    end
    V=Matrix{Float64}(I,n,n) #initializing the orthogonal transformation
    return A,copy(A),V
end
## One A returned will be stored to compare initial and final
{% endhighlight %}

Now on to the rotations!

 We don't always want to compute the eigenvectors, so those are in the optional entries slot.
Both tell the function to compute the vectors with `computeV=true`
and input the `V=V` after the semicolon.


{% highlight julia %}
function Rotate(A::Array,p::Int,q::Int; computeV=false, V::Array=Matrix{Float64}(I,1,1) )
    θ=(A[q,q]-A[p,p])/(2*A[p,q]);
    t=sign(θ)/(abs(θ)+sqrt(θ^2+1));
    
    c=1/sqrt(t^2+1)
    s=t*c
    τ=s/(1+c)
    
    l=size(A)[1]
    Ap=copy(A[:,p])
    Aq=copy(A[:,q])
    for r in 1:l
        A[r,p]=Ap[r]-s*(Aq[r]+τ*Ap[r]) 
        A[r,q]=Aq[r]+s*(Ap[r]-τ*Aq[r])
        
        A[p,r]=A[r,p]
        A[q,r]=A[r,q]
    end
    A[p,q]=0
    A[q,p]=0
    A[p,p]=Ap[p]-t*Aq[p]
    A[q,q]=Aq[q]+t*Aq[p]
    
    if computeV==true
        Vp=copy(V[:,p])
        Vq=copy(V[:,q])    
        for r in 1:l   
            V[r,p]=c*Vp[r]-s*Vq[r]
            V[r,q]=s*Vp[r]+c*Vq[r]
        end
        return A,V
    else
        return A;
    end
end
{% endhighlight %}

This function performs one sweep

{% highlight julia %}

function Sweep(A;compV=false,V=Matrix{Float64}(I,1,1))
    n=size(A)[1]
    for ii in 2:n
        for jj in 1:(ii-1) ## Just over one triangle
            if compV==false
                A=Rotate(A,ii,jj)
            else
                A,V=Rotate(A,ii,jj;computeV=true,V=V);
            end
        end
    end
    
    if compV==false
        return A
    else
        return A,V
    end 
end
{% endhighlight %}



Just creating some size of matrix

{% highlight julia %}
A,A0,V=makeA(5);
{% endhighlight %}


{% highlight julia %}
## keep evaluating for a couple iterations
## watch how it changes
A,V=Sweep(A;compV=true,V=V);
display(roundmatrix(A,1e-10))
display(A)
display(V)
display(convergence(A))
{% endhighlight %}


This output is after several sweeps

        5×5 Array{Float64,2}:
        2.75547   0.0     0.0      0.0        0.0    
        0.0      -1.0683  0.0      0.0        0.0    
        0.0       0.0     1.20611  0.0        0.0    
        0.0       0.0     0.0      0.734884   0.0    
        0.0       0.0     0.0      0.0       -1.48581

        5×5 Array{Float64,2}:
        2.75547       1.03898e-31  6.65473e-46   -1.20864e-52  -4.52913e-60 
        1.03898e-31  -1.0683       4.25018e-46   -2.9053e-69    1.37672e-80 
        6.65473e-46   4.25018e-46  1.20611        1.91399e-75   5.09273e-116
        -1.20864e-52  -2.9053e-69   1.91399e-75    0.734884      0.0         
        -4.52913e-60   1.37672e-80  5.09273e-116   0.0          -1.48581     

        5×5 Array{Float64,2}:
        0.702867    0.340672  -0.613045  0.0266825  -0.115695 
        -0.0750998   0.373435   0.239379  0.773343   -0.446705 
        0.597584   -0.208098   0.484961  0.308218    0.519041 
        -0.0642891  -0.754388  -0.459432  0.458861   -0.0716542
        -0.37296     0.363431  -0.347289  0.309317    0.715913 

        1.0794844600380612e-62



Compare the Optimized LAPLACK routine to your results

{% highlight julia %}
eigen(A0)
{% endhighlight %}

        Eigen{Float64,Float64,Array{Float64,2},Array{Float64,1}}
        eigenvalues:
        5-element Array{Float64,1}:
        -1.4858101513857376
        -1.068298321190819 
        0.7348837698415409
        1.206106933351782 
        2.755472216052164 
        eigenvectors:
        5×5 Array{Float64,2}:
        -0.115695   -0.340672  -0.0266825   0.613045   0.702867 
        -0.446705   -0.373435  -0.773343   -0.239379  -0.0750998
        0.519041    0.208098  -0.308218   -0.484961   0.597584 
        -0.0716542   0.754388  -0.458861    0.459432  -0.0642891
        0.715913   -0.363431  -0.309317    0.347289  -0.37296  




{% highlight julia %}
## A good check to make sure V is an orthonomal transformation
roundmatrix(V*A*transpose(V)-A0,1e-12)
{% endhighlight %}

          5x5 Array{Float64,2}:
           0.0  0.0  0.0  0.0  0.0
           0.0  0.0  0.0  0.0  0.0
           0.0  0.0  0.0  0.0  0.0
           0.0  0.0  0.0  0.0  0.0
           0.0  0.0  0.0  0.0  0.0


How long does it take to make a Sweep?
How much memory will the computation take?
This is dependent on how large the matrix is, and determines whether or not we
want to use this algorithm.

{% highlight julia %}

A,A0,V=makeA(10);
@time Sweep(A);
A,A0,V=makeA(20);
@time Sweep(A);
A,A0,V=makeA(100);
@time Sweep(A);
{% endhighlight %}

        0.000016 seconds (230 allocations: 32.594 KiB)
        0.000069 seconds (955 allocations: 196.188 KiB)
        0.008129 seconds (24.75 k allocations: 17.372 MiB, 27.48% gc time)


In addition to time per sweep, we need to know how many sweeps we need to run. So again we run it on a 10x10, 20x20, and 100x100. The efficiency of the algorithm would get a lot worse if we have to sweep the 100x100 a bunch of times.


{% highlight julia %}
A10,Ap10,V=makeA(10);
A20,Ap20,V=makeA(20);
A100,Ap100,V=makeA(100);
nsweep=collect(1:7);
conv10=zeros(7)
conv20=zeros(7)
conv100=zeros(7)
for i in nsweep
    A10=Sweep(A10)
    A20=Sweep(A20)
    A100=Sweep(A100)
    conv10[i]=convergence(A10)
    conv20[i]=convergence(A20)
    conv100[i]=convergence(A100)
end

[nsweep conv10/10 conv20/20 conv100/100]
{% endhighlight %}



        7×4 Array{Float64,2}:
        1.0  1.74923       2.64638      14.7488     
        2.0  0.0945499     0.422473      2.80609    
        3.0  0.000314227   0.0162891     0.399226   
        4.0  6.31792e-10   1.09268e-5    0.0356924  
        5.0  3.62048e-22   1.18607e-11   0.000598666
        6.0  9.38425e-48   7.94096e-26   3.28477e-7 
        7.0  1.14895e-112  1.23362e-55   6.11775e-13



Well, so we've seen how to do one form of exact diagonalization that works, but doesn't scale very well up to 100x100 matrices.  So stay tuned for the Householder method, hopefully coming up soon.

Until then, happy computing :)
