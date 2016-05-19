---
title: Jacobi Transformation of a Symmetric Matrix
layout: post
comments: True
category: Numerics
tags: [Exact Diagonalization]
author: Christina C. Lee
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


{% highlight MATLAB %}
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

{% highlight MATLAB %}
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




{% highlight MATLAB %}
## Here we create a random symmetric matrix
function makeA(n::Int)
    A=randn(n,n);
    for ii in 1:n
        A[ii,1:ii]=transpose(A[1:ii,ii])
    end
    V=eye(n) #initializing the orthogonal transformation
    return A,copy(A),V
end
## One A returned will be stored to compare initial and final
{% endhighlight %}

Now on to the rotations!

 We don't always want to compute the eigenvectors, so those are in the optional entries slot.
Both tell the function to compute the vectors with `computeV=true`
and input the `V=V` after the semicolon.


{% highlight MATLAB %}
function Rotate(A::Array,p::Int,q::Int; computeV=false, V::Array=eye(1))
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

{% highlight MATLAB %}

function Sweep(A;compV=false,V=eye(1))
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

{% highlight MATLAB %}
A,A0,V=makeA(5);
{% endhighlight %}


{% highlight MATLAB %}
## keep evaluating for a couple iterations
## watch how it changes
A,V=Sweep(A;compV=true,V=V);
roundmatrix(A,1e-10),A,V,convergence(A)
{% endhighlight %}


This output is after several sweeps

        (
        5x5 Array{Float64,2}:
         -1.59942  0.0       0.0       0.0      0.0
          0.0      1.03678   0.0       0.0      0.0
          0.0      0.0      -0.823094  0.0      0.0
          0.0      0.0       0.0       3.09433  0.0
          0.0      0.0       0.0       0.0      1.3409,

        5x5 Array{Float64,2}:
         -1.59942      5.1314e-30    2.32594e-36  -9.54088e-49  -1.22782e-53
          5.1314e-30   1.03678       2.65014e-38   9.13791e-56   6.64996e-67
          2.32594e-36  2.65014e-38  -0.823094     -9.56652e-61   2.08002e-92
         -9.54088e-49  9.13791e-56  -9.56652e-61   3.09433       0.0
         -1.22782e-53  6.64996e-67   2.08002e-92   0.0           1.3409     ,

        5x5 Array{Float64,2}:
          0.0537334   0.0599494  -0.735228  0.139      0.658511
          0.310018    0.612957   -0.14049   0.611348  -0.367001
          0.759653   -0.475834    0.264118  0.282571   0.216575
         -0.480405   -0.546544   -0.132383  0.644217  -0.194831
         -0.305189    0.30913     0.593648  0.33477    0.588905,

        2.6331310238375346e-59)



Compare the Optimized LAPLACK routine to your results

{% highlight MATLAB %}
eig(A0)
{% endhighlight %}

          ([-1.599424470672961,-0.8230937166650976,1.0367806031602211,
          1.3408963512476402,3.0943321944116593],
          5x5 Array{Float64,2}:
           -0.0537334   0.735228   0.0599494  -0.658511  -0.139
           -0.310018    0.14049    0.612957    0.367001  -0.611348
           -0.759653   -0.264118  -0.475834   -0.216575  -0.282571
            0.480405    0.132383  -0.546544    0.194831  -0.644217
            0.305189   -0.593648   0.30913    -0.588905  -0.33477 )




{% highlight MATLAB %}
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

{% highlight MATLAB %}

A,A0,V=makeA(10);
@time Sweep(A);
A,A0,V=makeA(20);
@time Sweep(A);
A,A0,V=makeA(100);
@time Sweep(A);
{% endhighlight %}

      0.000028 seconds (320 allocations: 30.469 KB)
      0.000099 seconds (1.33 k allocations: 187.266 KB)
      0.007413 seconds (34.66 k allocations: 17.448 MB, 14.20% gc time)


In addition to time per sweep, we need to know how many sweeps we need to run. So again we run it on a 10x10, 20x20, and 100x100. The efficiency of the algorithm would get a lot worse if we have to sweep the 100x100 a bunch of times.


{% highlight MATLAB %}
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




    7x4 Array{Float64,2}:
     1.0  1.10944       2.43759      14.6644
     2.0  0.105628      0.312076      2.87182
     3.0  0.000265288   0.017073      0.498082
     4.0  6.64324e-9    0.000119472   0.0390564
     5.0  4.05463e-18   3.56679e-11   0.00133833
     6.0  3.17274e-42   1.96318e-23   6.07661e-7
     7.0  6.76289e-110  4.07871e-49   3.98102e-13



Well, so we've seen how to do one form of exact diagonalization that works, but doesn't scale very well up to 100x100 matrices.  So stay tuned for the Householder method, hopefully coming up soon.

Until then, happy computing :)
