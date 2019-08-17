---
title: Atomic Orbitals Pt. 2
layout: post
comments: True
category: Prerequisites Required
tags: [Quantum]
author: Christina C. Lee
description: We try out the GLVisualize.jl package to view the hydrogen orbitals computed in the last post.
image: "https://albi3ro.github.io/M4/Images/Orbitals2/3p_surface.png"
---
<b>Prerequiresites:</b> Quantum Mechanics course

If you haven't read it already, check out [Atomic Orbitals Pt. 1]({{base.url}}/M4/Atomic-Orbitals.html).  Today, we try and make some prettier pictures. GLVisualize is quite a beautiful package, but not entirely the easiest to use at this point with some not so consistent documentation.

To add this package:


{% highlight julia %}
Pkg.add("GLVisualize")
{% endhighlight %}

and test with:


{% highlight julia %}
Pkg.test("GLVisualize")
{% endhighlight %}

But, other steps may be necessary to get the package working.  On a Mac, I was required to install the Homebrew.jl package.


{% highlight julia %}
#Pkg.update();
#Pkg.add("GSL");
using GSL;
using GLVisualize;
{% endhighlight %}


{% highlight julia %}
a0=1; #for convenience, or 5.2917721092(17)×10−11 m

# The unitless radial coordinate
ρ(r,n)=2r/(n*a0);

#The θ dependence
function Pmlh(m::Int,l::Int,θ::Real)
    return (-1.0)^m *sf_legendre_Plm(l,m,cos(θ));
end

#The θ and ϕ dependence
function Yml(m::Int,l::Int,θ::Real,ϕ::Real)
    return  (-1.0)^m*sf_legendre_Plm(l,abs(m),cos(θ))*e^(im*m*ϕ)
end

#The Radial dependence
function R(n::Int,l::Int,ρ::Real)
    if isapprox(ρ,0)
        ρ=.001
    end
     return sf_laguerre_n(n-l-1,2*l+1,ρ)*e^(-ρ/2)*ρ^l
end

#A normalization: This is dependent on the choice of polynomial representation
function norm(n::Int,l::Int)
    return sqrt((2/n)^3 * factorial(n-l-1)/(2n*factorial(n+l)))
end

#Generates an Orbital Funtion of (r,θ,ϕ) for a specificied n,l,m.
function Orbital(n::Int,l::Int,m::Int)
    if (l>n)    # we make sure l and m are within proper bounds
        throw(DomainError())
    end
    if abs(m)>l
       throw(DomainError())
    end
    psi(ρ,θ,ϕ)=norm(n, l)*R(n,l,ρ)*Yml(m,l,θ,ϕ);
    return psi
end

#We will calculate is spherical coordinates, but plot in cartesian, so we need this array conversion
function SphtoCart(r::Array,θ::Array,ϕ::Array)
    x=r.*sin(θ).*cos(ϕ);
    y=r.*sin(θ).*sin(ϕ);
    z=r.*cos(θ);
    return x,y,z;
end

function CarttoSph(x::Array,y::Array,z::Array)
    r=sqrt(x.^2+y.^2+z.^2);
    θ=acos(z./r);
    ϕ=atan(y./x);
    return r,θ,ϕ;
end

"Defined Helper Functions"
{% endhighlight %}



Here, create a square cube, and convert those positions over to spherical coordinates.


{% highlight julia %}
range=-10:.5:10
x=collect(range);
y=collect(range);
z=collect(range);
N=length(x);
xa=repeat(x,outer=[1,N,N]);
ya=repeat(transpose(y),outer=[N,1,N]);
za=repeat(reshape(z,1,1,N),outer=[N,N,1]);
println("created x,y,z")

r,θ, ϕ=CarttoSph(xa,ya,za);
println("created r,θ,ϕ")
{% endhighlight %}


{% highlight julia %}
Ψ=Orbital(3,2,-1)
Ψp=Orbital(3,1,0)
{% endhighlight %}




{% highlight julia %}
Ψv = zeros(Float32,N,N,N);
ϕv = zeros(Float32,N,N,N);
{% endhighlight %}



{% highlight julia %}
for nn in 1:N
    for jj in 1:N
        for kk in 1:N
            val=Ψ(ρ(r[nn,jj,kk],2),θ[nn,jj,kk],ϕ[nn,jj,kk]);
            #val+=Ψp(ρ(r[nn,jj,kk],2),θ[nn,jj,kk],ϕ[nn,jj,kk]);
            Ψv[nn,jj,kk]=convert(Float32,abs(val));
            ϕv[nn,jj,kk]=convert(Float32,angle(val));
        end
    end
end

mid=round(Int,(N-1)/2+1);
Ψv[mid,mid,:]=Ψv[mid+1,mid+1,:]; # the one at the center diverges
Ψv=(Ψv-minimum(Ψv))/(maximum(Ψv)-minimum(Ψv) );
{% endhighlight %}


{% highlight julia %}
w,r = glscreen()

robj=visualize(Ψv)

#choose this one for surfaces of constant of intensity
view(visualize(robj[:intensities],:iso))

#choose this for a block of 3D density
#view(visualize(Ψv))
r()
{% endhighlight %}



## 2p Orbital

{% include image.html img="M4/Images/Orbitals2/2p_spaceb.png" title="2p" caption="2p Orbital block showing the density of the wavefunction."%}

{% include image.html img="M4/Images/Orbitals2/2p_surface.png" title="2p" caption="2p Orbital shown via isosurface."%}

## 3d orbitals
{% include image.html img="M4/Images/Orbitals2/3d0_surface.png" title="3d0" caption="3dz2 Orbital shown via isosurface. This corresponds to $n=3$, $l=2$, $m=0$." %}

{% include image.html img="M4/Images/Orbitals2/3d-1_surface.png" title="3dm1" caption="A 3d Orbital shown via isosurface. This corresponds to $n=3$, $l=2$, $m=-1$. This is not one of the canonical images, but instead an $m$ shape." %}

{% include image.html img="M4/Images/Orbitals2/3d2-2_spaceb.png" title="3dxy" caption="3dxy (x2-y2) orbital shown in density.  This is the sum of an $m=-2$ and $m=2$ state, for $n=3,l=2$. "%}

## 3p
In order to get this 3p surface image to come out correctly, I used the square root of the values instead in order to be able to see the much fainter outer lobe.

{% include image.html img="M4/Images/Orbitals2/3p_surface.png" title="3p" caption="3p surface plot." %}

{% include image.html img="M4/Images/Orbitals2/3p_spaceb.png" title="3p" caption="3p space plot." %}
