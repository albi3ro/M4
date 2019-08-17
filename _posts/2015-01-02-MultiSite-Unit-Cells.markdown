---
title: Computationally Visualizing Crystals Pt. 2
layout: post
comments: True
category: General Audience
tags: [Lattices]
author: Christina C. Lee
description: Here we tile unit cells via Bravais lattice to create some common lattices like honeycomb and diamond.
image: "https://albi3ro.github.io/M4/Images/Projectile/pyrochloresvg2.png"
---

#### Christina C. Lee, github: albi3ro

<b>Prerequisites:</b> [Computationally Visualizing Crystals Pt. 1]({{base.url}}/M4/Crystal-Shapes.html)

Time to one-up the Bravais lattice from Part 1.  Many beautiful lattices don't adhere to the "every site the same" policy.  They still repeat, but just take a little bit longer to get around to doing so.

Take the Kagome Lattice below,


{% include image.html img="M4/Images/MultiUnit/kagomesvg.png" title="kagome" caption="A Kagome Lattice." %}

{% include image.html img="M4/Images/MultiUnit/kagomebasket.jpg" title="basket" caption="A basket woven in the Japanese kagome style. <sub>Wikimedia commons</sub>" %}

If we look at the stars at the center of triangles, we can recognize a point triangular Bravais lattice.  Now each of those stars stands for a grouping of three sites in a <i>unit cell</i>.  According to <a href="http://chemwiki.ucdavis.edu/Physical_Chemistry/Physical_Properties_of_Matter/Phases_of_Matter/Solids/Unit_Cell">Chem Wiki</a>, a unit cell is:

>A unit cell is the most basic and least volume consuming repeating structure of any solid. It is used to visually simplify the crystalline patterns solids arrange themselves in.

I chose these triangles to be be the unit cells above and in my computational representation below, but can you think of any other ways to represent the unit cell?

Turns out, there isn't a unique way.  We can go further and define the <i>Wigner-Seitz</i> unit cell, which uses the Bravais translations to pick out just ONE of the various possible definitions.

In my line of work though, we often use either the easiest to write down, or the one that has the symmetries we want.

### Introducing Some Lattice Options

You saw Kagome above.

The options I've put in now are:

* honeycomb
* kagome
* shuriken (a.k.a. square-Kagome)
* diamond
* pyrochlore


The ones implemented here, except for diamond, are frustrated lattices that I work with in my research. Honeycomb is well known in condensed matter physics for being the structure of graphene. This is an extremely important material right now, though I work with it in terms of the Kitaev spin model (to be discussed at a later date). Kagome and pyrochlore are also popular models within my community.  The shuriken lattice is more uncommon, but gaining ground in the frustration community.

{% include image.html img="M4/Images/MultiUnit/Shurikens.jpg" title="Shurikens" caption="Japanese Shurikens- a type of ninja fighting star. <sub>By kaex0r (http://www.flickr.com/photos/kaex0r414/191765028/) [CC BY 2.0 (http://creativecommons.org/licenses/by/2.0)], via Wikimedia Commons</sub>" %}

{% highlight julia %}
# importing our packages
using Plots;
plotlyjs()
{% endhighlight %}



{% highlight julia %}
lattice="shuriken";

Nx=3;
Ny=3;
Nz=1;
{% endhighlight %}


{% highlight julia %}
# A cell to just evaluate
# This one sets the unit vectors (a,b,c) for the different unit cells
# Can you guess what a lattice will look like by looking at the vectors?
if(lattice=="honeycomb") #also the graphite lattice
    d=2;
    Ncell=2;
    unit=[[0 0 0]
        [sqrt(3)/2 1/2 0]];
    a=[sqrt(3),0,0];
    b=[sqrt(3)/2,3/2,0];
    c=[0,0,1];
elseif(lattice=="kagome")
    d=2;
    Ncell=3;
    unit=[[0 0 0]
          [1 0 0]
        [.5 sqrt(3)/2 0]];
    a=[2,0,0];       #Look familiar? Checkout pt from Pt. 1
    b=[1, sqrt(3), 0];
    c=[0,0,1];
elseif(lattice=="shuriken")
    d=2;
    Ncell=6;
    unit=[[0 0 0]
          [.5 0 0]
          [0 .5 0]
          [.5 .5 0]
        [.5+.25*sqrt(3) .25 0]
        [.25 .5+.25*sqrt(3) 0]];
    a=[.5+.5*sqrt(3),0,0];
    b=[0,.5+.5*sqrt(3),0];
    c=[0,0,1];
elseif(lattice=="diamond")
    d=3;
    Ncell=2;
    unit=[[0 0 0]
          [.25 .25 .25]];
    a=[.5,.5,0];    #Look familiar? Checkout fcc from Pt.1
    b=[.5,0,.5];
    c=[0,.5,.5];
elseif(lattice=="pyrochlore")
    d=3;
    Ncell=4;
    unit=[[0 0 0]
        [.25 .25 0]
        [.25 0 .25]
        [0 .25 .25]];
    a=[.5,.5,0];
    b=[.5,0,.5];
    c=[0,.5,.5];

else
    println("Please have a correct lattice")
end
"Cell finished"
{% endhighlight %}



### Connections to Bravais Lattices

If you look at some of the comments above, and checkout the basis vectors from [<i>Crystal Shapes</i>]({{base.url}}/M4/General/Crystal-Shapes.html), like pt,
        \begin{equation}
            a=[1,0,0]\;\;\;\;\;\;\;\;\; b=[.5,\frac{\sqrt{3}}{2},0],
        \end{equation}
you'll notice they're the same except for a scaling factor.  This has to be true, since only 14 different patterns tile 3D space uniquely.


{% highlight julia %}
# Another cell to just evaluate
# Here we set up some numbers and matrices for our computation
N=Nx*Ny*Nz*Ncell;    #The total number of sites
aM=transpose(repeat(a,outer=[1,Ncell]));
bM=transpose(repeat(b,outer=[1,Ncell*Nx])); #these allow us to copy an entire row or layer at once
cM=transpose(repeat(c,outer=[1,Ncell*Nx*Ny]));

X=Array{Float64}(undef,N,3);  #where we store the positions
"Cell finished"
{% endhighlight %}



```julia
# Another cell to just evaluate
# Here we are actually calculating the positions for every site
for i in 1:Nx    #for the first row
    X[Ncell*i-Ncell+1:Ncell*i,:]=unit.+(i-1)*aM;
end

for j in 2:Ny    #copying the first row into the first layer
    X[Ncell*Nx*(j-1).+(1:Ncell*Nx),:]=X[1:Ncell*Nx,:].+(j-1)*bM;
end

for j in 2:Nz    #copying the first layer into the entire cube
    X[Ncell*Ny*Nx*(j-1).+(1:Ncell*Nx*Ny),:]=X[1:Ncell*Nx*Ny,:].+(j-1)*cM;
end
"Cell finished"
```




{% highlight julia %}
scatter(X[:,1],X[:,2])
# scatter(X[:,1],X[:,2],X[:,3])
plot!(title=lattice,xlabel="X",ylabel="Y",zlabel="Z",legend=false)
{% endhighlight %}


<iframe src="/M4/Images/MultiUnit/honeycomb.html"  style="border:none; background: #ffffff"
width="700px" height="450px"></iframe>

<iframe src="/M4/Images/MultiUnit/kagome.html"  style="border:none; background: #ffffff"
width="700px" height="450px"></iframe>

<iframe src="/M4/Images/MultiUnit/shuriken.html"  style="border:none; background: #ffffff"
width="700px" height="450px"></iframe>

<iframe src="/M4/Images/MultiUnit/diamond.html"  style="border:none; background: #ffffff"
width="700px" height="450px"></iframe>

<iframe src="/M4/Images/MultiUnit/pyrochlore.html"  style="border:none; background: #ffffff"
width="700px" height="450px"></iframe>


### Perdy Pictures
From these plots, some 3D structures like the pyrochlore are hard to visualize.  So here's a nice graphic I made that might help a little bit more.

{% include image.html img="M4/Images/MultiUnit/pyrochloresvg2.png" title="Pyrochlore" caption="Hopefully this pyrochlore is a little easier to visualize than the pyplot version.  Took me long enough to make in inkscape." %}

{% include image.html img="M4/Images/MultiUnit/honeycomb.png" title="honeycomb" caption="Tikz produced Honeycomb.  Coloring indicative of the lattice description of the Kitaev model.  " %}

The honeycomb, like several other lattices you see around here, is <i>bipartite</i>.
You can see in my image that black sites are only next to white sites, and vice versa.
  This property can make the system much easier to work with.  What lattices are bipartite, and which ones aren't?

If you keep reading, these lattices will keep cropping up again and again.  I'll probably throw in some new ones as well.

Anyway, we will move onto some Quantum Mechanics to look at atomic orbitals soon!
