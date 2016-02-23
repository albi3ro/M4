---
title: Computationally Visualizing Crystals
layout: post
comments: True
category: General
tags: [Lattices]
author: Christina C. Lee
---

#### Christina C. Lee, github: albi3ro

<b>Prerequisites:</b> none


In condensed matter, we find ourselves in the interesting middle ground of dealing with large numbers (e.g. $10^{23}$) of extremely small particles such as atoms, or electrons.

Luckily, the particles don't each do their own thing, but often come in nice, structured, repeated units.  <i>Lattices</i>.  As our first step into the field, we will look at the most basic type, a <i>Bravais lattice</i>.

In a Bravais lattice, every site looks like every other site. Mathematically, we use three vectors, $\vec{a},\vec{b},\vec{c}$ to express how we move from one site to a neighbor.

\begin{equation}
\mathbf{R}_{lmn}=l \vec{a} + m \vec{b} + n \vec{c}  \;\;\;\; \text{for } l,m,n \in \mathbb{N}
\end{equation}

For consistency, we have to put a constraint on these vectors; we cannot combine two of the vectors and obtain the third.  If we could, then we couldn't have sites in an entire three dimensional space.

Stay tuned for a later post where we explore more elaborate lattices.


```julia
# importing our packages
Pkg.add("PyPlot");
Pkg.update();
using PyPlot;
```


## Define The Relevant Variables

Choose the lattice you want to look at, and use that string for the lattice variable.
Current options:

* Simple Cubic = "sc"
* Plane triangular lattice = "pt"
* Body-Centered Cubic = "bcc"
* Face-Centered Cubic = "fcc"

Note: Square is Simple Cubic for Nz=1

14 distinct lattice types are possible, but these common four give the important ideas.

Also, input the size of lattice you want to look at.


```julia
lattice="sc";

Nx=3;
Ny=3;
Nz=3;
```


```julia
# A cell to just evaluate
# This one sets the unit vectors (a,b,c) for the different unit cells
# Can you guess what a lattice will look like by looking at the vectors?
if(lattice=="sc")
    d=3;
    a=[1,0,0];
    b=[0,1,0];
    c=[0,0,1];
elseif(lattice=="pt")
    d=2;
    a=[1,0,0];
    b=[.5,sqrt(3)/2,0];
    c=[0,0,1];
elseif(lattice=="bcc")
    d=3;
    a=[.5,.5,.5];
    b=[.5,.5,-.5];
    c=[.5,-.5,.5];
elseif(lattice=="fcc")
    d=3;
    a=[.5,.5,0];
    b=[.5,0,.5];
    c=[0,.5,.5];
else
    println("Please have a correct lattice")
end
```




```julia
# Another cell to just evaluate
# Here we set up some numbers and matrices for our computation
N=Nx*Ny*Nz;    #The total number of sites
aM=transpose(a);
bM=transpose(repeat(b,outer=[1,Nx])); #these allow us to copy an entire row or layer at once
cM=transpose(repeat(c,outer=[1,Nx*Ny]));

X=Array{Float64}(N,3);  #where we store the positions
```





```julia
# Another cell to just evaluate
# Here we are actually calculating the positions for every site
for i in 1:Nx    #for the first row
    X[i,:]=(i-1)*a;
end

for j in 2:Ny    #copying the first row into the first layer
    X[Nx*(j-1)+(1:Nx),:]=X[1:Nx,:]+(j-1)*bM;
end

for j in 2:Nz    #copying the first layer into the entire cube
    X[Ny*Nx*(j-1)+(1:Nx*Ny),:]=X[1:Nx*Ny,:]+(j-1)*cM;
end
```
<div class="progtip">
<h3 color="black"> Programming Tip:</h3>
 <p>In Julia, ranges, like <code>1:Nx</code>, are a special variable type that can be manipulated.  We can add numbers to them:
 <code>3+(1:3)=4:6</code>,
 or add a minus sign to force it to iterate in the opposite direction, though with different start/stop:
 <code>-(1:3)=-1:-1:-3</code></p>
<p>
 <span color="#000000">Danger!</span> Make sure to use the parentheses around the range if you are performing these operations.</p>
</div>


```julia
pygui(false);  #if true, launches new window with interactive capabilities

drawcube=true;  #gives lines for a cube, helps interpret the dots
ls=2;  # how many cubes to draw
if(drawcube==true)
    v=collect(0:ls);
    zed=zeros(v);
    for i in 0:ls
        for j in 0:ls
            plot3D(zed+i,v,zed+j)
            plot3D(zed+i,zed+j,v)

            plot3D(v,zed+i,zed+j)
            plot3D(zed+j,zed+i,v)

            plot3D(v,zed+j,zed+i)
            plot3D(zed+j,v,zed+i)
        end
    end
end

scatter3D(X[:,1],X[:,2],X[:,3],s=200*ones(X[:,1]),alpha=1)
```

{% include image.html img="M4/Images/CrystalShapes/sc.png" title="sc" caption="Simple Cubic: The easiest lattice out there short of the 1D chain." %}

{% include image.html img="M4/Images/CrystalShapes/pt.png" title="pt" caption="Point Triangular: A 2D lattice.  Plotted using scatter instead of scatter3D." %}

{% include image.html img="M4/Images/CrystalShapes/bcc.jpg" title="bcc" caption="Body Centered Cubic:  Notice how some sites fall on the cubic lattice, but others fall in between.  Generated with pygui(true) and then manipulating in 3D." %}

{% include image.html img="M4/Images/CrystalShapes/fcc.jpg" title="fcc" caption="Face Centered Cubic: Here the sites either fall on the the cubic corners of in the center of the sides.   Generated with pygui(true), ls=1, and then manipulating in 3D." %}

## Go Back and Fiddle!

As you might have noticed, this isn't just a blog where you read through the posts.  Interact with it.  Change some lines, and see what happens.  I choose body centered cubic to display first, but what do the other lattices look like?

Chose `pygui(true)` to pop open a window and manipulate the plot in 3D.

Look at different lattice sizes.

Can you hand draw them on paper?

{% include image.html img="M4/Images/CrystalShapes/handdraw2.jpg" title="handdrawn fcc" caption="A face centered cubic I decided to draw myself. " %}


Let me know what you think, and enjoy the sequel as well!

 [<center><i> Multi-site unit cells</i></center>]({{base.url}}/M4/general/MultiSite-Unit-Cells.html)
