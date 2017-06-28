---
title: Monte Carlo Calculation of pi
layout: post
comments: True
category: Numerics
tags: [Monte Carlo]
author: Christina C. Lee
description: An Introduction to Monte Carlo methods by calculating Pi
---

## Monte Carlo- Random Numbers to Improve Calculations

When one hears "Monte Carlo", most people might think of something like this:

<center>
<img src="../images/MonteCarlo/871.jpg" width="40%" style="align:center; margin: 20px"/></center>

Monte Carlo, Monaco: known for extremely large amounts of money, car racing, no income taxes,and copious gambling.

In addition to Monaco, Europe, Las Vegas decided to host a Monte Carlo-themed casino as well.  So during the Manhattan project, when the best minds in the United States were camped out in the New Mexican desert, they had plenty of inspiration from Las Vegas, and plenty of difficult problems to work on in the form of quantifying the inner workings of nuclei.  Enrico Fermi first played with these ideas, but Stanislaw Ulam invented the modern Monte Carlo Markov Chain later.

At the same time, these scientists now had computers at their disposal.  John von Neumann programmed Ulam's algorithm onto ENIAC, Electronic Numerical Integrator and Computer, the very first electronic, general purpose computer, even though it did still run on vacuum tubes.

That still doesn't answer, why do random numbers actually help us solve problems?

Imagine you are visiting a new city for the first time, maybe Monte Carlo. You only have a day or two, and you want to really get to know the city.  You have two options for your visit

* Hit the tourist sites you researched online
* Wander around.  Try and communicate with the locals.  Find an out-of-the-way restaurant and sample food not tamed for foreigners.  Watch people interact.  Get lost.

Both are legitimate ways to see the city.  But depending on what you want, you might choose a different option.  The same goes for exploring physics problems.  Sometimes you want to go in and study just everything you knew about beforehand, but sometimes you need to wander around, not hit everything, but get a better feeling for what everything might be like.

## Buffon's Needle: Calculation of π
Even back in the 18th century, Georges-Louis Leclerc, Comte de Buffon posed a problem in geometric probability.  Nowdays, we use a slightly different version of that problem to calculate π and illustrate Monte Carlo simulations.

Suppose we have a square dartboard, and someone with really bad, completely random aim, even though he/she always at least hits inside the dartboard.  We then inscribe a circle inside that dartboard. After an infinity number of hits, what is the ratio between hits in the circle, and hits in the square?

{% include image.html img="M4/Images/MonteCarlo/dartboard.png" title="dartboard" caption="Randomly thrown darts than can either be in the circle or not." %}


$$
f= \frac{N_{circle}}{N_{square}} =\frac{\text{Area of circle}}{\text{Area of square}} =\frac{\pi r^2}{4 r^2}= \frac{\pi}{4}
$$

$$
    \pi = 4 f
$$

## Onto the Code!


{% highlight julia %}
using Plots
using Colors
violet=parse(Colorant,"#6c71c4");
cyan=parse(Colorant,"#3aa198");

plotlyjs()
{% endhighlight %}

We will generate our random numbers on the unit interval.  Thus the radius in our circumstance is $.5$.

 Write a function `incircle(r2)` such that if `r2` is in the circle, it returns true, else, it returns false.  We will use this with the julia function `filter`.  Assume `r2` is the radius squared, and already centered around the middle of the unit circle


 ```julia
 function incircle(r2)
     if r2<.25
         return true
     else
         return false
     end
 end
 ```


```julia
#The number of darts we will throw at the board.  We will see how accurate different numbers are
N=round.(Int,10.^(1:.25:4))
# We will perform each number multiple times in order to calculate error bars
M=15;
```


```julia
πapprox=zeros(Float64,length(N),M);

for ii in 1:length(N)
    #initialize an array of proper size
    X=zeros(Float64,N[ii],2);
    for jj in 1:M

        #popular our array with random numbers on the unit interval
        rand!(X);

        #calculate their radius squared
        R2=(X[:,1]-.5).^2+(X[:,2]-.5).^2

        # 4*number in circle / total number
        πapprox[ii,jj]=4.*length(filter(incircle,R2))/N[ii];

    end
end
```


{% highlight julia %}
# Get our averages and standard deviations
πave=sum(πapprox,2)/M;
πstd=std(πapprox,2);
{% endhighlight %}

## Analysis

So that was a nice, short little piece of code.  Lets plot it now to see means.


```julia
plot(xlabel="N points",ylabel="π estimate",title="Monte Carlo Estimate of π")
scatter!(repmat(N,M),reshape(πapprox,195),xscale=:log10,label="Estimates")
plot!(N,πave,yerr=πstd,xscale=:log10,label="Average of $M")
```

<iframe src="/M4/Images/MonteCarlo/pie_estimation.html"  style="border:none; background: #ffffff"
width="700px" height="450px"></iframe>



When we have fewer numbers of points, our estimates vary much more wildly, and much further from 3.1415926 .
But, at least, the guesses from our different runs all seem equally distributed around the correct value, so it seems we have no systematic error.

As we get up to $10^4$, our estimate starts getting much more accurate and consistent.


```julia
plot(N,πstd,label="Error",xscale=:log10)
plot!(xlabel="N points",ylabel="std",title="Dependence of error on number of points")
```
<iframe src="/M4/Images/MonteCarlo/pie_error.html"  style="border:none; background: #ffffff"
width="700px" height="450px"></iframe>


So what we guessed in the first plot about dispersion in estimate, we quantify here in this plot.  When we only have 10 darts, the guesses vary by up to .3, but when we get down to 1,000 trials, we are starting to be consistent to .0002.

Now lets just make a graphical representation of what we've been doing this whole time.  Plot our points on unit square, and color the ones inside a circle a different color.


```julia
X=zeros(Float64,1000);
Y=zeros(Float64,1000);
rand!(X);
rand!(Y);
inside =( (X-.5).^2+(Y-.5).^2 ) .< .25;
```


```julia
scatter(X[inside],Y[inside],label="Inside",
    markercolor=violet)
scatter!(X[.!inside],Y[.!inside],label="Outside",
    markercolor=cyan)
plot!(title="The Dartboard")
```

<iframe src="/M4/Images/MonteCarlo/dartboard.html"  style="border:none; background: #ffffff"
width="700px" height="450px"></iframe>


That's all folks!
Now here's a picture of some pie to congratulate you on calculating π.
{% include image.html img="M4/Images/MonteCarlo/pie.jpg" title="pie" caption="By Scott Bauer, USDA ARS - This image was released by the Agricultural Research Service, the research agency of the United States Department of Agriculture, with the ID K7252-47 (next).This tag does not indicate the copyright status of the attached work. A normal copyright tag is still required. See Commons:Licensing for more information.English | français | македонски | +/−, Public Domain, https://commons.wikimedia.org/w/index.php?curid=264106" %}
