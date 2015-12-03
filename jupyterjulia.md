---
layout: page
title: Interacting with Julia
permalink: /jupyterjulia/
---

## Julia - A possible future for Scientific Computing

The blog will be written in the form of julia (<a href="http://julialang.org/">http://julialang.org/</a>) jupyter notebooks. This new language has a syntax similar to Matlab or Python, and is easily human readable.  I encourage active participation by manipulating the code.

You can either install Julia and jupyter notebook on your own machine, or use juliabox (<a href="https://www.juliabox.org/">https://www.juliabox.org/</a>) to view and interact with the files.

To use juliabox (the easiest option), log in with a google account, and then click on the sync tab.  Under the Github directory, input the clone <a href=https://github.com/albi3ro/M4.git">https://github.com/albi3ro/M4.git</a> , use branch master, and name the folder as you see fit.

Check out <a href="http://julialang.org/">http://julialang.org/</a> for documentation on this new and exciting language, and learn how to install it on your own machine.

### Introductory Note on Packages
Many times I will use extra packages developed for Julia, such as plotting functionality.  If you are using the package for the first time on your machine, you will need to evaluate `Pkg.add("...")` and `Pkg.update()` before `using ...`

If you are using JuliaBox, you will need to evaluate the those lines each time.

The first time a packaged is used on a machine by `using ...`, the computer has to compile the package, and that can take a little bit of time.  As the language develops, the developers intend to create a more efficient system, but for now, we need patience.
