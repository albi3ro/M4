---
layout: page
title: Interacting with Julia
permalink: /jupyterjulia/
---

## Julia - A possible future for Scientific Computing

The blog will be written in the form of julia (<a href="http://julialang.org/">http://julialang.org/</a>) jupyter notebooks. This new language has a syntax similar to Matlab or Python, and is easily human readable.  I encourage active participation by manipulating the code.

Check out <a href="http://julialang.org/">http://julialang.org/</a> for documentation on this new and exciting language, and learn how to install it on your own machine.  "Docs" on the right hand side provides an easy link to the documentation.  

You can either install Julia and jupyter notebook on your own machine, or use juliabox (<a href="https://www.juliabox.com/">https://www.juliabox.com/</a>) to view and interact with the files.

To use juliabox (the easiest option), log in with a google account, and then click on the sync tab.  Under the Github directory, input the clone <a href="https://github.com/albi3ro/M4.git">https://github.com/albi3ro/M4.git</a> , use branch master, and name the folder as you see fit.

### Introductory Note on Packages
Many times I will use extra packages developed for Julia, such as plotting functionality.  If you are using the package for the first time on your machine, you will need to evaluate `Pkg.add("...")` and `Pkg.update()` before `using ...`

If you are using JuliaBox, you will need to evaluate the those lines each time.

The first time a packaged is used on a machine by `using ...`, the computer has to compile the package, and that can take a little bit of time.  As the language develops, the developers intend to create a more efficient system, but for now, we need patience.

### My Packages

In some blogs, I use modules that I have personally written in order to reduce the amount of unnesseary code present in the notebook.  When one of these is present:

* Make sure the corresponding *.jl file is downloaded
* Make sure the *.jl file is on your path

You add the location of the file each time by typing `push!(LOAD_PATH,"\path\to\code")`.

You can permanently add a location to julia by adding the line `@everywhere push!(LOAD_PATH,"\path\to\code")` to the file `.juliarc.jl` in your home directory.
