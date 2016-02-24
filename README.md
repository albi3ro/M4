# M4
## Metals, Magnets, and Miscelaneous Materials: Blog Introduction to Computational Condensed Matter

### by Christina Lee,
### PhD student at Okinawa Institute of Science and Technology Graduate University, Japan.

### Check out http://albi3ro.github.io/M4/ to get the static blog website version of these jupyter notebooks.


You thought particle physics and astronomy were the only big frontiers in physics? Well, you thought wrong.
Follow this blog to start discovering some of the wonders of materials all around you, and some of the revolutions that are shocking the field.

I will vary the level of the content from undergraduate physics major to graduate student in the field.  However, I still encourage the curious layperson to follow along.  Later, I will develop a difficulty ranking scheme.

## Julia - A possible future for scientific computing

The blog will be written in the form of Julia (http://julialang.org/) Jupyter notebooks. This new language has a syntax similar to MATLAB or Python, and is easily human readable.  I encourage active participation by playing with the code.

You can either install Julia and Jupyter notebook on your own machine, or use JuliaBox (https://www.juliabox.org/) to view and interact with the files.

To use JuliaBox (the easiest option), log in with a Google account, and then click on the sync tab.  Under the Github directory, input the clone https://github.com/albi3ro/M4.git , use branch master, and name the folder as you see fit.

Check out http://julialang.org/ for documentation on this new and exciting language, and learn how to install it on your own machine.

### Introductory note on packages
Many times I will use extra packages developed for Julia, such as for plotting functionality. If you are using the package for the first time on your machine, you will need to evaluate `Pkg.add("...")` and `Pkg.update()` before `using ...`

If you are using JuliaBox, you will need to evaluate the those lines each time.

The first time a packaged is used on a machine by `using ...`, the computer has to compile the package, and that can take a little bit of time.  As the language develops, the developers intend to create a more efficient system, but for now, we need patience.

### Using my packages
I have written packages of my own for use in certain sections.  To add these to the path, download them to the folder you are using, and in your command line, type `push!(LOAD_PATH,"/path/to/my/code")`, where `/path/to/my/code` is the location of the package.  Then, you can type `using Lattices` in just the same way.    

## For further reference
Interested in learning more?

* Principles of Condensed Matter Physics - Chaikin & Lubensky
* Condensed Matter in a Nutshell - Mahan
* Inna Vishik on Quora https://www.quora.com/profile/Inna-Vishik
* https://thiscondensedlife.wordpress.com/

Know any good ones? Let me know!
