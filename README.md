# M4
## Metals, Magnets, and Miscelaneous Materials: Blog Introduction to Comptational Condensed Matter

### by Christina Lee,
### PhD student at the Okinawa Institute of Science and Technology

### Checkout http://albi3ro.github.io/M4/ to get the static blog website version of these jupyter notebooks.


You thought particle physics and astronomy were the only big fronteirs in physics? Well, you thought wrong.
Follow this blog to start discovering some of the wonders of materials all around you, and some of the revolutions that are shocking the field.

I will vary the level from undergraduate physics major to graduate student in the field.  I still encourage the curious layperson to peruse.  So you know what you are getting yourself into, I've developed a difficulty ranking scheme:

*<b>General:</b> Those with some background in programming, math, and physics should be able to follow.
*<b>Prerequisites:</b> I will spell out classes in the beginning that one should have taken
    in order to understand the physics of the post.  For example, in order to understand calculating
    the orbitals of a hydrogen-like atom, I will assume you analytically solved it in a Quantum Mechanics course.
    Even if you haven't taken  the classes, feel free to still try.
*<b>Research:</b> A general research technique.
*<b>In Field:</b> A specialized technique.  Unless you are doing research on the topic, or need to use this tool,
     you will probably not learn this technique, though researchers might hear about it. Consider where you want to spend your time.
*<b>Numerics:</b>Without a specific physics goal.  Sometimes I may just want to talk about something programming or
    numerical method related.
*<b>Fun:</b> Something off topic and silly for everyone!


## Julia - A possible future for Scientific Computing

The blog will be written in the form of julia (http://julialang.org/) jupyter notebooks. This new language has a syntax similiar to Matlab or Python, and is easily human readable.  I encourage active participation by manipulating the code.

You can either install Julia and jupyter notebook on your own machine, or use juliabox (https://www.juliabox.org/) to view and interact with the files.

To use juliabox (the easiest option), log in with a google account, and then click on the sync tab.  Under the Github directory, input the clone https://github.com/albi3ro/M4.git , use branch master, and name the folder as you see fit.

Check out http://julialang.org/ for documentation on this new and exciting language, and learn how to install it on your own machine.

### Introductory Note on Packages
Many times I will use extra packages developed for Julia, such as plotting functionality.  If you are using the package for the first time on your machine, you will need to evaluate `Pkg.add("...")` and `Pkg.update()` before `using ...`

If you are using JuliaBox, you will need to evaluate the those lines each time.

The first time a packaged is used on a machine by `using ...`, the computer has to compile the package, and that can take a little bit of time.  As the language develops, the developers intend to create a more efficient system, but for now, we need patience.

### Using My Packages

I have written packages of my own for use on certain sections.  To add these to the path, download them to the folder you are using, and in your command line, type `push!(LOAD_PATH,"/path/to/my/code")`, where `/path/to/my/code` is the location of the package.  Then, you can type `using Lattices` in just the same way.

## For further reference
Interested in learning more?

* Principles of Condensed Matter Physics- Chaikin & Lubensky
* Condensed Matter in a Nutshell- Mahan
* Inna Vishik on Quora https://www.quora.com/profile/Inna-Vishik
* https://thiscondensedlife.wordpress.com/

Know any good ones? Let me know!
