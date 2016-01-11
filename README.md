# M4
## Metals, Magnets, and Miscelaneous Materials: Blog Introduction to Comptational Condensed Matter

### by Christina Lee,
### PhD student at the Okinawa Institute of Science and Technology

### Checkout http://albi3ro.github.io/M4/ to get the static blog website version of these jupyter notebooks.


You thought particle physics and astronomy were the only big fronteirs in physics? Well, you thought wrong.
Follow this blog to start discovering some of the wonders of materials all around you, and some of the revolutions that are shocking the field.

I will vary the level from undergraduate physics major to graduate student in the field.  I still encourage the curious layperson to give some a try.  Later, I will develop a difficulty ranking scheme.

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


<sub><center>
The MIT License (MIT)

<sub><sub>
Copyright (c) 2015 Christina Colleen Lee

<sub><sub>
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

<<<<<<< HEAD
<sub><sub>
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
=======
<ul>
<li>Principles of Condensed Matter Physics- Chaikin & Lubensky
<li>Condensed Matter in a Nutshell- Mahan
<li>Inna Vishik on Quora https://www.quora.com/profile/Inna-Vishik
<li>https://thiscondensedlife.wordpress.com/
>>>>>>> af75f7d02cbf586c0e3b54b3ce3f1f895e6f839f

<sub><sub>
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
