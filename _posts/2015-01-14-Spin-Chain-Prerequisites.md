---
title: 1D Spin Chain Prerequisites
layout: post
comments: True
category: Graduate
author: Christina C. Lee
image: "https://albi3ro.github.io/M4/Images/SpinChainED/mag.png"
tags:
  Magnet
  Quantum
  ED
series: 1D_Spin_Chain
description: The first in a series on Exact Diagonalization of the Quantum Mechanical Heisenberg Spin Chain Model, this post covers theoretical background and how we will be manipulating binary numbers.
---

Today's just background :(  I found my post getting a bit bloated, so I decided to pull most of the talking into an introduction post for the series.   Keep reading to learn how to compute, diagonalize, and analyze the matrix corresponding to a many-body, strongly interacting quantum mechanical system.  

If you are not a physics graduate student, still give this series a read and try not to get too intimidated. If you are a physics graduate student, still give this series a read and try not to get too intimidated.  I've been working on this series for quite a while, so I hope it's worth my effort.

For further information, take a look at this set of lecture notes on the arXiv, [Computational Studies of Quantum Spin Systems](https://arxiv.org/abs/1101.3281) by Anders W. Sandvik.


_______________________
## So what is a 1D Heisenberg Spin $\frac{1}{2}$ Chain?
Let's break down each part of that phrase.  

#### Spin $\frac{1}{2} $
At each site, we can have a particle pointing up $| \uparrow \rangle$, down $|\downarrow \rangle$, or some super-position of the two.  

### Heisenberg
Our spin has three degrees of freedom and full $SU(2)$ symmetry.  $SU(2)$ is the mathematical group that describes a spin's degrees of freedom. Once we have solved the physics of the Heisenberg case, we also have solved the <i>XY model</i> ($J_z=0$) and the <i>Ising model</i> $J_x=J_y=0$. As we change the model, we also change the symmetry group first to $SO(2)$ for the XY model and then $\mathbb{Z}_2$ for the Ising model.

### 1D Chain
Spin only interacts with two neighbors.

{% include image.html img="M4/Images/SpinChainED/spinchain.svg" title="Spin Chain" caption="A spin chain." %}

## The Hamiltonians

Our Hamiltonian has the general form,

$$
{\cal H} = \sum_i J_x S_i^x S_{i+1}^x + J_y S_i^y S_{i+1}^y + J_z S_i^z S_{i+1}^z
$$

or we can restrict it to the case of $J_x=J_y$ to get the more convenient form,

$$
{\cal H} = \sum_i J_{XY} \left(S^+_i S^-_{i+1}+S^-_i S^+_{i+1} \right) +J_z S_i^z S_{i+1}^z.
$$

Here $S^{x,y,z}$ are our Pauli operators

$$
S^x=\begin{pmatrix}
0 & 1 \\
1 & 0 \\
\end{pmatrix}
\;\;\;\;\;\;\;
S^y=\begin{pmatrix}
0    &    -i    \\
i    &    0    \\
\end{pmatrix}
\;\;\;\;\;\;\;
S^z=\begin{pmatrix}
1    &    0    \\
0    &    -1    \\
\end{pmatrix},
$$

and $S^{\pm}$ are the ladder operators

$$
S^+=\frac{S^x+i S^y}{2} = \begin{pmatrix}
0 & 1 \\
0 & 0\\
\end{pmatrix}
\;\;\;\;\;\;\;
S^- =\frac{S^x-i S^y}{2}= \begin{pmatrix}
0 & 0 \\
1 & 0\\
\end{pmatrix}.
$$


Assuming we write our basis states in the $S^z_i$ basis, we can divide the terms from the restricted Hamiltonian into on-diagonal and off-diagonal terms.  The $S^z_i S^z_{i+1}$ terms compute the magnetization squared, $\vec{S} \cdot \vec{S} $, for a given state and a conserved quantity.  These also lie on the diagonal of the matrix corresponding to the Hamiltonian.

$$
| \Psi_{i} \rangle = H_{ii} |\Psi_i \rangle
= J_{z} \sum_m S_m^z S^z_{m+1}  |\Psi_i \rangle
$$

The ladder terms, when applied as an operator to the state, create a new state.  Thus they act as off-diagonal terms.

$$
| \Psi_{j} \rangle = H_{ij} |\Psi_i \rangle
= J_{XY} \sum_m \left( S^+_m S^-_{m+1} + S^-_m S^+_{m+1}\right) |\Psi_i \rangle
$$


## Analytical Solutions

    "The model is solvable, not the solution."
             --- someone


The Ising, XY, and Heisenberg cases all fall into the special class of problems which have exact solutions in the infinite size limit.  Very few exact solutions of quantum mechanics problems exist, so we try and get as much mileage as we can out of the couple of ones we have.  

Interestingly enough, once one solution to a problem comes along, someone figures out a different way to approach the same problem.

The 1D Quantum Ising Model is equivalent to the [2D Ising model of classical statistical mechanics](http://albi3ro.github.io/M4/Monte-Carlo-Ferromagnet.html), exactly solved in 1944 by Lars Onsager.  The solution is also equivalent to a description of free Majorana fermions. [1]

The Jordan-Wigner Transformation solves the 1D XY model by mapping spins to fermions.  This transformation only works in special 1D circumstances and the Kitaev Honeycomb model.  Since spins possess different anti-commutation relationships than fermions, we attach a string of operators stretching from infinity to each spin.  This series of operators changes the relationship between a spin and its neighbors to fermionic. After the transformation, we get a Hamiltonian that is quadratic in the fermionic momentum operators $d_k $ , $d^{\dagger}_k$, and we can see the $\cos (ka)$ dispersion relationship for the excitations,

$$
{\cal \tilde{H}}_{XY}=-J \sum_k \cos (ka) d^{\dagger}_k d_k.
$$

I might write a full article on this later.  

Performing a Jordan-Wigner transformation on the full Heisenberg model gives a four-operator scattering term.  The Bethe Ansatz, which I honestly don't know anything about, solves the full 1D Heisenberg Model, as well as some 1D Bose gas and Hubbard model problems.  Come back to me in many years, or ask a Russian mathematician.  

## Noether's Theorem: Symmetries and Conserved Quantities

While many people look to Marie Curie, I think Emily Noether is the best female scientist/ mathematician yet.  Her theorem is the single most beautiful piece of physics I have ever seen.  For each symmetry, there exists a conserved quantity.  Space translational symmetry gives us momentum conservation.  Time translational symmetry gives us energy conservation.  

Conserved quantities save us so much in classical mechanics but save us even more in condensed matter physics.  By working out the symmetries of the problem, we can break the huge Hilbert space into smaller chunks of conserved quantities and only have a smaller problem with which to work.  

 So what are our conserved quantities?
* Magnetization from 3*N Rotational symmetries
* Momentum from N Translational symmetries
* Parity from Time inversion symmetry

I will only demonstrate magnetization conservation since we won't have to alter our basis to accommodate it.


### One easy initialization cell.

`n`: The number of spins we will deal with.  
$2^n$: Dimension of our Hilbert Space.  Or more simply, the number of different eigenstates our system will have.

Try to keep `n` around `5` until you know what your computer can handle.  


```julia
n=4
nstates=2^n
```


Exact Diagonalization is often memory limited.  Thus we want to represent our states in the most compact format possible.  Luckily, if we are dealing with spin $\frac{1}{2}$, we can just use the `0`'s ,$\mid \downarrow \rangle $ , and `1`'s, $\mid \uparrow \rangle $,  of the machine.  If you are dealing with higher spin, you can use base 3, 4, etc...  Part of the reason I needed to create this separate post was to examine working with binary data.

We will keep our states stored as Int, but Julia has operations we can perform to look at the binary format and change the bits.


```julia
# psi is an array of all our wavefunctions
psi=convert.(Int8,collect(0:(nstates-1)) )

# Lets look at each state both in binary and base 10
println("binary form \t integer")
for p in psi
    println(bitstring(p)[end-n:end],"\t\t ",p)
end
```

      binary form 	 integer
      00000		 0
      00001		 1
      00010		 2
      00011		 3
      00100		 4
      00101		 5
      00110		 6
      00111		 7
      01000		 8
      01001		 9
      01010		 10
      01011		 11
      01100		 12
      01101		 13
      01110		 14
      01111		 15



Because we will be using the powers of $2$ frequently in our calculations, we will store all them in an array. They are our placeholders, like $1,10,100,...$


```julia
powers2=convert.(Int8, 2.0 .^collect(0:(n-1)) )

println("binary form \t integer")
for p in powers2
    println(bitstring(p)[end-n:end],"\t\t ",p)
end
```

      binary form 	 integer
      00001		 1
      00010		 2
      00100		 4
      01000		 8



### & And Operation and Computing Magnetization

Once we have magnetization for a state, a conserved quantity, we also have magnetization squared for the diagonals.

We could continue to look at the binary format of a number by calling `bin`, but that converts the number to an array of strings.  So instead we want to perform bitwise operations to determine what the binary format looks like in terms of numbers.

Julia supports bitwise <b>not, and, xor </b>/ exclusive or, logical shift right, arithmetic shift right, and logical/ arithmetic shift left.  For our purposes, we will only be interested in <b>and</b> and <b>xor</b> .

<b>and</b> takes in two inputs and produces one output, given by the following logic table:

<center>
<table>
  <tr>
    <th> a</th>
    <th> b</th>
    <th> a&b</th>
  </tr>
  <tr>
    <td> 0 </td>
    <td> 0 </td>
    <td> 0</td>
  </tr>
  <tr>
    <td> 1 </td>
    <td> 0 </td>
    <td> 0</td>
  </tr>
  <tr>
    <td> 0 </td>
    <td> 1 </td>
    <td> 0</td>
  </tr>
  <tr>
    <td> 1 </td>
    <td> 1 </td>
    <td> 1</td>
  </tr>
</table>
</center>
Julia's `&` is the bitwise operation and.  That means if I combine two numbers, it states the overlap between the two. 1 overlaps with 1; 2 overlaps with 2; 3 overlaps with 2 and 1.

We will use this to compute magnetization.


```julia
println(bitstring(1)[end-2:end],"\t", 
    bitstring(3)[end-2:end],"\t", bitstring(1&3)[end-2:end])
```

    001	011	001



```julia
#initializing the magnetization array
m=zeros(length(psi))

println("String \tp&powers2 \tNormed \t\tMagentization")
for i in 1:length(psi)
    
    #Writing the magnetization
    m[i]=sum(round.(Int,(psi[i].&powers2)./powers2))
    
    println(bitstring(psi[i])[end-n:end],"\t",psi[i].&powers2,"\t"
    ,round.(Int,(psi[i].&powers2)./powers2),"\t",m[i]) 
end
```

      String 	p&powers2 	Normed 		Magentization
      00000	Int8[0, 0, 0, 0]	[0, 0, 0, 0]	0.0
      00001	Int8[1, 0, 0, 0]	[1, 0, 0, 0]	1.0
      00010	Int8[0, 2, 0, 0]	[0, 1, 0, 0]	1.0
      00011	Int8[1, 2, 0, 0]	[1, 1, 0, 0]	2.0
      00100	Int8[0, 0, 4, 0]	[0, 0, 1, 0]	1.0
      00101	Int8[1, 0, 4, 0]	[1, 0, 1, 0]	2.0
      00110	Int8[0, 2, 4, 0]	[0, 1, 1, 0]	2.0
      00111	Int8[1, 2, 4, 0]	[1, 1, 1, 0]	3.0
      01000	Int8[0, 0, 0, 8]	[0, 0, 0, 1]	1.0
      01001	Int8[1, 0, 0, 8]	[1, 0, 0, 1]	2.0
      01010	Int8[0, 2, 0, 8]	[0, 1, 0, 1]	2.0
      01011	Int8[1, 2, 0, 8]	[1, 1, 0, 1]	3.0
      01100	Int8[0, 0, 4, 8]	[0, 0, 1, 1]	2.0
      01101	Int8[1, 0, 4, 8]	[1, 0, 1, 1]	3.0
      01110	Int8[0, 2, 4, 8]	[0, 1, 1, 1]	3.0
      01111	Int8[1, 2, 4, 8]	[1, 1, 1, 1]	4.0


## Masks and Permuting Indices

The off diagonal, ladder, elements of the Hamiltonian are the permutation of two neighboring elements in the array.  We can permute two indices by combining a mask number with a bitwise XOR `$` .


```julia
mask=Int8(3)
testp=1
println("Mask \ttest \tmasked test")
println(bitstring(mask)[end-2:end],'\t',bitstring(testp)[end-2:end]
    ,'\t',bitstring(testp⊻mask)[end-2:end])
```

    Mask 	test 	masked test
    011	001	010


But the mask 3 aka 11 only switches the spins in the first two positions.  I need to switch spins in any two adjacent locations.  I create this by summing together padded powers of two in order to get the 11 in the correct location.


```julia
mask=convert.(Int8,[0;powers2]+[powers2;0])
mask=mask[2:end-1]

println("Mask base10 \tMask Binary \tSummed from")
for i in 1:length(mask)
    println(mask[i],"\t\t",bitstring(mask[i])[end-n:end],"\t\t",
        bitstring(powers2[i])[end-n:end],"\t",
        bitstring(powers2[i+1])[end-n:end])
end
```

    Mask base10 	Mask Binary 	Summed from
    3		00011		00001	00010
    6		00110		00010	00100
    12		01100		00100	01000



So now lets test how the first of our three masks behaves:
We know that if the mask changes a 01 for a 10, or vice versa, that the overall magnetization will not be changed.  So, we test is our mask is successful by comparing the remaining magnetization.  The rows offset by two spaces have matching magnetizations.


```julia
println("Psi \tPsi \tMasked \t Masked\t \tmPsi  \tmMasked")
for p in psi
    if m[p+1]==m[p.⊻mask[1]+1]
        println("  ",p,"\t  ",bitstring(p)[end-n:end],"\t  ",
            p.⊻mask[1],"\t  ",bitstring(p.⊻mask[1])[end-n:end],
            "\t\t  ",m[p+1],"\t  ",m[p.⊻mask[1]+1]) 
    else
        println(p,'\t',bitstring(p)[end-n:end],'\t',
            p.⊻mask[1],'\t',bitstring(p.⊻mask[1])[end-n:end],
            "\t\t",m[p+1],"\t",m[p.⊻mask[1]+1]) 
    end
end
```

      Psi 	Psi 	Masked 	 Masked	 	mPsi  	mMasked
      0	00000	3	00011		0.0	2.0
        1	  00001	  2	  00010		  1.0	  1.0
        2	  00010	  1	  00001		  1.0	  1.0
      3	00011	0	00000		2.0	0.0
      4	00100	7	00111		1.0	3.0
        5	  00101	  6	  00110		  2.0	  2.0
        6	  00110	  5	  00101		  2.0	  2.0
      7	00111	4	00100		3.0	1.0
      8	01000	11	01011		1.0	3.0
        9	  01001	  10	  01010		  2.0	  2.0
        10	  01010	  9	  01001		  2.0	  2.0
      11	01011	8	01000		3.0	1.0
      12	01100	15	01111		2.0	4.0
        13	  01101	  14	  01110		  3.0	  3.0
        14	  01110	  13	  01101		  3.0	  3.0
      15	01111	12	01100		4.0	2.0

Now let's try the same thing, but for the second mask.  This changes the second and third spin.


```julia
println("Psi \tPsi \tMasked \t Masked\t \tmPsi  \tmMasked")
for p in psi
    if m[p+1]==m[p.⊻mask[2]+1]
        println("  ",p,"\t  ",bitstring(p)[end-n:end],"\t  ",
            p .⊻ mask[2],"\t  ",bitstring(p.⊻mask[2])[end-n:end],
            "\t\t  ",m[p+1],"\t  ",m[p.⊻mask[2]+1]) 
    else
        println(p,'\t',bitstring(p)[end-n:end],'\t',
            p .⊻ mask[2],'\t',bitstring(p.⊻mask[2])[end-n:end],
            "\t\t",m[p+1],"\t",m[p.⊻mask[2]+1]) 
    end
end
```

      Psi 	Psi 	Masked 	 Masked	 	mPsi  	mMasked
      0	00000	6	00110		0.0	2.0
      1	00001	7	00111		1.0	3.0
        2	  00010	  4	  00100		  1.0	  1.0
        3	  00011	  5	  00101		  2.0	  2.0
        4	  00100	  2	  00010		  1.0	  1.0
        5	  00101	  3	  00011		  2.0	  2.0
      6	00110	0	00000		2.0	0.0
      7	00111	1	00001		3.0	1.0
      8	01000	14	01110		1.0	3.0
      9	01001	15	01111		2.0	4.0
        10	  01010	  12	  01100		  2.0	  2.0
        11	  01011	  13	  01101		  3.0	  3.0
        12	  01100	  10	  01010		  2.0	  2.0
        13	  01101	  11	  01011		  3.0	  3.0
      14	01110	8	01000		3.0	1.0
      15	01111	9	01001		4.0	2.0



That's it for background.  

In the next few posts, I'll first cover breaking the matrix up according to the magnetization symmetry, and then the sorting and the searching that accompanies that change.  Then we will actually diagonalize the matrix and looks at the results.  

In the meantime, go back and check out my posts on <b>"Jacobi Transformation of a Symmetric Matrix"</b> and <b>"Julia with MKL on OSX"</b> for some information on this topic.  

Happy physicsing :)


[1]  [Field Theories of Condensed Matter Physics](https://www.amazon.com/Field-Theories-Condensed-Matter-Physics/dp/0521764440) by Fradkin.  


<table class="tagtable" style="table-layout: fixed">
    <tr><th class="cattype"  colspan="4">1D Spin Chain</th></tr>
        <tr>
          <th class="thprop" style="width: 1em">#</th>
          <th class="thtitle" style="width: auto">Title</th>
          <th class="thprop" style="width: 4em">Level</th>
          <th class="thprop" style="width: 5em">Tags</th>
        </tr>
    <tbody>
    <tr class="trtag">
        <td class="tagprop">1</td>
        <td class="tagtitle"><h3><a href="{{base.url}}/M4/Jacobi-Transformation.html">Jacobi Transformation</a></h3></td>
        <td class="tagprop"><p class="post-meta">Numerics</p></td>
        <td class="tagprop"><p class="post-meta">ED</p></td>
    </tr>
    <tr class="trtag">
        <td class="tagprop">2</td>
        <td class="tagtitle"><h3><a href="{{base.url}}/M4/Spin-Chain-Prerequisites.html">1D Spin Chain Prerequisites</a></h3></td>
        <td class="tagprop"><p class="post-meta">Graduate</p></td>
          <td class="tagprop"><p class="post-meta">Magnet Quantum ED</p></td>
    </tr>
    <tr class="trtag">
        <td class="tagprop">3</td>
        <td class="tagtitle"><h3><a href="{{base.url}}/M4/Spin-Chain-pt2.html">1D Spin Chain Values and Vectors</a></h3></td>
        <td class="tagprop"><p class="post-meta">Graduate</p></td>
            <td class="tagprop"><p class="post-meta">Magnet Quantum ED</p></td>
    </tr>
</tbody>
</table>