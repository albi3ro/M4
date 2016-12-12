# Exact Diagonalization Prerequisites

## Christina Lee

## Category: Grad

I was going to try "Exact Diagonalization of a 1D Heisenburg Spin Chain" all in one post. Then I realized, even for a grad level post, the post started to get a bit hefty. So here is <b>Part 1</b> of a currently four part series. I haven't posted in a while, partially because I've been needing to work on my thesis research, partially because this hasn't been an easy post to write.

If you are not a physics graduate student, still give this series a read and try to not get too intimidated. If you are a physics graduate student, still give this series a read and try not to get too intimidated.

For further information, checkout this set of lecture notes on the arxiv, [Computational Studies of Quantum Spin Systems](https://arxiv.org/abs/1101.3281) by Anders W. Sandvik.

Today's just theory. Sorry, no programming :(

## So what is a 1D Heisenburg Spin $\frac{1}{2}$ Chain? 
Let's break down each part of that phrase.  

#### Spin $\frac{1}{2} $
At each site, we can have a particle pointing up $| \uparrow \rangle$, down $|\downarrow \rangle$, or some super-position of the two.  

### Heisenburg
Our spin has three degrees of freedom and full $SU(2)$ symmetry, the most general case. Once we can solve the physics of the Heisenburg case, we also have solved the <i>XY model</i> ($J_z=0$) and the <i>Ising model</i> ($J_x=J_y=0$).

### 1D Chain
Spin only couples to two neighbors.

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
S^x=\begin{bmatrix}
0 & 1 \\
1 & 0 \\
\end{bmatrix}
\;\;\;\;\;\;\;
S^y=\begin{bmatrix}
0	&	-i	\\
i	&	0	\\
\end{bmatrix}
\;\;\;\;\;\;\;
S^z=\begin{bmatrix}
1	&	0	\\
0	&	-1	\\
\end{bmatrix},
$$
and $S^{\pm}$ the ladder operators
$$
S^+=\frac{S^x+i S^y}{2} = \begin{bmatrix}
0 & 1 \\
0 & 0\\
\end{bmatrix}
\;\;\;\;\;\;\;
S^- =\frac{S^x-i S^y}{2}= \begin{bmatrix}
0 & 0 \\
1 & 0\\
\end{bmatrix}.
$$
## Analytical Solutions

    "The model is solvable, not the solution"

The Ising, XY, and Heisenburg cases all fall into the special class of problems where they have exact solutions in the infinite size limit.  Very few exact solutions of Quantum problems exist, so we try and get as much milage as we can out of the couple of ones we have.  

Interestingly enough, once one solution to a problem comes along, someone figures out a different way to approach the same problem.

The 1D Quantum Ising Model is equivalent to the 2D Ising model of classical statistical mechanics, solved in 1944 by Lars Onsager.  The solution is also equivalent to a description of free Majorana fermions.  For more information, see [Field Theories of Condensed Matter Physics](https://www.amazon.com/Field-Theories-Condensed-Matter-Physics/dp/0521764440) by Fradkin.  

The Jordan-Wigner Transformation solves the 1D XY model by mapping spins to fermions.  This transformation only works in special 1D circumstances.  Since spins possess different anti-commutation relationships than fermions, we attach a string of operators stretching to infinity to each spin, thus correcting the problem.  After the transformation, we get a Hamiltonian that is quadratic in the momentum fermionic operators $d_k$, $d^{\dagger}_k$, and we can see the $\cos (ka)$ dispersion relationship for the excitations,
$$
{\cal \tilde{H}}_{XY}=-J \sum_k \cos (ka) d^{\dagger}_k d_k.
$$
I might write a full article on this later.  Performing this transformation on the full Heisenburg model gives a four-operator scattering term, which cannot easily be dealt with.  

The Bethe Ansatz, which I honestly don't know anything about, solves the full 1D Heisenburg Model, as well as some 1D Bose gas and Hubbard model problems.  Come back to me in many years, or ask a Russian mathematician.  

## Noether's Theorem: Symmetries and Conserved Quantities

 If you meet anyone who doesn't believe women can be great at mathematics or science, point them to Emily Noether, the genius behind perhaps the most beautiful idea in modern physics.  For each symmetry, there exists a conserved quantity.  Space translational symmetry gives us momentum conservation.  Time translational symmetry gives us energy conservation.  
 
 These save us so much in classical mechanics, but saves us even more in condensed matter physics.  By working out the symmetries of the problem, we can break the huge Hilbert space into smaller chunks of conserved quantities and only have a smaller problem to deal with.  
 
 So what are our conserved quantities? Magnetization and momentum 























