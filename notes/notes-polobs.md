# Notes on extracting polarization observables

*	[11-14-13](#11-14-13)
	*	Formalism

<h2 id="11-14-13">11-14-13</h2>

##Formalism

$\left(\frac{d\sigma}{dX^{ij}d\phi^{j}}\right)^{h}
\doteq
f^{h}(X^{ij},\phi^{j}) = A^{ij} +
						B^{ij}\cos\phi^{j} +
						C^{ij}\cos2\phi^{j} +
						hPD^{ij}\sin\phi^{j}$

where

*	ij = index over Varset,Variable (3x5 matrix)
*	$R2^{ij}_{\alpha} \doteq 
	[A^{ij},B^{ij},C^{ij},D^{ij}] \equiv 
	[R_{T}+\epsilon_{L}R_{L}, R_{LT}, R_{TT}, R_{LT'}]$
	*	$R2^{ij}_{\alpha} = f(Q^{2},W,X^{ij})$	

##Event Selection

1. `eid`
2. `efid`
3. `momcorr`
4. `MM cut`

##`R2` Extraction Method

Of the methods listed earlier:

1.	Fit $f^{h}(X^{ij},\phi^{j})$ to extract `R2`
2.	Calculate Asymmetry $\doteq$ $f^{h=+}-f^{h=-}$ and then extract $D^{ij}$
3.	$\int f^{h}(X^{ij},\phi^{j}) * (\cos\phi/\cos 2\phi/\sin\phi)d\phi$ to extract $B^{ij}/C^{ij}/D^{ij}$

Method 3. is used, which even at the level of algorithmic detail is listed below. $\color{red} \text{NOTE that when multiplying by $\sin\phi$, the sign of the polarization is explicity used}$

For every `q2wbin`:

1.	`h5[pol]` where `pol` $\in$ {POS,NEG,UNP,AVG}; `pol` $\neq$ AVG
2.	`h5m[pol,pob]` = `h5[pol]`$\cdot$ `h5f[pob]` 
	*	`pob` $\in$ {A,B,C,D}; `pol` $\neq$ AVG
	*	`h5f[pob]`:
		*	 For every bin `i`,  `h5f[pob](i)` = `f[pob](i)`
		*	 `f[pob]` $\in$ {N.A.,$\cos\phi$,$\cos 2\phi$,$\color{red}{\text{sign(pol)}}$ $\sin\phi$}
3.	`hR2_Xij[pol,pob]` = `h5m[pol,pob]` `Project` on to $X^{ij}$; `pol` $\neq$ AVG
4. 	`hR2_Xij[pol=AVG,pob]` = (`hR2_Xij[pol=POS,pob]` + `hR2_Xij[pol=NEG,pob]`)/2

##Notes on current Observations

1. `R2_Xij[UNP,D]` $\neq$ 0 for Simulation. Why?
2. `R2_Xij[UNP,D]` $\neq$ 0 for Experiment.






This is a normal paragraph:

This is also

	This is a code block

Use `printf`




