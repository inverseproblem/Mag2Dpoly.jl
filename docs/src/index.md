


# Contents

```@contents
Pages = ["index.md"]
Depth = 4
```

## User guide

**Mag2Dpoly** *is a Julia package conceived for forward magnetic anomaly calculation due to two-dimensional polygonal bodies with uniform arbitrary polarization*. 

The formulations implemented in this package are that of Talwani & Heirtzler (1962, 1964), Won & Bevis (1987) and revised Kravchinsky et al. (2019).

If you use this code for research or else, please cite the related paper:

Alessandro Ghirotto, Andrea Zunino, Egidio Armadillo & Klaus Mosegaard (2021). **Magnetic Anomalies Caused by 2D Polygonal Structures with Uniform Arbitrary Polarization: new insights from analytical/numerical comparison among available algorithm formulations**. *Geophysical Research Letters, 48*(7), e2020GL091732.

The specific procedures for each formulation, the analytical/numerical results derived from their comparison and the rectification made to Kravchinsky et al. (2019) algorithm are describde in detail in the paper above.

## Documentation

```@meta
Author = "Andrea Zunino"
Author = "Alessandro Ghirotto"
```

### Installation

To install the package simple enter into the package manager mode in Julia by typing "`]`" at the 
REPL prompt and then use `add`, i.e.,
```
(v1.5) pkg> add Mag2Dpoly
```
The package will be automatically downloaded from the web and installed.

!!! warning
    At the moment the package is not yet registered in the official Julia registry, so, 
    to install it run the following in package mode:
		
    ```julia
    (v1.5) pkg> add https://github.com/inverseproblem/Mag2Dpoly.jl
    ```
	
Alternatively, use the path where the directory of the package is located, be it local or remote (Github):
```
(v1.5) pkg> add /path/to/Mag2Dpoly.jl
```

### Theoretical Background

For a theoretical explanation, let us consider a three-dimensional non-magnetic 
space in which a body infinitely extended in the ``y`` direction is immersed. 

The common aim of all formulations is the calculation of the magnetic field of 
this body upon an observation point ``(x_0,z_0)`` located along a profile aligned to 
the ``x`` direction (the positive ``z`` axis is assumed pointing downward).

The starting assumption is that our body can be considered as discretized by an 
infinite number of uniformly-magnetized elementary volumes with infinitesimal dimensions ``dx``, ``dy``, ``dz``.

Within this assumption, the magnetic field associated to the body can be mathematically 
expressed in terms of a line integral around its periphery, represented in two dimensions 
as its polygonal cross-section (in red).
![](images/intro.svg)

### Tutorial
First load the module and define some magnetization vectors,
```@example ex1
using Mag2Dpoly 

# induced magnetization
Jind = MagnetizVector(mod=4.9,Ideg=90.0,Ddeg=45.0)
# remanent magnetization
Jrem = MagnetizVector(mod=3.1,Ideg=45.0,Ddeg=0.0)
nothing # hide
```
and then define some observation points 
```@example ex1
# angle with the North axis
northxax = 90.0

# number of observatoin 
N=101
xzobs = hcat(LinRange(0.0,100.0,N), -1.0*ones(N))
nothing # hide
```
Finally the general list of vertices of the poligonal bodies and the relative indices mapping each body to its vertices:
```@example ex1
# vertices of the poligonal bodies
vertices  = [35.0 50.0;
             65.0 50.0;
             80.0 35.0;
             65.0 20.0;
             35.0 20.0;
             20.0 35.0]
			 
# indices of vertices for the body
ind1 = collect(1:6)
bodyindices = [ind1]
# construct the poligonal body object
pbody = MagPolygBodies2D(bodyindices,vertices)
```

At this point the total field can be computed. We select `"talwani"` as the forward type:
```@example ex1
# type of forward algorithm
forwardtype = "talwani"
# arrays of magnetization vectors
Jinds = [Jind]
Jrems = [Jrem]
# compute total field 
tmag = tmagpolybodies2Dgen(xzobs,Jinds,Jrems,northxax,pbody,forwardtype)
```

Now we can plot the results:

    using PyPlot
    xmi=minimum(xzobs[:,1]) 
    xma=maximum(xzobs[:,1])
    figure()
    subplot(211)
    plot(xzobs[:,1],tmag,".-") 
    xlim(xmi,xma)
    subplot(212)
    x = [pbody.bo[1].ver1[:,1]...,pbody.bo[1].ver2[end,1]]
    y = [pbody.bo[1].ver1[:,2]...,pbody.bo[1].ver2[end,2]]
    plot(x,y,"o-")
    xlim(xmi,xma)
    gca().invert_yaxis()

![](images/plotex1.svg)

## Public API
```@docs
Mag2Dpoly
```

### Data structures
```@docs
BodySegments2D
```

!!! warning 
    Vertices of the polygonal bodies must be provided 
    counterclockwise to the function `BodySegments2D`
    to perform magnetic anomaly calculation using the
    functions in the next section **Forward functions**


```@docs
MagPolygBodies2D
MagnetizVector
```

### Forward functions
#### Single polygonal body
```@docs
tmagpolybodies2D
tmagpolybodies2Dgen
```
#### Multiple polygonal bodies
```@docs
tmagpoly2D
tmagpoly2Dgen
```
#### Forward algorithms alone
!!! note
    These functions are not exported. To call them 
    type `Mag2Dpoly.` before the name of the functions.
	
```@docs
Mag2Dpoly.tmagtalwani
Mag2Dpoly.tmagtalwanired
Mag2Dpoly.tmagkrav
Mag2Dpoly.tmagwonbev
```

 
### Useful functions
!!! note
    These functions are not exported. To call them
    type `Mag2Dpoly.` before the name of the functions.
	
```@docs
Mag2Dpoly.convert_H_to_B_nT
Mag2Dpoly.convert_B_nT_to_H
Mag2Dpoly.magcomp
Mag2Dpoly.checkanticlockwiseorder
```
