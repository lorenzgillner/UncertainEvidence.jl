# UncertainEvidence.jl

This package provides methods for handling epistemic uncertainty, primarily implementing the [Dempster-Shafer Theory of Evidence (DST)](https://en.wikipedia.org/wiki/Dempster%E2%80%93Shafer_theory). It is, to some extend, a port of the [IPP toolbox](https://www.uni-due.de/il/ipptoolbox.php) to the Julia language.

This package is still in an early development phase, but the basic principles of the Dempster-Shafer framework are already present.

## Usage

Let's consider a simple introductory example. Two experts, $X_1$ and $X_2$, state their opinions on the possible causes (*focal elements*; $A$, $B$ and $C$) of a situation. Using UncertainEvidence.jl, 

```julia
using UncertainEvidence
```

these opinions can be formulated in the following way:

```julia
# expert one
X1 = BPA(
    'A' => 0.5,
    'B' => 0.3,
    'C' => 0.2
)

# expert two
X2 = BPA(
    'A' => 0.6,
    'B' => 0.1,
    'C' => 0.3
)
```

The basic probability assessment (`BPA`) is the fundamental data structure for computations. To generate a general assessment of the situation, `BPA`s can be combined with Dempster's rule of combination (more combination rules coming soon):

```julia
X12 = combine_dempster(X1, X2)
```

Based on this combined structure, the lower (*belief*) and upper bound (*plausibility*) of the likelihood of a cause can be queried:

```julia
# lower bound
bel('A', X12)

# upper bound
pls('A', X12)
```

The masses of a `BPA` should always sum up to 1. To ensure this, use the  method `bpa` instead:

```julia
X = bpa(
    Set('A') => 0.1,
    Set('B') => 0.2,
    Set('C') => 0.3
)

# Dict{Set{Char}, Float64} with 4 entries:
#   Set(['C', 'A', 'B']) => 0.4
#   Set(['B'])           => 0.2
#   Set(['A'])           => 0.1
#   Set(['C'])           => 0.3
```

Note the use of `Set`s in this example; if the sum of masses is lower than 1, the remaining mass is automatically assigned to $\Omega$, the set of all focal elements. This would fail, if all other focal elements would be of type `Char`.