# Itovsn3: A Mathematica package for Symbolic Ito calculus

This repository contains the source code for the
[`FernandoDuarte/Itovsn3`](https://resources.wolframcloud.com/PacletRepository/resources/FernandoDuarte/Itovsn3/)
paclet.

## Installation

`FernandoDuarte/Itovsn3` is available for installation from the [Wolfram Paclet Repository][https://resources.wolframcloud.com/PacletRepository/]:


To install this paclet in your Wolfram Language environment, evaluate this code:

```wolfram
PacletInstall[ResourceObject["FernandoDuarte/Itovsn3"]]
```

To load the code after installation, evaluate this code:

```wolfram
Needs["FernandoDuarte`Itovsn3`"]
```

## Documentation

* [Main guide](https://resources.wolframcloud.com/PacletRepository/resources/FernandoDuarte/Itovsn3/guide/Main.html)
* Tech notes
  - [Stochastic Integration](https://resources.wolframcloud.com/PacletRepository/resources/FernandoDuarte/Itovsn3/tutorial/StochasticIntegration.html)
  - [Bessel processes](https://resources.wolframcloud.com/PacletRepository/resources/FernandoDuarte/Itovsn3/tutorial/Bessel.html)
  - [Derivation of the Black-Scholes formula](https://resources.wolframcloud.com/PacletRepository/resources/FernandoDuarte/Itovsn3/tutorial/BlackScholes.html)
  - [Distribution of the Levy stochastic area](https://resources.wolframcloud.com/PacletRepository/resources/FernandoDuarte/Itovsn3/tutorial/ItoArea.html)
  - [Mardia-Dryden distribution](https://resources.wolframcloud.com/PacletRepository/resources/FernandoDuarte/Itovsn3/tutorial/MardiaDryden.html)
  - [Coupled pairs of Brownian motions reflecting off a half-plane](https://resources.wolframcloud.com/PacletRepository/resources/FernandoDuarte/Itovsn3/tutorial/Reflect.html)

# Original Source and Attributions
`Itovsn3` is a package created by [Wilfrid S. Kendall](http://www.warwick.ac.uk/go/WSK), emeritus professor in the
Statistics department at the University of Warwick.

The original source code for version 3.70 of the package is in the `itovsn3-mma`
folder of this repository. The rest of the files and folders are the paclet version of the same
code with minor modifications to comply with paclet requirements, put together
by [Fernando Duarte](https://fernandoduarte.github.io), Department of Economics,
Brown University.

## References
Wilfrid S. Kendall: "Doing stochastic calculus with Mathematica", in Economic and Financial Modeling with Mathematica, edited by H. Varian, Springer-Verlag, New York, pp 214-238 (1993).