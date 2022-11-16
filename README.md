# Anharmonic oscillator

## Introduction 

This program calculates the lowest $10$ eigenvalues of an anharmonic oscillator system perturbatively up to $100$-th perturbation order. The Hamiltonian is $H(x)=p^2+(x^2-\epsilon x^4)$ where $\epsilon$ is the perturbation parameter. The perturbation is summed with Padé approximant. 

Alternative to the analytical approach, direct solvers based on diagonalization in energy representation and shooting method are included. 

## Prerequisite 

1. [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)
2. [GMP](https://gmplib.org/)

## How to use 

### Generate perturbation coefficient 

1. Compile as follows to calculate for the $n$-th eigenvalues, replace the ``n`` in option ``-DOn`` with ``0``, ``1``,..., or ``9``. Only $n=0,1, \ldots, 9$ are supported, otherwise compile error will arise. 

```bash
g++ energy.cpp -fdiagnostics-color=always -std=c++17 -lm -lgmp -lgmpxx -DOn -DNDEBUG -o PERTURBATION_EXECTUABLE_NAME
```

Compilation usually takes $3$ minutes to complete. 

2. Run the output executable as usual. You can redirect the output to file. To help you check answers, in ``results/`` coefficients for all $10$-orders are stored.  

```bash
./PERTURBATION_EXECTUABLE_NAME > PERTURBATION_FILE_NAME
```

### Generate Padé approximant 

1. Run ``pade_generator.py`` that transforms the coefficients generated from previous steps into code. 

```bash
python3 pade_generator.py
```

2. Compile ``pade.cpp``. 

```bash
g++ pade.cpp -fdiagnostics-color=always -std=c++17 -lm -lgmp -lgmpxx -O2 -DNDEBUG -o PADE_EXECTUABLE_NAME
```

3. Run the executable. 
```bash
./PADE_EXECTUABLE_NAME
```
The output would be calculated eigenvalues with $\epsilon$ from $-10$ to $10$ at a step of $0.1$. 

### Numerical calculation with  diagnalization 

1. Compile ``bf_numerical.cpp``. 

```bash
g++ bf_numerical.cpp -fdiagnostics-color=always -std=c++17 -lm -lgmp -lgmpxx -O2 -DNDEBUG -o DIAG_EXECTUABLE_NAME
```

2. Run ``DIAG_EXECTUABLE_NAME`` with one argument. This argument is the value of $\epsilon$. To give an example, 

```bash
./DIAG_EXECTUABLE_NAME 0.1
```
The output would be $10$ numbers representing the lowest $10$ eigenvalues.

### Numerical calculation with shooting method 

**N.B.** Current implementation works only with $\epsilon < 0$.  

1. Compile ``rk4.cpp``. Depending on the parity of the unperturbed state, compile with 

```bash
g++ rk4.cpp -fdiagnostics-color=always -std=c++17 -lm -lgmp -lgmpxx -O2 -DNDEBUG -DBVC_SYM -o SHOOT_EXECUTABLE_NAME
```

for symmetric state or


```bash
g++ rk4.cpp -fdiagnostics-color=always -std=c++17 -lm -lgmp -lgmpxx -O2 -DNDEBUG -DBVC_ASYM -o SHOOT_EXECUTABLE_NAME
```

for anti-symmetric state. 

2. Run the executable with $3$ arguments: $\epsilon$, lower and upper guess of eigenvalue. Suppose you compiled with ``BVC_SYM``, then an example call would be 

```bash
./SHOOT_EXECUTABLE_NAME -0.1 1.0 1.1
```

and the output is ``1.065285873413086``. 