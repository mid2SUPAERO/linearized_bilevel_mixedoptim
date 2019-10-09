Works with
- python 2.7
- numpy 1.11.3
- scipy 0.19.0

This folder contains :

- main.py :
        Call to solver, definition of input parameters.
- bilevel.py :
        Bi-level solver procedure.
- categorical.py :
        Categorical data : each choice corresponds to a material and a stiffener,
        It is possible to change / add new catalogs.
- lower_pb.py :
        Implementation of Psi(c) function (lower level problem),
        Definition of constraints, objective, call to scipy continuous solver.
- truss_3b :
        Implementation of tress and displacements computations for 3-bar truss structure,
        Implementation of weight computation
        Implementation of buckling limit stress computation (to be adapted if change of stiffener profile)