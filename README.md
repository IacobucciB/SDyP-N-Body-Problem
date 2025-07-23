# SDyP-N-Body-Problem

## Overview
This repository contains a comprehensive implementation of the N-Body Problem, a fundamental problem in physics and computer science that involves predicting the motion of multiple celestial bodies under mutual gravitational interaction. The problem is computationally intensive and serves as a benchmark for evaluating numerical methods, performance optimization, and parallel computing.

The project includes:

+ A sequential baseline implementation.
+ A parallel implementation using POSIX Pthreads for shared-memory systems.
+ A hybrid MPI + Pthreads implementation for distributed-memory systems.
+ A visual simulator using Raylib to graphically display the simulation in real time.
+ A set of validation tools to ensure the correctness of parallel implementations against the sequential reference.

This project was developed as part of the Parallel and Distributed Systems course at the National University of La Plata.

## Repository Structure

The repository is organized into several directories, each serving a specific purpose in the development, execution, and validation of the N-Body simulation.

* **`/src`** Contains the main source code files written in C:
  * `sequential.c`: Sequential implementation of the N-Body simulation. Serves as the baseline for parallelization and validation.
  * `pthreads.c`: Shared-memory parallel implementation using POSIX Pthreads.
  * `mpi.c`: Hybrid implementation using MPI for distributed memory and Pthreads for multithreading on each node.

* **`/run-local`** Includes Bash scripts for compiling and running the implementations in the `src` directory locally. These scripts allow quick configuration changes such as the number of bodies, threads, and simulation steps for testing and benchmarking purposes.

* **`/example`** Contains an example provided by the university as a reference when solving the problem. This serves as a guideline or comparison point.

* **`/scripts`** Contains Bash scripts used to submit and run simulations on the university’s computing cluster. This directory also stores the output results from those cluster executions.

* **`/validation`** Includes Python scripts that validate the output of the parallel programs (`pthreads` and `mpi`) by comparing their results to the sequential implementation. These scripts help ensure numerical accuracy and correctness.

* **`/simulation`** Contains a real-time simulation built with Raylib and Pthreads. This component provides a graphical visualization of the system evolution.

  > 📌 *Raylib is used under its open-source license. Please see [Raylib's website](https://www.raylib.com/) for more information.*

## Dependencies
+ C Compiler (GCC or Clang)
+ POSIX Pthreads
+ OpenMPI (for mpi.c)
+ Python 3 (for validation scripts)
+ Raylib (for the 3D simulation)