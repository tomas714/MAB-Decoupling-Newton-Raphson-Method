# MAB Decoupling - Newton-Raphson method

This project implements a **decoupling algorithm for a Multi-Active Bridge (MAB) converter** using the **Newton-Raphson method**. The code is written in **C** and uses only **standard C libraries**.

## Table of Contents

- [Project Description](#project-description)
- [Requirements](#requirements)
- [Usage](#usage)
- [License](#license)
- [Author](#author)

## Project Description

The goal of this project is to numerically decouple a MAB converter, which is a crucial step in the design of control algorithms for power electronic systems. The Newton-Raphson method is used to solve the nonlinear equations describing the system's steady-state behavior.

## Requirements

To build and run this project, you will need the following tools installed on your system:

- `gcc` – GNU Compiler Collection for compiling C code [[Home page](https://www.gnu.org/software/gcc/)]
- `make` – Build automation tool to manage compilation [[Home page](https://www.gnu.org/software/make/)]

### Build project
```bash
make all
```

### Clean project
```bash
make clean
```

## Usage

After building, an executable is generated. Input parameters and usage instructions can be found in the `main.c` file or relevant source code comments.

## License

This project is licensed under the terms of the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html).

## Author

**Tomáš Basarik** \
Department of Electrical Engineering and Mechatronics \
Technical University of Košice \
Košice, Slovakia