# numerical-methods-extension-odes

MATLAB repository featuring implementations of advanced and implicit numerical methods for solving ordinary differential equations (ODEs), symbolic order-of-convergence analysis, construction of Runge-Kutta collocation formulas, and practical application examples in physical, biological, and ecological systems.

---

## Description

This project extends the classical study of solving Initial Value Problems (IVPs) for ODEs. It includes:
* **Implicit and Multistep Numerical Methods**: Implementations of Adams-Bashforth, Adams-Moulton, Milne-Simpson, Nyström, Implicit Runge-Kutta, and Implicit Taylor methods[cite: 137].
* **Symbolic Order Analysis**: Automated functions based on rooted tree analysis and Taylor series to symbolically determine the order of accuracy of numerical methods[cite: 137].
* **Collocation Methods**: Tools for constructing Runge-Kutta methods using polynomial collocation[cite: 137].
* **Practical Modeling Examples**: Real-world case studies applied to mechanical systems, ecological compartment dynamics, and tumor growth models[cite: 137].

---

## Repository Structure

Below is a detailed overview of each file and module within the repository[cite: 137]:

| File / Module | Description |
| :--- | :--- |
| `ABashforth1.m` | Implementation of the 1st-order explicit Adams-Bashforth multistep method[cite: 137]. |
| `AMoulton1.m` | Implementation of the 1st-order implicit Adams-Moulton multistep method[cite: 137]. |
| `MilneSimpson1.m` | Implementation of the 1st-order Milne-Simpson multistep method[cite: 137]. |
| `Nystrom1.m` | Implementation of the 1st-order explicit Nystrom multistep method[cite: 137]. |
| `RungeKuttaImplicito1.m` | Scheme for the 1st-order Implicit Runge-Kutta method[cite: 137]. |
| `TaylorImplicito1.m` | Implementation of the 1st-order Implicit Taylor numerical method[cite: 137]. |
| `calcularCoefRKOrden.m` | Analytical/symbolic calculation of Runge-Kutta coefficients for a target order[cite: 137]. |
| `calcularFRK.m` | Calculation of the system of implicit equations in Runge-Kutta methods[cite: 137]. |
| `calcularJRK.m` | Jacobian matrix computation for implicit Runge-Kutta solvers[cite: 137]. |
| `calcularPsiTaylor.m` / `calcularJyPsiTaylor.m` | Evaluation of the auxiliary function $\Psi$ and its Jacobian for the implicit Taylor method[cite: 137]. |
| `coefRKImplicito.m` | Definition of Butcher tableaux and coefficients for implicit RK schemes[cite: 137]. |
| `comprobarOrdenRKImp.m` | Analytical/symbolic verification of the accuracy order achieved by an implicit RK method[cite: 137]. |
| `comprobarOrdenTaylorImp.m` | Symbolic verification of the order of accuracy in implicit Taylor methods[cite: 137]. |
| `construirRKColocacion.m` | Automatic generation of Runge-Kutta methods using polynomial collocation nodes[cite: 137]. |
| `crearArbol.m` / `combinarHijos.m` / `generarArbolesHastaOrden.m` | Generation and manipulation of rooted trees for order analysis in Runge-Kutta methods[cite: 137]. |
| `comprobarDimensiones.m` | Dimension validation function for state vectors and matrices[cite: 137]. |
| `esCeroOrdenRK.m` / `esCeroSimbolico.m` | Auxiliary functions for evaluating zero residuals in symbolic expressions[cite: 137]. |
| `extraerValoresSolucion.m` | Extraction and formatting of solution data during numerical integration[cite: 137]. |
| `generarPhik.m` | Generation of $\Phi_k$ coefficients associated with order conditions[cite: 137]. |
| `leerOpciones.m` | Management of input options and tolerances for numerical solvers[cite: 137]. |
| `ejemplo1_masa_resorte.m` | Simulation script for the Mass-Spring dynamic system[cite: 137]. |
| `ejemplo2_compartimentos_lagos.m` | Simulation script for the lake compartmental pollution model[cite: 137]. |
| `ejemplo3_tumor_linealizado.m` | Simulation script for the linearized tumor growth model[cite: 137]. |
| `main.pdf` | Explanatory supporting document covering theoretical details and derivations[cite: 137]. |
| `LICENSE` | Project distribution license[cite: 137]. |

---

## Requirements and Installation

1. **MATLAB** (R2020a or newer recommended)[cite: 137].
2. **Symbolic Math Toolbox** (required for symbolic order checking functions and tree analysis)[cite: 137].

### Installation:
Clone this repository or download the compressed source code[cite: 137]:
```bash
git clone [https://github.com/your-username/numerical-methods-extension-odes.git](https://github.com/your-username/numerical-methods-extension-odes.git)
cd numerical-methods-extension-odes
```
## Running the Examples

You can execute any of the example scripts directly in the MATLAB Command Window:

```matlab
% Mass-Spring system simulation
run('ejemplo1_masa_resorte.m')

% Lake pollution compartment model simulation
run('ejemplo2_compartimentos_lagos.m')

% Linearized tumor growth model simulation
run('ejemplo3_tumor_linealizado.m')
```
