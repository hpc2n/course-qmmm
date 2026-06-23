# QM/MM Methods for Biomolecular Simulations


# Learning Objectives

By the end of this workshop participants should be able to:

1. Explain why QM/MM methods are needed.
2. Distinguish between QM, MM, and coarse-grained models.
3. Understand the QM/MM Hamiltonian.
4. Describe common QM/MM embedding schemes.
5. Select an appropriate QM region.
6. Recognize practical issues in QM/MM calculations.
7. Evaluate the advantages and limitations of QM/MM simulations.


**Duration:** 2 hours  
**Instructor:** Pedro Ojeda-May, Application Expert, HPC2N  
*Based on slides by Kwangho Nam, University of Texas at Arlington*


## Table of Contents

1. [Motivation: Why QM/MM?](#1-motivation-why-qmmm)
2. [Overview of Computational Methods](#2-overview-of-computational-methods)
3. [The Hybrid QM/MM Potential](#3-the-hybrid-qmmm-potential)
4. [Solving the QM Hamiltonian](#4-solving-the-qm-hamiltonian)
5. [Schemes for the QM/MM Hamiltonian](#5-schemes-for-the-qmmm-hamiltonian)
6. [QM/MM Coupling: Embedding Schemes](#6-qmmm-coupling-embedding-schemes)
7. [Electrostatic Embedding in Detail](#7-electrostatic-embedding-in-detail)
8. [Practical Issue I: QM Region and Method Selection](#8-practical-issue-i-qm-region-and-method-selection)
9. [Practical Issue II: QM/MM Boundary Treatment](#9-practical-issue-ii-qmmm-boundary-treatment)
10. [Practical Issue III: Periodic vs. Non-Periodic Boundary](#10-practical-issue-iii-periodic-vs-non-periodic-boundary)
11. [Free Energy Methods Combined with QM/MM](#11-free-energy-methods-combined-with-qmmm)
12. [Speed and Parallelization](#12-speed-and-parallelization)
13. [Speeding up AI-QM/MM: Multiscale Approach](#13-speeding-up-ai-qmmm-multiscale-approach)
14. [Summary and Practical Recommendations](#14-summary-and-practical-recommendations)


## 1. Motivation: Why QM/MM?

Enzymes are extraordinarily efficient catalysts. A central goal in computational biochemistry is to understand *how* they achieve such remarkable rate enhancements. The Michaelis–Menten scheme captures the key steps:

$$E + S \underset{k_2}{\stackrel{k_1}{\rightleftharpoons}} ES \xrightarrow{k_\text{cat}} E + P$$

The **catalytic efficiency** is defined as $k_\text{cat}/k_\text{uncat}$, where $k_\text{uncat}$ is the rate of the uncatalyzed reaction. Some examples:

| Enzyme | $k_\text{cat}/k_\text{uncat}$ |
|---|---|
| OMP decarboxylase | $\sim 10^{17}$ |
| $\beta$-Amylase | $\sim 10^{17}$ |
| Fumarase | $\sim 10^{15}$ |
| Carbonic anhydrase | $\sim 10^{7}$ |

Enzymes lower the activation free energy $\Delta G_C$ relative to the uncatalyzed barrier $\Delta G_U$, shifting the transition state from $E+S^\ddagger$ to $ES^\ddagger$.

> **Key question:** What molecular interactions drive this rate enhancement? Answering this requires an atomistic description of bond breaking and formation — which purely classical force fields cannot provide.



## 2. Overview of Computational Methods

### 2.1 Potential Energy Functions

| Category | Methods | Cost scaling |
|---|---|---|
| Molecular Mechanics (MM) | CHARMM, AMBER, MM2 | $\mathcal{O}(N \log N)$ with PME |
| Coarse-Grained (CG) | Go-potential, Martini | < MM |
| Semi-empirical QM | AM1, PM3, SCC-DFTB, xTB, EVB | $\mathcal{O}(N^3)$, very fast |
| DFT | B3LYP, PBE, M06-2X | $\mathcal{O}(N^3)$–$\mathcal{O}(N^4)$ |
| Ab initio MO | HF, MP2, CCSD(T), FCI | $\mathcal{O}(N^4)$–$\mathcal{O}(N^7)$ |
| Hybrid | QM/MM, MM/CG, QM/MM/CG | intermediate |

**Relative cost:** Full QM $\gg$ QM/MM $\gg$ MM $>$ CG

### 2.2 Statistical Simulation Methods

- **Molecular Dynamics (MD):** Newtonian (Verlet, Velocity Verlet), Langevin, Car–Parrinello; ensembles NPT, NVT, NVE.
- **Monte Carlo:** Metropolis method, Gibbs ensemble MC.

### 2.3 Searching the Potential Energy Surface

- Energy minimizations: steepest descent (SD), conjugate gradient (CG), Newton–Raphson (NR).
- Identification of reactant, transition state, and product configurations.
- Conformational search.



## 3. The Hybrid QM/MM Potential

The central idea: treat the chemically active region with QM and the environment with the cheaper MM force field.

- **QM where needed** — small reactive region.
- **Low cost** — the QM region typically contains fewer than 200 atoms.
- **Bond formation/breaking** is naturally described at the QM level.

The effective Hamiltonian is partitioned as:

$$\hat{H}_\text{eff} = \hat{H}_\text{QM} + \hat{H}_\text{QM/MM} + \hat{H}_\text{MM}$$

where the QM/MM coupling term is further decomposed into electrostatic, van der Waals, and boundary contributions:

$$\hat{H}_\text{QM/MM} = \hat{H}_\text{QM/MM}^\text{elec} + \hat{H}_\text{QM/MM}^\text{vdW} + \hat{H}_\text{QM/MM}^\text{boundary}$$

- $\hat{H}_\text{QM}$: full quantum mechanical Hamiltonian for the QM region.
- $\hat{H}_\text{MM}$: classical force field energy for the MM region.

> **Reference:** Warshel & Levitt, JMB (1976); Field, Bash & Karplus, JCC (1990).



## 4. Solving the QM Hamiltonian

### 4.1 The Full Electronic Problem

Under the **Born–Oppenheimer approximation**, nuclear and electronic motions are separated. The total QM wavefunction factorizes:

$$\Psi_\text{QM}^\text{tot}(\mathbf{r}^N, \mathbf{R}^M) = \Psi_\text{QM}^\text{nu}(\mathbf{R}^M)\,\Psi_\text{QM}^\text{elec}(\mathbf{r}^N; \mathbf{R}^M)$$

The electronic Schrödinger equation is:

$$\hat{H}_\text{QM}^\text{elec}(\mathbf{r}^N; \mathbf{R}^M)\,\Psi_\text{QM}^\text{elec} = E_\text{QM}^\text{elec}(\mathbf{R}^M)\,\Psi_\text{QM}^\text{elec}$$

with the electronic Hamiltonian (in atomic units):

$$\hat{H}_\text{QM}^\text{elec} = -\frac{1}{2}\sum_{i=1}^{N}\nabla_i^2 - \sum_{i=1}^{N}\sum_{A=1}^{M}\frac{Z_A}{r_{iA}} + \sum_{i=1}^{N}\sum_{j>i}^{N}\frac{1}{r_{ij}}$$

The total QM energy (including nuclear repulsion) is:

$$E_\text{QM}(\mathbf{R}^M) = \langle\Psi_\text{QM}^\text{elec}|\hat{H}_\text{QM}^\text{elec}|\Psi_\text{QM}^\text{elec}\rangle + \sum_{A=1}^{M}\sum_{B>A}^{M}\frac{Z_A Z_B}{R_{AB}}$$

### 4.2 The Slater Determinant and LCAO-MO

The many-electron wavefunction is approximated as a **Slater determinant** of molecular orbitals (MOs):

$$\Psi_\text{QM}^\text{elec} = \frac{1}{\sqrt{N!}}\begin{vmatrix} \phi_1(1) & \phi_2(1) & \cdots & \phi_N(1) \\ \phi_1(2) & \phi_2(2) & \cdots & \phi_N(2) \\ \vdots & \vdots & \ddots & \vdots \\ \phi_1(N) & \phi_2(N) & \cdots & \phi_N(N) \end{vmatrix}$$

Each MO is expanded as a **linear combination of atomic orbitals (LCAO)** — the basis set:

$$\phi_i = \sum_{\mu} C_{i\mu}\,\chi_\mu$$

### 4.3 The SCF Procedure (Hartree–Fock)

Applying the variational principle,

$$E_\text{QM}^\text{trial} \geq E_\text{QM}^\text{exact}, \qquad \frac{\partial E_\text{QM}}{\partial C_{\mu i}} = 0$$

leads to the **Roothaan–Hall equations**:

$$\mathbf{F}_\text{QM}^\text{elec}\,\mathbf{C} = \mathbf{S}\,\mathbf{C}\,\mathbf{E}_\text{QM}^\text{elec}$$

where $\mathbf{F}$ is the Fock matrix, $\mathbf{S}$ is the overlap matrix, $\mathbf{C}$ is the matrix of MO coefficients, and $\mathbf{E}$ contains the MO energies. These are solved iteratively in the **self-consistent field (SCF)** cycle:

```
Input Geometry
    → Compute electron integrals
    → Initial density matrix
    → Build Fock matrix
    → Diagonalize Fock matrix
    → Update density matrix
    → Converged? → Exit
         ↑________________|  (SCF cycle)
```

### 4.4 Computational Scaling

| Method | Scaling |
|---|---|
| MM with PME | $\mathcal{O}(N \log N)$ |
| HF / DFT | $\mathcal{O}(N^3)$–$\mathcal{O}(N^4)$ |
| MP2 | $\mathcal{O}(N^5)$ |
| CCSD(T) | $\mathcal{O}(N^7)$ |



## 5. Schemes for the QM/MM Hamiltonian

Two strategies exist for combining QM and MM energies.

### 5.1 Subtractive Scheme (e.g., ONIOM)

$$E_\text{QM/MM} = E_\text{MM}(\text{whole}) + E_\text{QM}(\text{QM region}) - E_\text{MM}(\text{QM region})$$

The QM region energy at the MM level is subtracted to avoid double counting.

### 5.2 Additive Scheme

$$E_\text{QM/MM} = E_\text{MM}(\text{MM region}) + E_\text{QM}(\text{QM region}) + E_\text{QM/MM}^\text{interaction}$$

Here the MM region explicitly excludes the QM atoms, and the coupling term $E_\text{QM/MM}^\text{interaction}$ is computed explicitly.

> The **additive scheme** is more common in modern biomolecular QM/MM, since it allows the MM environment to polarize the QM density (electrostatic embedding).



## 6. QM/MM Coupling: Embedding Schemes

How the QM/MM interaction is treated determines the accuracy and cost of the calculation.

### 6.1 Mechanical Embedding

- MM point charges do **not** polarize the QM electron density; the QM calculation is effectively done in vacuum.
- QM/MM electrostatics and vdW interactions are evaluated entirely at the MM level.
- Simple but can produce spurious effects (e.g., artificial hydrogen transfer).

### 6.2 Electrostatic Embedding *(most common)*

- MM point charges are included directly in the QM Hamiltonian, polarizing the QM density.
- vdW interactions between QM and MM atoms are still handled at the MM level (require parameterization).
- Used in the vast majority of QM/MM implementations.

### 6.3 Polarization Embedding

- MM atoms carry polarizable dipoles; QM and MM charges are **mutually polarized**.
- Requires a double SCF (micro-iteration): converge both MM dipoles and QM density.
- Highest accuracy, highest cost.

> **References:** Thiel et al., JPC (1996); Morokuma et al., JCTC (2006).



## 7. Electrostatic Embedding in Detail

In the electrostatic embedding scheme the electronic Hamiltonian is augmented with the electrostatic potential of the MM charges $\{q_M\}$ at positions $\mathbf{R}_{MM}$:

$$\hat{H}_\text{QM+QM/MM}^\text{elec} = \hat{H}_\text{QM}^\text{elec} - \sum_{i=1}^{N}\sum_{M=1}^{N_{MM}}\frac{q_M}{r_{iM}}$$

The corresponding Schrödinger equation is:

$$\hat{H}_\text{QM+QM/MM}^\text{elec}(\mathbf{r}^N; \mathbf{R}^M, \mathbf{R}_{MM}^{N_{MM}})\,\Psi_\text{QM}^\text{elec} = E_\text{QM+QM/MM}^\text{elec}(\mathbf{R}^M, \mathbf{R}_{MM}^{N_{MM}})\,\Psi_\text{QM}^\text{elec}$$

The total effective energy is:

$$E_\text{eff} = \langle\Psi_\text{QM}^\text{elec}|\hat{H}_\text{QM+QM/MM}^\text{elec}|\Psi_\text{QM}^\text{elec}\rangle + \sum_{A=1}^{M}\sum_{B>A}^{M}\frac{Z_A Z_B}{R_{AB}} + \sum_{A=1}^{M}\sum_{M=1}^{N_{MM}}\frac{Z_A q_M}{R_{AM}} + E_\text{QM/MM}^\text{vdW} + E_\text{QM/MM}^\text{boundary} + E_\text{MM}$$

The MM charges enter the Fock matrix build step of the SCF cycle, so the QM wavefunction self-consistently responds to the electrostatic environment.

> **Reference:** Field, JCC (1990).



## 8. Practical Issue I: QM Region and Method Selection

### 8.1 Which atoms belong in the QM region?

- At minimum: the substrate(s) and key catalytic residues directly involved in chemistry.
- For metal-containing enzymes: include metals and their first coordination shell.
- Typical QM region size: **< 200 atoms**.

### 8.2 Which QM method and basis set?

| Method | Accuracy | Speed | Parallel efficiency |
|---|---|---|---|
| HF | Poor | Slow $\mathcal{O}(N^4)$ | Reasonable–Good |
| DFT | Good | Slow $\mathcal{O}(N^4)$ | Reasonable–Good |
| MP2 | Very good | Very slow $\mathcal{O}(N^5)$ | Poor |
| Semi-empirical (AM1, PM3, SCC-DFTB) | System-dependent | Very fast ($\sim 10^3\times$ HF/DFT) | Good |

- For HF/DFT: use at least a **6-31G(d)** basis set; check for basis set convergence.
- SE-QM methods can be improved by reparameterization for specific systems.

> **Rule of thumb:** The choice is driven by the required accuracy and total simulation cost. Most long QM/MM MD studies use semi-empirical QM; ab initio QM/MM is increasingly feasible for energy corrections.



## 9. Practical Issue II: QM/MM Boundary Treatment

When the QM/MM boundary cuts through a covalent bond (common in proteins), the dangling valence on the QM side must be capped.

### 9.1 H-Link Atom Approach

- A hydrogen atom is placed along the broken bond to cap the QM region.
- Simple to implement; widely used.
- Care needed: the H-link atom should not interact with nearby MM charges.

### 9.2 Double-Link Atom Approach (Brooks et al.)

- Two link atoms are used (one in MM, one in QM).
- Gaussian smearing applied to MM atoms near the boundary to reduce overpolarization.

### 9.3 Pseudo-Bond Method

- The boundary C–C bond is replaced by a pseudo-atom with optimized parameters.
- Yang and coworkers: DFT-based parameterization.
- Thiel and coworkers: connection bond for SE-QM methods.

### 9.4 Local Self-Consistent Field (LSCF) Method (Rivail et al.)

- A frozen localized orbital spans the QM/MM boundary.
- Difficult to implement; transferability is not guaranteed.

### 9.5 Generalized Hybrid Orbital (GHO) Method (Gao et al.)

- The boundary MM atom carries hybrid orbitals; one "active" orbital enters the QM calculation.
- Primarily designed for SE-QM methods.

### 9.6 Frozen Orbital Approximation (Friesner et al.)

- A pre-optimized frozen orbital caps the QM region.
- Implemented in Schrödinger's QSite/Jaguar.

> **Practical advice:** Always cut along a **C–C single bond**, at least several bonds away from the chemically active region. Be aware of polarization artifacts near the cutting bond, especially with the H-link approach.



## 10. Practical Issue III: Periodic vs. Non-Periodic Boundary

### 10.1 The Problem

Most enzyme QM/MM simulations are solvated in a periodic box. Long-range electrostatics must be handled consistently across QM, QM/MM, and MM regions.

- Most QM/MM implementations historically used a **cutoff** or **no-cutoff** scheme under stochastic boundary conditions — which introduces artifacts.

### 10.2 QM/MM with Particle Mesh Ewald (PME)

A rigorous extension of PME to QM/MM partitions the total energy as:

$$E_\text{tot}^\text{PME} = E_\text{QM}^\text{RS}[\rho] + E_\text{QM/MM}^\text{RS}[\rho] + \Delta E_\text{QM}^\text{PME}[Q] + E_\text{MM}^\text{PME}[q]$$

where RS denotes a real-space short-range contribution and $Q$, $q$ are QM multipole and MM charges, respectively.

**Advantages:**

- Full PBC with no cutoff artifacts.
- Balanced treatment of QM, QM/MM, and MM long-range electrostatics.
- Produces stable MD trajectories.
- Available in CHARMM, AMBER, Q-Chem, and others.

> **References:** Nam, JCTC (2005, 2014). Alternative formulations include ambient-potential composite Ewald (York et al.), multipole moments (Rothlisberger et al.), ESP/ChElPG charges (Herbert et al.), Gen-Ew (Thiel et al.), and augmented charges (Shao et al.).



## 11. Free Energy Methods Combined with QM/MM

Converged free energy calculations require **> 1 ns of QM/MM MD**, which is only feasible with efficient SE-QM/MM methods.

### 11.1 Potential of Mean Force (PMF) by Umbrella Sampling

The reaction free energy along a distinguished coordinate $z$ (e.g., $r_{AH} - r_{HB}$ for a proton transfer):

$$W_{CM}(z) = -RT\ln\rho(z) + C$$

where $\rho(z)$ is the probability density sampled along $z$, and $C$ is a constant. Harmonic bias potentials ("umbrella" windows) are applied to sample the full range of $z$; WHAM or similar methods then reconstruct the unbiased PMF.

### 11.2 Alchemical Free Energy by Thermodynamic Integration

For processes such as solvation free energies or pKa calculations, a coupling parameter $\lambda$ interpolates between states A and B:

$$\Delta G_{A \to B} = \int_{\lambda'=0}^{\lambda'=1}\left.\frac{\partial G}{\partial\lambda}\right|_{\lambda'}\!d\lambda' = \int_{\lambda'=0}^{\lambda'=1}\left\langle\frac{dH(\lambda)}{d\lambda}\right\rangle_{\!\lambda'}\!d\lambda'$$

### 11.3 Path-Based Methods

The finite temperature string method optimizes a reaction pathway $\alpha$ in a collective variable space:

$$F_{\alpha^*} - F_0 = \int_0^{\alpha^*}\sum_{j=1}^{N}\frac{\partial F}{\partial\theta_j^*}\cdot\frac{d\theta_j^*}{d\alpha}\,d\alpha$$



## 12. Speed and Parallelization

### 12.1 Typical Simulation Speeds

| Method | Speed |
|---|---|
| MM (GPU) | $\sim 100$ ns/day |
| MM (CPU) | $\sim 40$ ns/day |
| SE-QM/MM | $\sim 0.2$ ns/day |
| AI-QM/MM | $\sim 0.001$ ns/day |

Since free energy convergence requires > 1 ns, **SE-QM/MM is currently the practical workhorse** for QM/MM free energy calculations.

> **References:** Nam, JCTC (2013, 2014); Ojeda-May, JCTC (2017).

### 12.2 Bottlenecks in SE-QM/MM

The two dominant costs are:

**1. SCF iteration — Fock matrix diagonalization, $\mathcal{O}(N^3)$**

Accelerators to reduce the number of SCF iterations:

- **DIIS** (Direct Inversion in the Iterative Subspace)
- **DXL-BOMD** (Extended Lagrangian Born–Oppenheimer MD with dissipation): an auxiliary density matrix $\mathbf{P}$ is propagated as a dynamical variable:

$$\mathbf{P}(t+\delta t) = 2\mathbf{P}(t) - \mathbf{P}(t-\delta t) + \delta t^2\,\varpi^2\left[\mathbf{D}(t) - \mathbf{P}(t)\right]$$

- **ELMD** (Extended Lagrangian MD): treats density matrix elements as dynamical variables satisfying the idempotency condition $\mathbf{PP} = \mathbf{P}$; highly parallelizable.

ELMD Lagrangian:

$$L = L_\text{QM/MM} + \frac{1}{2}m_P\sum_{i<j}\left(\frac{dP_{ij}}{dt}\right)^2 - \text{Tr}\left[\Lambda(\mathbf{PP}-\mathbf{P})\right]$$

> **References:** Niklasson, JCP (2009); Schlegel et al., JCP (2001); Herbert & Head-Gordon, JCP (2004); Nam, JCTC (2013); Ojeda-May et al., JCTC (2017).

**2. Long-range QM/MM electrostatics — handled via QM/MM-PME** (see Section 10.2).



## 13. Speeding up AI-QM/MM: Multiscale Approach

Ab initio QM/MM (AI-QM/MM) is up to 100,000× slower than MM. A two-level approach decomposes the total energy:

$$E_\text{tot} = E_\text{AI-QM/MM-PME}(\rho:\text{MM}^\text{PBC})$$

$$\approx E_\text{SE-QM/MM-PME}(\rho:\text{MM}^\text{PBC}) + \left[E_\text{AI-QM/MM}(\rho:\text{MM}^\text{RS}) - E_\text{SE-QM/MM}(\rho:\text{MM}^\text{RS})\right]$$

**Strategy:**

- Long-range electrostatics are handled at the inexpensive SE-QM/MM level with full PBC/PME.
- The AI-QM/MM correction (short-range) is applied less frequently using the **rRESPA multiple time step algorithm** — analogous to the ONIOM energy interpolation.

This reduces the AI-QM/MM cost by roughly an order of magnitude without sacrificing accuracy for the chemically important region.

> **Reference:** Nam, JCTC (2014).



## 14. Summary and Practical Recommendations

### QM Region and Method

- QM/MM is **not a black box**: apply chemical knowledge of the reaction mechanism.
- Transition metals are particularly challenging — be cautious with method and region selection.
- Choose between DFT and semi-empirical QM based on the required accuracy and timescale.
- Expect to reparameterize vdW parameters for SE-QM atoms.
- A larger QM region is not always better: check basis set convergence and theory-level limitations.
- Most X-ray structures lack natural substrates — substrate placement via modeling or docking is often required.

### QM/MM Boundary

- Always cut along a **C–C single bond**, at least several bonds away from the reactive site.
- Be aware of polarization artifacts near the cutting bond, especially with the H-link approach.
- Several boundary methods exist; the choice depends on the QM level and available software.

### Periodic Boundary Conditions

- Use **PBC with PME** whenever possible — not all QM packages or theories support this.
- If PBC is unavailable, use **no-cutoff** for both QM/MM and MM interactions.
- Do **not** mix cutoff for MM with no-cutoff for QM/MM interactions.





## Further Reading

* Warshel & Levitt (1976) – Original QM/MM formulation
* Field, Bash & Karplus (1990) – Modern QM/MM implementation
* Senn & Thiel (2009) – QM/MM methods review
* van der Kamp & Mulholland (2013) – Biomolecular QM/MM applications
