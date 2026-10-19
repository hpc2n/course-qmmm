# From Classical Force Fields to Machine-Learning Potentials

## Introduction

Molecular dynamics simulations of minerals, mineral–water interfaces, Biomolecules, and electrolyte solutions have traditionally relied on classical force fields. 
These force fields have been extremely successful because they allow us to simulate systems containing hundreds of thousands or even millions of 
atoms over relatively long time scales.

For a classical force-field developer, the potential energy is normally written as a sum of physically motivated analytical contributions, schematically,

$$
E_{\mathrm{FF}}=E_{\mathrm{bonded}}+E_{\mathrm{vdW}}+E_{\mathrm{Coulomb}}.
$$

The electrostatic contribution, for example, is commonly represented by Coulomb interactions,

$$
E_{\mathrm{Coulomb}}=\sum_{i<j}\frac{q_iq_j}{4\pi\epsilon_0 r_{ij}},
$$

For periodic ionic systems, this pair expression requires an appropriate periodic electrostatic summation rather than a simple finite sum.

Short-range dispersion and repulsion can be represented by a Lennard-Jones potential,

$$
E_{\mathrm{LJ}}=\sum_{i<j}4\epsilon_{ij}\left[\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12}-\left(\frac{\sigma_{ij}}{r_{ij}}\right)^6\right].
$$

For mineral systems, one may additionally introduce harmonic bonds, angle terms, Buckingham potentials, polarization models, three-body terms, or other 
specialized interactions.

The important characteristic of this approach is that **we choose the mathematical form of the potential before fitting it**. Force-field development 
then consists largely of determining appropriate parameters such as charges, equilibrium distances, force constants, and van der Waals parameters.

Machine-learning potentials change this philosophy. Instead of asking:

> What analytical potential should I use to represent an Fe–O or Si–O interaction?

one can ask:

> Can we learn the potential-energy surface directly from electronic-structure calculations?

This is the central idea behind machine-learning interatomic potentials.

## From analytical force fields to learned potential-energy surfaces

Consider a system containing $N$ atoms. One can describe the configuration using the atomic coordinates
$\mathbf R =\{\mathbf r_1,\mathbf r_2,\ldots,\mathbf r_N\},$ the atomic numbers $\mathbf Z =\{Z_1,Z_2,\ldots,Z_N\},$
and, for a periodic system, the simulation-cell matrix $\mathbf H.$

An ML potential attempts to construct a function $E_\theta=E_\theta(\mathbf R,\mathbf Z,\mathbf H),$
where $\theta$ represents the parameters of the neural network or another machine-learning model.

We train these parameters using quantum-mechanical reference calculations. Ideally,

$$
E_\theta(\mathbf R)\approx E_{\mathrm{DFT}}(\mathbf R).
$$

But energies alone are generally not sufficient for molecular dynamics. We also want the potential-energy 
surface to have the correct derivatives. The force acting on atom $i$ is

$$
\mathbf F_i=-\frac{\partial E}{\partial\mathbf r_i}.
$$

Therefore, a useful ML potential should satisfy

$$
-\frac{\partial E_\theta}{\partial\mathbf r_i}\approx \mathbf F_i^{\mathrm{DFT}}.
$$

This is particularly important for molecular dynamics because the forces, rather than the absolute energies, 
determine the evolution of the atomic coordinates. In this sense, an ML potential is still a force field. 
What changes is the representation of the potential-energy surface.

## Locality and atomic environments

Many modern ML potentials make an important approximation: the total energy can be decomposed into atomic contributions,

$$
E_{\mathrm{ML}}=\sum_{i=1}^{N} E_i.
$$

The energy assigned to atom $i$ depends on its local atomic environment,

$$
E_i =f_\theta(\mathcal N_i),
$$

where the neighborhood can be defined as

$$
\mathcal N_i=\left\{j\ne i: r_{ij}<r_c\right\}.
$$

Here $r_c$ is the cutoff radius. The model uses the neighbors' species and relative positions, not just their indices; periodic images must also be 
included where appropriate. Message-passing models can have an effective receptive field extending over several neighbor hops.

This construction has several useful properties. It makes the computational cost approximately linear with system size at fixed density, cutoff, and model 
complexity, allows models trained on relatively small configurations to be applied to larger systems, and naturally represents local chemical environments.

For a mineral, the network may therefore encounter environments such as

$$
\mathrm{Si-O-Si},
\qquad
\mathrm{Si-O-Al},
\qquad
\mathrm{Fe-O},
\qquad
\mathrm{Al-OH},
\qquad
\mathrm{Na-O_{water}},
$$

without us explicitly writing separate analytical terms for each one.

The network learns how those environments contribute to the energy from the training data.

This does **not**, however, mean that all physics is automatically learned correctly. Long-range electrostatics, for example, deserves special attention in 
mineral systems, as we will discuss later.

## ANI as an introduction to ML potentials

ANI provides a relatively intuitive example of an ML potential.

ANI represents the environment around each atom using descriptors traditionally called atomic environment vectors. These descriptors encode radial and angular
information about neighboring atoms. Conceptually, we can write

$$
E =\sum_i NN_{Z_i}(G_i),
$$

where $G_i$ describes the environment around atom $i$, and $NN_{Z_i}$ is a neural network associated with its chemical element.

ANI is particularly convenient as an introduction to GROMACS machine-learning potentials because TorchANI is based on PyTorch, and ANI models can be packaged as
TorchScript models.

For example, ANI-2x can be wrapped so that GROMACS provides coordinates and atomic numbers, the wrapper converts the units expected by ANI, and the network 
returns an energy.

However, ANI-2x has an extremely important limitation for our application.

Its supported elements are $\mathrm{H,\ C,\ N,\ O,\ F,\ S,\ Cl}$; see the [TorchANI element definitions][torchani].

This chemical space is appropriate for many organic molecules, but not for typical mineral systems.

Consider, for example, a hydrated mineral containing $\mathrm{Fe,\ Si,\ Al,\ O,\ H,\ Na}.$

ANI-2x recognizes O and H but does not provide models for Fe, Si, Al, or Na.

Therefore, ANI-2x cannot simply be applied to the complete mineral–water–sodium system.

This illustrates an important distinction:

> The capabilities of an ML architecture and the capabilities of a particular pretrained model are not the same thing.

An architecture might theoretically represent many elements, but a pretrained model is only meaningful within the chemical space represented during its training.

## MACE for mineral systems

For materials simulations, architectures such as MACE offer a flexible way to represent many-body interactions.

MACE stands for **Multi-Atomic Cluster Expansion** and combines ideas from atomic cluster expansions with equivariant message-passing neural networks.

The word *equivariant* is important. Physical energies should be invariant under rotation,

$$
E(\mathbf R)=E(\mathcal R\mathbf R),
$$

where $\mathcal R$ is a rotation. For periodic systems, the cell must be rotated consistently; any external fields must also be transformed.

Forces, however, are vectors and must rotate with the system,

$$
\mathbf F_i(\mathcal R\mathbf R)=\mathcal R\mathbf F_i(\mathbf R).
$$

Modern equivariant neural networks build these symmetry properties into the architecture rather than requiring the network to learn them purely from data.

MACE constructs increasingly sophisticated representations of the environment surrounding each atom using information exchanged between neighboring atoms. 
Very schematically, one can think of an atomic energy as containing contributions associated with increasingly complex many-body environments,

$$
E_i\sim E_i^{(0)}+\sum_j E_{ij}^{(1)}+\sum_{jk}E_{ijk}^{(2)}+\cdots.
$$

The actual MACE formulation is considerably more sophisticated than this expression, but the important physical idea is that the network can represent many-body correlations.

For minerals this is attractive because the energy of an oxygen atom, for example, does not depend merely on the distance to one silicon atom. Its environment 
can depend on coordination geometry, neighboring Al or Fe substitutions, hydroxylation, water molecules, ions, and many other factors.

## Training our own mineral ML potential

Suppose our target simulation contains

$$
\boxed{\mathrm{Fe,\ Si,\ Al,\ O,\ H,\ Na}}.
$$

Rather than trying to force ANI-2x to describe this chemical space, we could train a model such as MACE specifically for it.

The first task is not actually neural-network training.

The first task is **constructing the reference dataset**.

We need atomic configurations representative of the states that the eventual molecular-dynamics simulation will explore. For a mineral–water system, that might 
include:

- bulk mineral structures;
- relaxed and distorted lattices;
- different Al/Si substitutions;
- Fe-containing environments;
- hydroxylated surfaces;
- mineral–water interfaces;
- different water coverages;
- solvated Na ions;
- Na approaching the mineral surface;
- surface-bound Na;
- defects;
- strained structures;
- thermally distorted configurations.

For each configuration, we perform a consistent electronic-structure calculation, typically using DFT, obtaining

$$
\mathbf R \longrightarrow \left\{E_{\mathrm{DFT}},\mathbf F_{\mathrm{DFT}},\boldsymbol{\sigma}_{\mathrm{DFT}}\right\},
$$

where $\boldsymbol{\sigma}$ represents stress information when required.

For Fe-containing systems, the reference calculations must also treat relevant spin states and oxidation environments consistently. A coordinate-only model cannot
distinguish multiple electronic states of an identical geometry without additional state information.

The same electronic-structure methodology should normally be used throughout the dataset. Otherwise the ML potential is being asked to reproduce an inconsistent 
potential-energy surface.

## Training is still force-field fitting

For classical force-field developers, it is useful to think of ML training as another parameterization procedure.

Instead of optimizing parameters such as

$$
q_i,\quad \sigma_{ij},\quad \epsilon_{ij},\quad k_b,
$$

we optimize thousands or millions of neural-network parameters,

$$
\theta = \{\theta_1,\theta_2,\ldots,\theta_M\}.
$$

A typical loss function might combine errors in energy, forces, and stress:

$$
\mathcal L = w_E\mathcal L_E +w_F\mathcal L_F +w_\sigma\mathcal L_\sigma.
$$

For example,

$$
\mathcal L_E =\left(E_{\mathrm{ML}}-E_{\mathrm{DFT}}\right)^2,
$$

while a force loss could be

$$
\mathcal L_F =\frac{1}{3N}\sum_{i=1}^{N}\left|\mathbf F_i^{\mathrm{ML}}-\mathbf F_i^{\mathrm{DFT}} \right|^2.
$$

For molecular dynamics, force accuracy is particularly important because even relatively small systematic errors in forces can affect structure and dynamics.

But there is a very important warning here:

$$
\boxed{\text{Low test-set RMSE does not necessarily mean a good force field.}}
$$

A potential could reproduce randomly selected configurations from the same dataset extremely well while failing catastrophically when MD encounters a new atomic environment.

## The extrapolation problem

Suppose we train exclusively using equilibrium configurations of a mineral at 300 K.

The resulting model may be extremely accurate close to those configurations.

But during an MD simulation an atom might enter an unusual geometry because of thermal fluctuations, surface rearrangement, ion adsorption, or simply numerical noise.

If the model has never seen anything remotely similar, its prediction is an extrapolation.

Unlike many classical potentials, a neural network does not necessarily become strongly repulsive simply because two atoms approach each other too closely unless this behavior is represented by its architecture or training data.

The network can therefore produce an unphysical low-energy configuration or enormous forces.

A training set intended to support stable dynamics should consequently contain more than equilibrium structures. We want something closer to

$$
\text{training space}=\text{equilibrium}+\text{thermal distortions}+\text{strained configurations}+\text{interfaces}+\text{unusual local environments}.
$$

This is where **active learning** becomes particularly useful.

We can train an initial potential, perform MD, use an uncertainty indicator such as ensemble disagreement to select configurations, calculate those structures with DFT, add them to the training set, and retrain.

The cycle becomes

$$
\mathrm{DFT} \rightarrow \mathrm{training}\rightarrow \mathrm{MD} \rightarrow \mathrm{uncertain\ configurations}\rightarrow \mathrm{DFT}\rightarrow \cdots
$$

This is analogous to discovering during classical force-field development that a particular structural environment is poorly represented and adding new reference data to constrain it.

## Connecting an ML potential to GROMACS

The next question is how we actually use such a model during molecular dynamics.

GROMACS 2026 provides the **NNPot** interface for neural-network potentials and requires a compatible LibTorch build. See the [GROMACS NNPot documentation][gromacs-nnpot].

Conceptually, during each MD step GROMACS has the coordinates $\mathbf R(t).$

These coordinates, together with information such as atomic numbers and the periodic cell, are passed to a PyTorch/TorchScript model:

$$
\{\mathbf R,\mathbf Z,\mathbf H,\mathrm{PBC}\}\longrightarrow E_{\mathrm{ML}}.
$$

The force is then obtained from the energy gradient,

$$
\mathbf F_i=-\frac{\partial E_{\mathrm{ML}}}{\partial\mathbf r_i}.
$$

GROMACS can then integrate Newton's equations,

$$
m_i\frac{d^2\mathbf r_i}{dt^2}=\mathbf F_i,
$$

just as it would with a classical potential.

The integrator itself does not need to know whether the forces came from a Lennard-Jones expression, PME electrostatics, or a neural network.

## TorchScript as the interface

A convenient feature of the GROMACS implementation is that the ML model can be provided as a **TorchScript model**.

For ANI, for example, we can construct a wrapper that accepts GROMACS inputs:

```python
def forward(
    self,
    positions,
    atomic_numbers,
    box=None,
    pbc=None
):
    ...  # Interface sketch; model evaluation is omitted.
```

One important responsibility of this wrapper is **unit conversion**.

GROMACS coordinates are normally expressed in nanometers, whereas many ML chemistry models operate in Ångström. Therefore,

$$
1\ \mathrm{nm}=10\ \mathring{\text{A}}.
$$

For ANI, energies are naturally associated with Hartree, while GROMACS uses

$$
\mathrm{kJ\,mol^{-1}}.
$$

Thus we use approximately

$$
1\ E_h =2625.5\ \mathrm{kJ\,mol^{-1}}.
$$

The wrapper therefore performs transformations such as

$$
\mathbf R_{\mathring{\text{A}}}=10\mathbf R_{\mathrm{nm}},
$$

converts the periodic cell by the same length factor, evaluates the network, and converts its energy back to GROMACS units. Keeping these operations in the differentiable computation ensures that energy gradients carry the correct force conversion. If forces are returned explicitly, their units must be converted 
separately.

Finally,

```python
torch.jit.script(model).save("model.pt")
```

exports a scriptable wrapper. Its input order and output format must match NNPot, and the PyTorch/LibTorch versions must be compatible.

For a model trained in electronvolts, as is common in materials ML, the relevant energy conversion is instead approximately

$$
1\ \mathrm{eV}
=
96.4853\ \mathrm{kJ\,mol^{-1}}.
$$

Unit consistency is absolutely critical because an unnoticed factor of ten in distance or a factor of 96 in energy will obviously produce meaningless dynamics.

## Neighbor lists and GROMACS 2026

Local ML potentials require neighbors satisfying $r_{ij}<r_c.$

But GROMACS already has highly optimized machinery for constructing neighbor lists.

This raises an obvious question:

> Why should an ML potential independently search the complete system for neighboring atoms?

The GROMACS 2026 NNPot interface can provide atom-pair information and corresponding periodic shifts to the model; see the [model-input options][gromacs-mdp]. 
This is particularly relevant for graph/message-passing architectures such as MACE.

Conceptually,

$$
\text{GROMACS neighbor search}
\longrightarrow
\text{pair list}
\longrightarrow
\text{ML model}.
$$

This is an important development because efficient neighbor construction is essential for scaling ML potentials to realistic condensed-phase systems.

## Machine-learning potentials in LAMMPS

LAMMPS provides another route for using learned interatomic potentials in mineral simulations. The scientific workflow is the same: prepare reference data, fit or
fine-tune a model, validate it, and then use the model to calculate forces during MD. The software interface and model export format differ from those used by GROMACS.

### Model families and interfaces

LAMMPS exposes many ML potentials through `pair_style` commands. Despite the name, a pair style can implement a many-body potential; it does not imply that the 
energy is a sum of independent pair interactions.

| Model or interface | LAMMPS pair style | Main role |
| --- | --- | --- |
| SNAP | `snap` | Uses bispectrum descriptors with fitted linear or quadratic energy models; requires the ML-SNAP package. |
| Atomic cluster expansion (ACE) | `pace` | Evaluates supported ACE potentials; requires the ML-PACE package. |
| High-dimensional neural-network potentials | `hdnnp` | Connects to n2p2 models; requires the ML-HDNNP package and its library dependency. |
| General ML interface | `mliap` | Supports model/descriptor combinations or a unified Python model; requires ML-IAP and dependencies appropriate to the selected model. |

These are examples, not an exhaustive list. Each interface expects its own model files and build configuration. A `.pt` extension alone does not establish compatibility. See the official documentation for [SNAP][snap], [PACE][pace], [HDNNP][hdnnp], and [ML-IAP][mliap].

### Example: a MACE model through ML-IAP

MACE documents both an original LAMMPS interface and a newer ML-IAP route. For the latter, its guide describes exporting a trained model with `create_lammps_model.py` and the `--format=mliap` option. The documented GPU workflow requires compatible MACE, Python, CUDA, cuEquivariance, and LAMMPS/Kokkos installations. Follow the [MACE installation and export guide][mace-lammps] for the versions being used.

For a model trained for Fe, Si, Al, O, H, and Na, an illustrative input fragment is:

```lammps
# Requires a compatible MACE/ML-IAP installation and exported model.
units          metal
atom_style     atomic
boundary       p p p
newton         on

# The data file must define the box, coordinates, types, and masses.
read_data      mineral_water.data

pair_style     mliap unified mineral.model-mliap_lammps.pt 0
pair_coeff     * * Fe Si Al O H Na

# 0.0005 ps = 0.5 fs; validate this timestep for the actual system.
timestep       0.0005
thermo         100
thermo_style   custom step temp pe ke etotal press

# Evaluate the starting configuration before attempting dynamics.
run            0
```

This fragment assumes a fully periodic system and six atom types ordered as Fe, Si, Al, O, H, and Na. The `pair_coeff` element list maps **LAMMPS type IDs**, in
numerical order, to chemical species. Change it to match the data file; an incorrect mapping can assign the wrong chemistry without changing the coordinates. The 
final `0` on the `pair_style` line is the unified interface's ghost-neighbor flag, following the MACE example. It should not be changed without consulting the 
interface requirements.

This is a setup and single-point check, not a complete production script. It needs the actual model and structure files. Subsequent dynamics also require an 
integration fix and appropriate initial velocities or a restart state. The `atomic` atom style suits this full-ML example; explicitly charged or bonded hybrid 
models may require a different atom style.

### Units and validation

With `units metal`, LAMMPS uses ångströms for distance, electronvolts for energy, picoseconds for time, eV/Å for force, and bar for pressure. With `units real`, 
energy is expressed in kcal/mol and time in femtoseconds. The chosen unit system must agree with the model interface; do not assume that switching `units` 
automatically converts an arbitrary ML model. See the [LAMMPS unit definitions][units].

Before running long trajectories, compare the same configurations in LAMMPS and the model's original evaluation environment. Check total energies, forces, and, 
when required, stresses or virials, accounting for sign, volume, and unit conventions. Then test timestep convergence, energy conservation in NVE, and the 
structural properties relevant to the mineral–water system.

LAMMPS integration does not extend the model's training domain. It also does not supply missing long-range physics automatically. Adding a Coulomb term to a model 
fitted to total DFT energies can double-count interactions. A combined ML/classical model therefore needs a deliberate energy decomposition and consistent training 
targets, followed by validation of the coupled system.

## Full ML versus hybrid simulations

There are two conceptually different ways we might use an ML potential.

The first is to describe the complete system with ML:

$$
E_{\mathrm{total}}=E_{\mathrm{ML}}.
$$

For example, if we had trained a model covering $\mathrm{Fe,\ Si,\ Al,\ O,\ H,\ Na}$, we could potentially apply it to the entire mineral–water–ion system.

For a wrapper accepting positions, atomic numbers, cell, and periodicity in that order, the relevant GROMACS options are:

```ini
nnpot-active       = true
nnpot-modelfile    = mineral.pt
nnpot-input-group  = System
nnpot-model-input1 = atom-positions
nnpot-model-input2 = atom-numbers
nnpot-model-input3 = box
nnpot-model-input4 = pbc
```

These options are only the NNPot portion of the input. The topology, model, and remaining MD settings must also be consistent.

But replacing an established mineral force field completely is not necessarily the only interesting application.

A second possibility is a **hybrid ML/classical approach**.

For example, we might retain the classical description of a large mineral substrate while using ML for a chemically interesting region around a surface, defect, 
adsorbed species, or reactive site.

Conceptually,

$$
E_{\mathrm{total}}=E_{\mathrm{ML}}+E_{\mathrm{MM}}+E_{\mathrm{ML-MM}}.
$$

GROMACS documents additive mechanical embedding, with classical coupling interactions and link atoms where needed. Such coupling must be validated for the chosen 
boundary and chemistry; selecting a region does not by itself guarantee a physically accurate interface.

This possibility is particularly interesting for developers of classical mineral force fields because ML does not necessarily have to replace decades of classical 
force-field development. It may instead complement it.

## The particular problem of electrostatics in minerals

There is one issue that deserves special attention in mineral systems: **long-range electrostatics**.

A local ML potential assumes approximately that $E_i =E_i(\mathcal N_i;r_c).$

For a strictly local descriptor, atoms outside the cutoff do not explicitly enter that environment. Message passing can enlarge the effective range, but a finite 
receptive field still does not provide a general long-range electrostatic treatment.

Coulomb interactions, however, behave as

$$
E_{\mathrm{Coulomb}}
\propto
\frac{q_iq_j}{r_{ij}},
$$

and therefore do not suddenly disappear at 5 or 6 Å.

This is particularly relevant for systems containing:

- charged mineral surfaces;
- isomorphic substitutions;
- counterions;
- electrolyte solutions;
- electric double layers;
- protonated and deprotonated sites.

A local model can implicitly reproduce electrostatic effects represented within its training configurations, but this does not automatically mean that it has learned 
the correct long-range behavior under different conditions.

For example, suppose we train only at one surface charge and one Na concentration. We should not automatically assume that the resulting local model will correctly predict behavior at a very different ionic strength.

This is an area where classical force-field expertise becomes especially valuable. Classical models already have well-developed treatments for long-range electrostatics,
including Ewald summation and PME.

Some modern ML architectures therefore combine learned short-range interactions with explicit electrostatic or charge models.

For mineral applications, deciding how to represent electrostatics may be just as important as choosing between MACE, ANI, NequIP, or another neural-network 
architecture.

## Validation: think like a force-field developer

One of the most important messages is that ML potential validation should not stop with a machine-learning metric.

Suppose someone reports

$$
\mathrm{RMSE}_{F}=30\ \mathrm{meV/\mathring{\text{A}}}.
$$

That number is useful, but it does not answer whether the potential is appropriate for a mineral simulation.

For the bulk mineral, we might examine lattice parameters,

$$
a,\ b,\ c,\ \alpha,\beta,\gamma,
$$

as well as density, elastic properties, structural distributions, defect energies and vibrational properties.

For the mineral–water interface, we can examine density profiles and radial distribution functions such as

$$
g_{\mathrm{O_{surface}-O_{water}}}(r).
$$

We should test water orientation, adsorption energies, hydration structures and surface relaxation.

For sodium we might examine

$$
g_{\mathrm{Na-O_{water}}}(r),
$$

coordination numbers, hydration-shell structure, residence times, diffusion, and adsorption at surface sites.
If chemical reactions are expected, such as proton transfer, we need explicit validation of those reactions.
We should also perform stability tests.

An NVE simulation is particularly informative. In the absence of numerical and model problems,

$$
E_{\mathrm{total}}=E_{\mathrm{kinetic}}+E_{\mathrm{potential}}
$$

should be approximately conserved, and therefore

$$
\frac{dE_{\mathrm{total}}}{dt}\approx 0.
$$

Significant systematic drift can indicate problems with the integration timestep, model precision, neighbor treatment, or the consistency of the predicted energy 
and forces.

## Classical force fields and ML potentials are complementary

It is tempting to describe ML potentials as replacements for classical force fields, but I think that is the wrong way to view the problem, particularly for 
mineral simulations.

Classical mineral force fields embody substantial physical and chemical knowledge:

- how electrostatics should be represented;
- which structural observables matter;
- what transferability means;
- which mineral environments are difficult;
- how substitutions affect local structure;
- how water behaves at mineral surfaces;
- what experimental data provide meaningful validation.

Machine learning does not make this expertise obsolete.

Instead, it changes where some of the modeling decisions occur.

In a traditional force field we spend considerable effort designing and parameterizing quantities such as

$$
q_i,\quad \sigma_{ij},\quad \epsilon_{ij},\quad k_b,\quad r_0.
$$

With an ML potential, some of that effort moves toward designing

$$
\boxed{
\text{reference data}
+
\text{model architecture}
+
\text{training strategy}
+
\text{validation domain}
}.
$$

The training dataset effectively becomes part of the definition of the force field.

## A practical workflow for mineral systems

A realistic development workflow could begin with an existing classical mineral force field.

We first use that force field, experimental structures and chemical intuition to identify the relevant configuration space:

$$
\text{bulk}+\text{surface}+
\text{water}+\text{ions}+\text{defects}.
$$

Representative structures are evaluated with DFT to obtain energies and forces,

$$
\{\mathbf R_k\}\xrightarrow{\mathrm{DFT}}\{E_k,\mathbf F_k,\boldsymbol{\sigma}_k \}.
$$

We then divide these data into training, validation and genuinely independent test sets and train an ML potential such as MACE.
The resulting model is first validated independently of the MD engine against electronic-structure calculations.
If satisfactory, it is exported into a form compatible with GROMACS NNPot or the chosen LAMMPS interface and used for short MD simulations.
Those simulations may reveal configurations that were poorly represented in the original dataset.
We identify these configurations, calculate new DFT reference values, add them to the dataset and retrain.

The development cycle therefore becomes

$$
\boxed{\mathrm{DFT}\rightarrow \mathrm{ML\ training}\rightarrow \mathrm{validation}\rightarrow
\mathrm{MD}\rightarrow \mathrm{new\ configurations}\rightarrow \mathrm{DFT}}.
$$

Eventually, when the potential is sufficiently robust, we can move toward longer production simulations.

This is not fundamentally different from iterative classical force-field development. What has changed is the mathematical representation of the potential and 
the amount of reference data we can incorporate.

## Final perspective

ANI provides an excellent demonstration of how a neural-network potential can communicate with GROMACS. It is relatively straightforward to wrap a TorchANI model, 
convert units, export it through TorchScript, and allow GROMACS to obtain energies and forces.
But ANI-2x is not a general mineral potential because its elemental coverage does not include elements such as Si, Al, Fe or Na.
For mineral applications, models such as MACE are much more interesting because the architecture can be 
trained specifically for the chemistry we need.

If our system contains $\mathrm{Fe,\ Si,\ Al,\ O,\ H,\ Na},$

But the fundamental scientific question has not changed.

For classical force fields we ask:

$$
\boxed{
\text{Does my analytical potential reproduce the physics I care about?}
}
$$

For machine-learning potentials we must ask exactly the same thing:

$$
\boxed{
\text{Does my learned potential reproduce the physics I care about?}
}
$$

A small force RMSE is not the final answer. A sophisticated neural-network architecture is not the final answer. And simply being able to run the model inside GROMACS 
or LAMMPS is not the final answer.

For mineral simulations, successful model development combines classical force-field expertise, representative reference data, a suitable ML architecture, and rigorous 
physical validation.

## Software documentation

- [GROMACS 2026 NNPot][gromacs-nnpot]
- [GROMACS 2026 molecular-dynamics parameters][gromacs-mdp]
- [LAMMPS ML-IAP][mliap]
- [MACE in LAMMPS with ML-IAP][mace-lammps]
- [LAMMPS SNAP][snap], [PACE][pace], and [HDNNP][hdnnp]
- [LAMMPS units][units]

[gromacs-nnpot]: https://manual.gromacs.org/documentation/2026.0/reference-manual/special/nnpot.html
[gromacs-mdp]: https://manual.gromacs.org/documentation/2026.0/user-guide/mdp-options.html
[mliap]: https://docs.lammps.org/pair_mliap.html
[mace-lammps]: https://mace-docs.readthedocs.io/en/latest/guide/lammps_mliap.html
[snap]: https://docs.lammps.org/pair_snap.html
[pace]: https://docs.lammps.org/pair_pace.html
[hdnnp]: https://docs.lammps.org/pair_hdnnp.html
[units]: https://docs.lammps.org/units.html

[torchani]: https://aiqm.github.io/torchani/api_autogen/torchani.utils.html
