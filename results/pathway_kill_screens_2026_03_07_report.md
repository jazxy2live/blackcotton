# Pathway Kill Screens

These screens answer three kill questions using the current BlackCotton model stack.

## Reference

- Config: `active production config`
- Reference candidate source: `results/final_lab_top3.json candidate 2`
- Reference candidate: mat/scw `38/34`, k `0.10`, eff `1.60`, ret `0.75`

## Simulation 1: Cotton Toxicity Threshold

- Kill metric: ROS proxy exceeds antioxidant capacity before 90% of simulated final pigment load.
- Current pathway status: `PASS`
- `balanced_current_pathway`: reference ROS ratio `0.549`, safe through melA `2.125`, first kill `2.25`, kill melanin `0.007644669340014237`
- `support_starved_pathway`: reference ROS ratio `0.604`, safe through melA `1.875`, first kill `2.0`, kill melanin `0.0006222641072688878`

## Simulation 2: CSC Metamaterial Mutation

- Kill metric: No morphology that reaches total absorption keeps tensile strength above 28 g/tex.
- Structural pathway status: `DEAD`
- Best absorber found: absorption `0.978`, strength `8.04 g/tex`
- Best safe absorber: absorption `0.225`, strength `28.23 g/tex`
- Best strength at absorption >= `0.90`: `15.868628593686347`
- Best strength at absorption >= `0.95`: `12.40519385384609`
- Best strength at absorption >= `0.97`: `9.423311844410117`
- Best strength at absorption >= `0.99`: `None`

## Simulation 3: Bacterial Co-Culture Alternative

- Kill metric: If hydrogen-bond retention falls below 0.80 at the target melanin load, the architecture needs binding-protein rescue.
- Target melanin/cellulose ratio: `0.24`
- `co_culture` no-linker: affinity `0.252`, hbond `0.742`, wet strength `15.85`
- `co_culture` rescue requirement: binding protein factor `0.3`
- `single_strain` no-linker: affinity `0.320`, hbond `0.761`, wet strength `16.26`
- `single_strain` rescue requirement: binding protein factor `0.2`

## Notes

- Simulation 1 is grounded in the existing ODE stack.
- Simulations 2 and 3 are coarse-grained surrogate screens, not atomistic MD or full fermentation process models.
