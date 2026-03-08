# Plant vs Bacterial Decision Report

Decision checkpoint for the BlackCotton primary path and fallback path.

## Recommendation

- Primary path: `plant_chemical`
- Secondary path: `bacterial_single_strain_binder_engineering`
- Deprioritized path: `csc_structural_metamaterial`
- Reason: Plant path still clears cotton-like strength with an active chemical gate, while the bacterial fallback remains materially below cotton spinning strength even after binder rescue.

## Plant Path

- Reference candidate: `mat/scw 38/34`, `k 0.10`, `eff 1.60`, `ret 0.75`
- Reference performance: success `0.858`, robust `0.864`, strength `29.00 g/tex`, yield `0.902`
- Chemical gate: `116` / `116` refined survivors remain chemically safe
- Chemical ROS: best refined ratio `0.519`, reference ratio `0.549`, first forced-kill point `melA 2.25x`
- Worst correlated scenario: `strict_high_noise_correlated` with best min success `0.362`
- Current search frontier: strength up to `29.13 g/tex`, yield up to `0.913`, ROS ratio down to `0.519`
- Chemical-gated hardening: strict+high-noise top success `0.445` -> `0.827` (delta `+0.382`)
- Hardening target status: `False` against target `0.850`
- Chemical-safe seed retention during hardening: `260` / `260`
- Best hardened candidate: `mat/scw 39/34`, `k 0.03`, `eff 1.70`, `ret 0.85`

## Bacterial Fallback

- Architecture: `single_strain`
- Required binder at target ratio `0.24`: `0.20`
- Best focused-window design: ratio `0.30`, binder `0.40`, hbond `0.805`, wet strength `17.08`
- Best safe strength seen in focused window: `17.78`

## Decision

- Keep the plant chemical path as the main program.
- Keep the bacterial path as a backup branch focused on binder engineering, not as the lead path.
- Keep the CSC structural path deprioritized under the current screen.

## Notes

- Plant chemical safety is tied to the current ODE-based ROS proxy and enforced as a hard gate.
- Bacterial results are coarse-grained surrogate outputs, not fermentation-scale wet-lab predictions.
