# TP413_CabEnv_Rosenbauer_Panther_6x6
Project for the course TP413 with Col. Jansens and Mr. Degreef. A quantative analysis on the cabin environment conditions needed in the airport firefighting vehicule Panther 6x6 constructed by Rosenbauer.

## What the code does

It computes, for one cooling case and one heating case, the steady-state energy balance of the
crew cab of the Panther 6x6, with the method of practical work 4 (PW4):

1. energy balance of every wall and window (sun, outside convection, infrared radiation,
   conduction, inside convection);
2. heat and water vapour released by the occupants (method of PW5);
3. energy and water balance of the cabin air → required temperature of the ventilation air;
4. if that temperature is not acceptable: fixed supply temperature, larger flow, recirculation;
5. power of the heating / cooling heat exchanger.

## Files

| File | Content |
|---|---|
| `panther_6x6.jl` | **The file to run.** Vehicle data (geometry, walls, glass), occupants, the two design cases, and the printing of the results. All inputs are here. |
| `cabin_model.jl` | The physics: moist air, convection correlations, surface balance, occupants, ventilation. Contains no vehicle data. |
| `validation_pw4.jl` | Runs the same model on PW4 Ex 4.2 and compares with the worked solution. |
| `Project.toml` | Julia packages needed (CoolProp, Printf). |

## How to run

```
julia --project=. -e "using Pkg; Pkg.instantiate()"    # once, installs CoolProp
julia --project=. panther_6x6.jl
julia --project=. validation_pw4.jl
```

In VS Code: open `panther_6x6.jl` and use "Julia: Execute active File in REPL".

## How the model works

Sign convention: a heat flow is positive when it goes **towards the cabin**. Units are SI
(temperatures in kelvin inside the code, typed and printed in °C).

**One surface** (`solve_surface`). Unknowns: the outside face temperature `T_e` and the inside
face temperature `T_si`. The cabin temperature is imposed, so every surface is solved on its own.

```
outside node :  Q_sun + Q_conv_o + Q_rad = Q_k
inside node  :  Q_k = Q_i

Q_sun    = α·A·G                                   G = DNI·cos θ + G_dif·F_sky
Q_conv_o = h_o·A·(T_amb − T_e)
Q_rad    = ε·A·σ·[F_sky·(T_sky⁴ − T_e⁴) + F_sur·(T_sur⁴ − T_e⁴)]
Q_k      = (T_e − T_si)/R_k                        R_k = Σ d/(k·A)
Q_i      = h_i·A·(T_si − T_cab)
glass only:  Q_trans = τ(θ)·A·G                    sunlight that goes straight into the cabin
```

`h_o` and `h_i` depend on the surface temperatures, so the solution is iterative: compute the
coefficients with the current temperatures, solve the balance, update, repeat. The printed
`residual` is what is left of the two balances at the end (it must be ≈ 0).

**Convection coefficients** (`natural_convection`, `forced_convection`, `outside_convection`).

| Situation | Correlation | Length |
|---|---|---|
| Vertical or inclined wall, natural | Churchill–Chu, AHTT eq. 8.13b (g → g·sin(tilt)) | height |
| Horizontal, unstable (hot side up / cold side down) | AHTT eq. 8.35 (Ra ≥ 10⁷), 8.37a below | area / perimeter |
| Horizontal, stable (hot side down / cold side up) | AHTT eq. 8.38 | short side |
| Forced, flat plate, Re < 5·10⁵ | St = 0.664·Re^(−1/2)·Pr^(−2/3) | along the flow |
| Forced, flat plate, Re ≥ 5·10⁵ | St = (0.037·Re^(−1/5) − 871/Re)·Pr^(−2/3) | along the flow |

Natural or forced outside is decided with the Richardson number Ri = Gr/Re²
(Ri > 10 natural, Ri < 0.1 forced, in between h³ = h_forced³ + h_natural³).
Inside the cabin natural convection is used, as in PW4. Air properties come from CoolProp at the
film temperature. The program prints, per surface and per side, the correlation used, Re, Ra, Ri,
Nu and whether the validity range is respected.

**Sky temperature**: Martin–Berdahl clear-sky model of the course notes.

**Occupants** (`occupant_heat`): Dubois area, skin temperature, sensible heat through the clothing
and by respiration, the rest of the metabolic heat is latent (PW5).

**Cabin and ventilation** (`ventilation`, `air_handling`).

```
energy :  ṁ·cp·(T_vent − T_cab) + P_sens + Σ Q_cabin = 0
water  :  x_cab = x_supply + ṁ_w/ṁ          ṁ_w = P_lat / 2501.6 kJ/kg
power  :  Q_hx = ṁ·[h(T_vent, x_supply) − h_mix]
```

First `T_vent` is computed with the hygienic minimum of fresh air. If it falls outside the allowed
range, `T_vent` is fixed at the limit, the flow is increased and everything above the fresh-air
minimum is recirculated. If the mixed air holds more water than saturated air at `T_vent`, the
surplus condenses on the cooling coil.

## Inputs

All inputs are at the top of `panther_6x6.jl`, each marked `[SPEC]` (Rosenbauer documents),
`[NORM]` (standard) or `[ASSUMED]` (estimate to justify in the report).

| Input | Value | Source |
|---|---|---|
| Cab box | 2.6 × 3.0 × 2.0 m | width from the drawing, rest estimated from the drawing |
| Glass | 9.75 m² (windshield 5.0, sides 2 × 1.9, roof 0.95) | total from brochure (105 ft²), split estimated |
| Opaque wall | 3 mm aluminium + 25 mm foam (0.040 W/mK) + 3 mm lining | materials from spec, thicknesses estimated |
| Windshield | laminated, solar transmittance 0.35 | spec (solar control film) |
| Other glass | tinted tempered, τ = 0.50, α = 0.43 | estimated |
| Cooling case | 43.5 °C, RH 15 %, 1000 W/m² sun at 62.6°, standing still, cabin 22 °C | 43.5 °C = upper design limit of the spec |
| Heating case | 0 °C, RH 80 %, night, 80 km/h, cabin 20 °C | 0 °C = lower design limit of the spec; 2/3 of top speed and 20 °C from MIL-STD-1472F 5.12.6.1 |
| Occupants | 6 (all seats) | brochure |
| Fresh air | 0.57 m³/min per person | MIL-STD-1472F 5.12.6.2 |
| Supply air limits | 10 °C to 40 °C | 40 °C as in PW4, 10 °C estimated |

## Results with these inputs

| | Cooling case | Heating case |
|---|---|---|
| Heat through walls and glass | +3068 W | −1205 W |
| of which sunlight through the glass | 1547 W | 0 W |
| Occupants + equipment, sensible | 430 W | 665 W |
| `T_vent` with fresh air only | −32.8 °C (not acceptable) | 27.3 °C (acceptable) |
| Supply temperature used | 10 °C | 27.3 °C |
| Supply flow | 14.0 m³/min | 3.8 m³/min |
| Recirculation | 78 % | 0 % |
| Cabin relative humidity | 51 % | 30 % |
| **Heat exchanger power** | **5.4 kW cooling** | **2.0 kW heating** |

## Verification

`validation_pw4.jl` reproduces PW4 Ex 4.2 (armoured vehicle, −5 °C, 50 km/h, 8 soldiers):

| | This code | PW4 solution |
|---|---|---|
| Heat loss through the walls | 4696 W | 4665 W |
| `T_vent` with the hygienic flow | 64.2 °C | 64.1 °C |
| Flow with `T_vent` = 40 °C | 0.1932 kg/s | 0.1927 kg/s |

The 0.7 % difference on the wall loss comes from the air properties (CoolProp instead of tables)
and the exponent of Pr in the turbulent flat-plate formula.

## Limitations

- The cab is a box; the curved nose and the wheel arches are not modelled.
- Radiation between the inside surfaces is neglected (as in PW4) and the cabin air is perfectly mixed.
- Inside convection is natural convection; the ventilation jets are not modelled. Imposing
  `h_inside = 8.0` in a `Case` changes the result to 6.5 kW cooling and 3.0 kW heating.
- The front of the cab is treated as a flat plate in parallel flow, as in PW4.
- Steady state only: no warm-up or cool-down time.
- Ideal cooling coil: the air leaves it saturated at the supply temperature when water condenses.

## Still to do

- Sensitivity analysis on the three most uncertain vehicle properties.
- Duct system with at least four outlets: diameters, flows, pressure drops, fan power.
