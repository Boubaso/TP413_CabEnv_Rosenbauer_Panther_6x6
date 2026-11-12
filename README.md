# TP413_CabEnv_Rosenbauer_Panther_6x6
Project for the course TP413 with Col. Jansens and Mr. Degreef. A quantative analysis on the cabin environment conditions needed in the airport firefighting vehicule Panther 6x6 constructed by Rosenbauer.

## What the code does

It computes, for one cooling case and one heating case, the steady-state energy balance of the
crew cab of the Panther 6x6:

1. energy balance of every wall and window (sun, outside convection, infrared radiation, conduction, inside convection)
2. heat and water vapour released by the occupants
3. energy and water balance of the cabin air → required temperature of the ventilation air
4. if that temperature is not acceptable: fixed supply temperature, larger flow, recirculation
5. power of the heating / cooling heat exchanger

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

## Theory 

Since an exhaustive revision of the theory is not demanded for this project, we will refer to the report itself which will refer in his turn to the course material of TP413.

