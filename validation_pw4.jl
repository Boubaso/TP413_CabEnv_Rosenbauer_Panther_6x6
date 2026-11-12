#=
validation_pw4.jl
=================
Check of `cabin_model.jl` against an exercise with a known answer: PW4, Ex 4.2
(armoured vehicle in winter, 50 km/h, 8 soldiers, cabin at 20 °C).

Reference values of the worked solution:
    heat lost through the walls        ≈ 4665 W
    T_vent with the hygienic flow      ≈ 64 °C
    flow when T_vent is limited to 40 °C ≈ 0.193 kg/s

Small differences (a few %) are expected: the worked solution uses tabulated air properties
and a slightly different exponent of Pr in the turbulent flat-plate formula.

Run:  julia validation_pw4.jl
=#

include("cabin_model.jl")

# --- vehicle of PW4: box 6.5 x 2.6 x 1.9 m, 1 cm steel + 1.5 cm spall liner ---
lv, wv, hv = 6.5, 2.6, 1.9
armour = Wall(name = "steel + spall liner",
              layers = [Layer("steel", 0.010, 48.0), Layer("spall liner", 0.015, 0.51)],
              α_sol = 0.65, ε_ir = 0.90)

pw4_surfaces = [
    Surface(name = "Roof",  wall = armour, A = lv * wv, tilt = 0.0,   face = (lv, wv)),
    Surface(name = "Floor", wall = armour, A = lv * wv, tilt = 180.0, face = (lv, wv), sheltered = true),
    Surface(name = "Front", wall = armour, A = wv * hv, tilt = 90.0,  face = (hv, wv), L_vert = hv),
    Surface(name = "Rear",  wall = armour, A = wv * hv, tilt = 90.0,  face = (hv, wv), L_vert = hv),
    Surface(name = "Side E", wall = armour, A = lv * hv, tilt = 90.0, face = (lv, hv)),
    Surface(name = "Side W", wall = armour, A = lv * hv, tilt = 90.0, face = (lv, hv)),
]

# --- conditions of Ex 4.2: -5 °C, RH 80 %, fully clouded (diffuse light 50 W/m²), 50 km/h ---
pw4_case = Case(name = "PW4 Ex 4.2", T_amb = degC(-5.0), RH_amb = 0.80, v = 50 / 3.6,
                G_dif = 50.0,
                T_sky = 257.76,             # cloudy sky, value computed in the worked solution
                ground_radiation = false,   # neglected in the statement
                T_cab = degC(20.0), T_vent_max = degC(40.0))

amb     = ambient(pw4_case)
walls   = [solve_surface(s, pw4_case, amb) for s in pw4_surfaces]
Q_walls = sum(w.Q_cabin for w in walls)
print_surfaces(walls)

# 8 soldiers, 100 W sensible and 50 W latent each, 0.5 m³/min of fresh air per person
vent = ventilation(pw4_case, amb, Q_walls, 8 * 100.0, 8 * 50.0, 8 * 0.5 / 60)
print_ventilation(pw4_case, amb, vent)

println("\n  Comparison with the worked solution of PW4")
@printf("    %-36s %10s %10s %8s\n", "", "this code", "PW4", "diff")
for (label, value, reference) in (
        ("wall heat flow [W]",              Q_walls,                 -4665.0),
        ("roof heat flow [W]",              walls[1].Q_cabin,        -1574.8),
        ("floor heat flow [W]",             walls[2].Q_cabin,         -311.6),
        ("side heat flow [W]",              walls[5].Q_cabin,         -989.6),
        ("fresh air flow [kg/s]",           vent.ṁ_fresh,             0.08747),
        ("T_vent, hygienic flow [degC]",    in_degC(vent.T_vent_hyg), 64.1),
        ("flow at T_vent = 40 degC [kg/s]", vent.ṁ_total,             0.1927))
    @printf("    %-36s %10.4g %10.4g %7.1f%%\n", label, value, reference, 100 * (value - reference) / abs(reference))
end
