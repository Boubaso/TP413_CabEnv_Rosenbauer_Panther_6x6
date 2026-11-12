#=
panther_6x6.jl
==============
TP413 project - cabin environment of the Rosenbauer Panther 6x6 (airport fire-fighting vehicle).

This is the file to RUN. It contains the DATA of the vehicle and of the two design cases
(one cooling case, one heating case). The physics is in `cabin_model.jl`.

Run:  julia panther_6x6.jl

Every value marked  [SPEC]  comes from the Rosenbauer documents,
every value marked  [NORM]  comes from a standard,
every value marked  [ASSUMED]  is an engineering estimate: to justify (or change) in the report.
=#

include("cabin_model.jl")

# ==========================================================================================
# 1. THE CABIN: geometry
# ==========================================================================================
# The crew cab is simplified to a box (as in PW4). Dimensions estimated from the dimensional
# drawing of the brochure (vehicle width 3.0 m, height 3.65 m).

L_cab = 2.6     # length, in the driving direction [m]                         [ASSUMED from drawing]
W_cab = 3.0     # width  [m]                                                   [SPEC: vehicle width]
H_cab = 2.0     # height floor -> roof [m]                                     [ASSUMED from drawing]

# Glazing. The brochure gives about 105 ft² = 9.75 m² of glass in total; split estimated:
A_windshield  = 5.0     # [m²]                                                 [ASSUMED]
A_side_window = 1.9     # per side [m²]                                        [ASSUMED]
A_roof_window = 0.95    # [m²]                                                 [ASSUMED]
tilt_windshield = 75.0  # the windshield leans backwards: 90° = vertical       [ASSUMED]

heading = 180.0         # compass direction the nose points to (180° = south = towards the sun at noon)

# ==========================================================================================
# 2. THE CABIN: wall constructions
# ==========================================================================================
# [SPEC] welded aluminium box sections with aluminium sheets, insulated; laminated windshield,
#        tempered side and roof windows; solar-control film: the windshield lets ≈ 35 % of the
#        solar radiation through; the other windows are tinted.
# Thicknesses and conductivities are typical values                            [ASSUMED]

paint_α = 0.65      # solar absorptance of the paint (red/yellow)              [ASSUMED]
paint_ε = 0.90      # infrared emissivity of paint                             [ASSUMED]
glass_ε = 0.84      # infrared emissivity of uncoated glass

cab_wall = Wall(name = "aluminium + insulation + lining",
                layers = [Layer("aluminium sheet",   0.003, 160.0),
                          Layer("insulation foam",   0.025, 0.040),
                          Layer("interior lining",   0.003, 0.20)],
                α_sol = paint_α, ε_ir = paint_ε)

cab_floor = Wall(name = "aluminium + insulation + rubber mat",
                 layers = [Layer("aluminium sheet",  0.004, 160.0),
                           Layer("insulation foam",  0.025, 0.040),
                           Layer("rubber floor mat", 0.005, 0.16)],
                 α_sol = paint_α, ε_ir = paint_ε)

# Glass: sunlight is partly transmitted (τ), partly absorbed (α), the rest is reflected.
windshield_glass = Wall(name = "laminated glass + solar film",
                        layers = [Layer("glass", 0.004, 1.0), Layer("PVB foil", 0.00076, 0.20),
                                  Layer("glass", 0.004, 1.0)],
                        τ_sol = 0.35,          # [SPEC]
                        α_sol = 0.40,          # [ASSUMED] (so 25 % is reflected)
                        ε_ir = glass_ε)

tinted_glass = Wall(name = "tinted tempered glass",
                    layers = [Layer("glass", 0.005, 1.0)],
                    τ_sol = 0.50, α_sol = 0.43,   # [ASSUMED] typical grey-tinted glass
                    ε_ir = glass_ε)

# ==========================================================================================
# 3. THE CABIN: list of heat-transfer surfaces                             (project item 5.a)
# ==========================================================================================
# face = (length along the driving direction, other length) of the WHOLE face of the box.
# The floor and the rear wall are "sheltered": they get no sun and do not see the sky.
#   - floor: sees the shaded ground/chassis under the cab, the driving wind passes under it
#   - rear wall: faces the superstructure (pump/tank body) -> no wind either

roof_face  = (L_cab, W_cab)
side_face  = (L_cab, H_cab)
front_face = (H_cab, W_cab)     # the air flows upwards over the front: flow length = height

surfaces = [
    Surface(name = "Roof",          wall = cab_wall,     A = L_cab * W_cab - A_roof_window,
            tilt = 0.0, face = roof_face),
    Surface(name = "Roof window",   wall = tinted_glass, A = A_roof_window,
            tilt = 0.0, face = roof_face),
    Surface(name = "Floor",         wall = cab_floor,    A = L_cab * W_cab,
            tilt = 180.0, face = roof_face, sheltered = true),
    Surface(name = "Front panel",   wall = cab_wall,     A = W_cab * H_cab - A_windshield,
            tilt = 90.0, azimuth = heading, face = front_face, L_vert = H_cab),
    Surface(name = "Windshield",    wall = windshield_glass, A = A_windshield,
            tilt = tilt_windshield, azimuth = heading, face = front_face, L_vert = H_cab),
    Surface(name = "Left side",     wall = cab_wall,     A = L_cab * H_cab - A_side_window,
            tilt = 90.0, azimuth = heading - 90, face = side_face),
    Surface(name = "Left window",   wall = tinted_glass, A = A_side_window,
            tilt = 90.0, azimuth = heading - 90, face = side_face),
    Surface(name = "Right side",    wall = cab_wall,     A = L_cab * H_cab - A_side_window,
            tilt = 90.0, azimuth = heading + 90, face = side_face),
    Surface(name = "Right window",  wall = tinted_glass, A = A_side_window,
            tilt = 90.0, azimuth = heading + 90, face = side_face),
    Surface(name = "Rear wall",     wall = cab_wall,     A = W_cab * H_cab,
            tilt = 90.0, azimuth = heading + 180, face = front_face, L_vert = H_cab,
            sheltered = true, wind = 0.0),
]

# ==========================================================================================
# 4. OCCUPANTS
# ==========================================================================================
# [SPEC] seating for up to 6 persons -> design with a full cab.
# Metabolic rates from the table of the course notes:
#   seated, quiet 60 W/m² ; driving a heavy vehicle 185 W/m².
# Clothing: insulated coveralls 1.37 clo, f_cl = 1.26 (defaults of `Person`).

seated_crew_member = Person(M = 60.0)
driver_driving     = Person(M = 185.0)

# ==========================================================================================
# 5. DESIGN CASES                                                           (project item 2-4)
# ==========================================================================================
# Fresh air: 0.57 m³/min per person   [NORM: MIL-STD-1472F 5.12.6.2]  (default of `Case`)
# Supply air: at most 40 °C (as in PW4), at least 10 °C                         [ASSUMED]

v_max = 120 / 3.6      # top speed [m/s]                                           [SPEC]

# --- COOLING: hottest day of the design range, vehicle on stand-by on the apron in full sun ---
cooling = Case(
    name   = "COOLING CASE - 43.5 degC, full sun, standing still, 6 occupants",
    T_amb  = degC(43.5),            # upper limit of the design temperature range  [SPEC]
    RH_amb = 0.15,                  # hot and dry                                  [ASSUMED]
    v      = 0.0,                   # standing still, no wind: worst case
    hour   = 12.0,
    DNI    = 1000.0,                # clear day (course notes)
    sun_elevation = 62.6,           # solar noon on 21 June at 50.8° N: 90 - 50.8 + 23.4
    sun_azimuth   = 180.0,          # the sun is in the south at solar noon
    T_ground = degC(43.5 + 15.0),   # sunlit concrete, warmer than the air         [ASSUMED]
    T_cab  = degC(22.0),            # upper end of the comfort range 20-22 °C
    crew   = fill(seated_crew_member, 6),
    P_equipment = 300.0,            # radios, displays, electronics                [ASSUMED]
)

# --- HEATING: coldest day of the design range, at night, driving ---
# [NORM: MIL-STD-1472F 5.12.6.1] heater sized at the minimum design temperature with the
# vehicle moving at 2/3 of its maximum speed; cabin above 20 °C.
heating = Case(
    name   = "HEATING CASE - 0 degC, night, driving at 80 km/h, 6 occupants",
    T_amb  = degC(0.0),             # lower limit of the design temperature range  [SPEC]
    RH_amb = 0.80,                  #                                              [ASSUMED]
    v      = 2 / 3 * v_max,         # 80 km/h                                      [NORM]
    hour   = 4.0,                   # night: no sun, clear (cold) sky
    T_cab  = degC(20.0),            #                                              [NORM]
    crew   = [driver_driving; fill(seated_crew_member, 5)],
    P_equipment = 0.0,              # no help from equipment: worst case
)

# ==========================================================================================
# 6. SOLVE AND PRINT
# ==========================================================================================

println("Rosenbauer Panther 6x6 - crew cab model")
@printf("  box %.1f x %.1f x %.1f m, volume %.1f m3, glass %.2f m2, opaque %.2f m2\n",
        L_cab, W_cab, H_cab, L_cab * W_cab * H_cab,
        sum(s.A for s in surfaces if s.wall.τ_sol > 0), sum(s.A for s in surfaces if s.wall.τ_sol == 0))
for w in (cab_wall, cab_floor, windshield_glass, tinted_glass)
    @printf("  %-38s R_cond = %.4f m2K/W\n", w.name, R_cond(w))
end

results = [solve_case(surfaces, case) for case in (cooling, heating)]
foreach(print_report, results)

# --- Check against the insulation requirement of the vehicle specification ---
# [SPEC] average heat loss of the cab ≤ 0.24 BTU/(h ft² °F) = 1.36 W/m²K (opaque walls).
println("\n", "="^120)
println(" Checks")
println("="^120)
U_spec = 0.24 * 5.678
for r in results
    opaque = [w for (w, s) in zip(r.walls, surfaces) if s.wall.τ_sol == 0]
    U_mean = sum(w.U * w.A for w in opaque) / sum(w.A for w in opaque)
    @printf("  %-8s mean U of the opaque walls = %.2f W/m2K  (specification: max %.2f)\n",
            split(r.case.name)[1], U_mean, U_spec)
end

# --- Summary ---
V_cab = L_cab * W_cab * H_cab
println("\n", "="^120)
println(" Summary")
println("="^120)
@printf("  %-34s %14s %14s\n", "", "cooling case", "heating case")
row(label, f) = @printf("  %-34s %14.1f %14.1f\n", label, f(results[1]), f(results[2]))
row("heat through walls + glass [W]",  r -> r.Q_walls)
row("  of which sun through glass [W]", r -> sum(w.Q_trans for w in r.walls))
row("occupants + equipment, sens. [W]", r -> r.P_sens)
row("occupants, latent [W]",            r -> r.P_lat)
row("T_vent with fresh air only [degC]", r -> in_degC(r.vent.T_vent_hyg))
row("supply temperature used [degC]",   r -> in_degC(r.vent.T_vent))
row("supply flow [m3/min]",             r -> 60 * r.vent.V_supply)
row("fresh air flow [m3/min]",          r -> 60 * r.vent.V_fresh)
row("recirculation [%]",                r -> 100 * r.vent.recirculation)
row("air changes per hour [1/h]",       r -> 3600 * r.vent.V_supply / V_cab)
row("cabin relative humidity [%]",      r -> 100 * r.vent.RH_cab)
row("HEAT EXCHANGER POWER [W]",         r -> r.vent.Q_hx)
