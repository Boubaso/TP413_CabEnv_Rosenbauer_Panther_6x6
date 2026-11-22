#=
panther_6x6.jl
==============
TP413 project - cabin environment of the Rosenbauer Panther 6x6 (airport fire-fighting vehicle).

This is the file to RUN. It contains the DATA of the vehicle and of the two design cases
(one cooling case, one heating case). The physics is in `cabin_model.jl`.

Run:  julia panther_6x6.jl
=#

include("cabin_model.jl")

# ==========================================================================================
# 1. THE CABIN: geometry
# ==========================================================================================


L_cab = 2.6     # length at floor level, in the driving direction [m]          
W_cab = 3.0     # width  [m]                                                   
H_cab = 2.0     # height floor -> roof [m]                                     

A_windshield  = 5.0     # [m²]                                                 
A_side_window = 1.9     # per side [m²]                                    
A_roof_window = 0.95    # [m²]                                                 

rake_windshield = 40.0                     # [°]                               
tilt_windshield = 90.0 - rake_windshield   # angle between its normal and the vertical (model convention)

# Everything else
slant_windshield = A_windshield / W_cab                 # length of the glass along the slope [m]
rise_windshield  = slant_windshield * cosd(rake_windshield)   # height it covers [m]
run_windshield   = slant_windshield * sind(rake_windshield)   # how far it leans back [m]
h_sill = H_cab - rise_windshield       # height of the vertical front panel under the glass [m]
L_roof = L_cab - run_windshield        # the roof is shorter than the floor [m]
@assert h_sill > 0 && L_roof > 0 "the windshield does not fit in the cab: check the geometry"

# A side wall is the rectangle L_cab x H_cab minus the triangle cut off by the windshield:
A_side_face = L_cab * H_cab - run_windshield * rise_windshield / 2    # [m²]
L_side      = A_side_face / H_cab      # length of the rectangle with the same height and area [m]
V_cab       = A_side_face * W_cab      # volume of the prism = cross-section x width [m³]

heading = 180.0         # compass direction the nose points to (180° = south = towards the sun at noon)

# ==========================================================================================
# 2. THE CABIN: wall constructions
# ==========================================================================================

paint_α = 0.65      # solar absorptance of the paint (red/yellow)             
paint_ε = 0.90      # infrared emissivity of paint                             
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
                        τ_sol = 0.35,          
                        α_sol = 0.40,          # (so 25 % is reflected)
                        ε_ir = glass_ε)

tinted_glass = Wall(name = "tinted tempered glass",
                    layers = [Layer("glass", 0.005, 1.0)],
                    τ_sol = 0.50, α_sol = 0.43,   # typical grey-tinted glass
                    ε_ir = glass_ε)

# ==========================================================================================
# 3. THE CABIN: list of heat-transfer surfaces                             
# ==========================================================================================

roof_face       = (L_roof, W_cab)
floor_face      = (L_cab, W_cab)
side_face       = (L_side, H_cab)               # rectangle with the same height and area
panel_face      = (h_sill, W_cab)               # the air flows upwards over the front
windshield_face = (slant_windshield, W_cab)     # ... and then up along the glass
rear_face       = (H_cab, W_cab)

surfaces = [
    Surface(name = "Roof",          wall = cab_wall,     A = L_roof * W_cab - A_roof_window,
            tilt = 0.0, face = roof_face),
    Surface(name = "Roof window",   wall = tinted_glass, A = A_roof_window,
            tilt = 0.0, face = roof_face),
    Surface(name = "Floor",         wall = cab_floor,    A = L_cab * W_cab,
            tilt = 180.0, face = floor_face, sheltered = true),
    Surface(name = "Front panel",   wall = cab_wall,     A = W_cab * h_sill,
            tilt = 90.0, azimuth = heading, face = panel_face, L_slope = h_sill),
    Surface(name = "Windshield",    wall = windshield_glass, A = A_windshield,
            tilt = tilt_windshield, azimuth = heading, face = windshield_face,
            L_slope = slant_windshield),
    Surface(name = "Left side",     wall = cab_wall,     A = A_side_face - A_side_window,
            tilt = 90.0, azimuth = heading - 90, face = side_face),
    Surface(name = "Left window",   wall = tinted_glass, A = A_side_window,
            tilt = 90.0, azimuth = heading - 90, face = side_face),
    Surface(name = "Right side",    wall = cab_wall,     A = A_side_face - A_side_window,
            tilt = 90.0, azimuth = heading + 90, face = side_face),
    Surface(name = "Right window",  wall = tinted_glass, A = A_side_window,
            tilt = 90.0, azimuth = heading + 90, face = side_face),
    Surface(name = "Rear wall",     wall = cab_wall,     A = W_cab * H_cab,
            tilt = 90.0, azimuth = heading + 180, face = rear_face, L_slope = H_cab,
            sheltered = true, wind = 0.0),
]

# ==========================================================================================
# 4. OCCUPANTS
# ==========================================================================================
#  seating for up to 6 persons -> design with a full cab.
# Metabolic rates from the table of the course notes:

seated_crew_member = Person(M = 60.0)
driver_driving     = Person(M = 185.0)

# ==========================================================================================
# 5. DESIGN CASES                                                          
# ==========================================================================================

v_max = 120 / 3.6      # top speed [km/h] --> [m/s]                                          

# --- COOLING: hottest day of the design range, vehicle on stand-by on the apron in full sun ---
cooling = Case(
    name   = "COOLING CASE - 43.5 degC, full sun, standing still, 6 occupants",
    T_amb  = degC(43.5),            # upper limit of the design temperature range  
    RH_amb = 0.15,                  # hot and dry                                  
    v      = 0.0,                   # standing still, no wind: worst case
    hour   = 12.0,
    DNI    = 1000.0,                # clear day (course notes)
    sun_elevation = 62.6,           # solar noon on 21 June at 50.8° N: 90 - 50.8 + 23.4
    sun_azimuth   = 180.0,          # the sun is in the south at solar noon
    T_ground = degC(43.5 + 15.0),   # sunlit concrete, warmer than the air       
    T_cab  = degC(22.0),            # upper end of the comfort range 20-22 °C
    crew   = fill(seated_crew_member, 6),
    P_equipment = 300.0,            # radios, displays, electronics               
)

# --- HEATING: coldest day of the design range, at night, driving 2/3 speed max ---
heating = Case(
    name   = "HEATING CASE - 0 degC, night, driving at 80 km/h, 6 occupants",
    T_amb  = degC(0.0),             # lower limit of the design temperature range  
    RH_amb = 0.80,                  #                                              
    v      = 2 / 3 * v_max,         # 80 km/h                                     
    hour   = 4.0,                   # night: no sun, clear (cold) sky
    T_cab  = degC(20.0),            #                                             
    crew   = [driver_driving; fill(seated_crew_member, 5)],
    P_equipment = 0.0,              # no help from equipment: worst case
)

# ==========================================================================================
# 6. SOLVE AND PRINT
# ==========================================================================================

println("Rosenbauer Panther 6x6 - crew cab model")
@printf("  cab %.1f x %.1f x %.1f m with an inclined windshield, volume %.1f m3, glass %.2f m2, opaque %.2f m2\n",
        L_cab, W_cab, H_cab, V_cab,
        sum(s.A for s in surfaces if s.wall.τ_sol > 0), sum(s.A for s in surfaces if s.wall.τ_sol == 0))
@printf("  windshield: %.0f deg from the vertical, slant %.2f m (rise %.2f m, run %.2f m)\n",
        rake_windshield, slant_windshield, rise_windshield, run_windshield)
@printf("  front panel height %.2f m, roof length %.2f m, side face %.2f m2 each\n",
        h_sill, L_roof, A_side_face)
for w in (cab_wall, cab_floor, windshield_glass, tinted_glass)
    @printf("  %-38s R_cond = %.4f m2K/W\n", w.name, R_cond(w))
end

results = [solve_case(surfaces, case) for case in (cooling, heating)]
foreach(print_report, results)

# --- Check against the insulation requirement of the vehicle specification ---
# average heat loss of the cab ≤ 0.24 BTU/(h ft² °F) = 1.36 W/m²K (opaque walls).
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
