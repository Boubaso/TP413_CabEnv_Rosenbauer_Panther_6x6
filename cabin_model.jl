#=
cabin_model.jl
==============
Steady-state thermal model of a vehicle cabin.

This file contains ONLY the physics (functions and data types). It knows nothing
about the Panther: the vehicle data and the design cases are in `panther_6x6.jl`.

Structure of the file (same order as the calculation):
    1. Constants and small helpers
    2. Moist air (psychrometrics)
    3. Data types: Layer, Wall, Surface, Person, Case
    4. Ambient conditions (humidity, sky temperature)
    5. Convection coefficients (natural / forced / mixed)
    6. Energy balance of ONE surface                (project item 5.a - 5.d)
    7. Heat released by the occupants               (PW5)
    8. Cabin balance, ventilation, recirculation    (project item 5.e - 5.h and 6)
    9. Solve a complete case
   10. Printing of the results

UNITS: everything inside the model is in plain SI:
    temperature [K], length [m], area [m²], power [W], pressure [Pa],
    mass flow [kg dry air / s], humidity ratio x [kg water / kg dry air].
Temperatures are typed in °C in the input file with `degC(...)` and printed in °C.

SIGN CONVENTION (same as PW4): a heat flow is POSITIVE when it flows TOWARDS THE CABIN.
=#

using CoolProp      # air and moist-air properties
using Printf        # formatted printing

# ------------------------------------------------------------------------------------------
# 1. Constants and small helpers
# ------------------------------------------------------------------------------------------

const σ_SB  = 5.670374e-8   # Stefan-Boltzmann constant [W/m²K⁴]
const g0    = 9.81          # gravity [m/s²]
const T_0C  = 273.15        # 0 °C in kelvin

# Moist-air constants of the course notes: h = 1003·t + x·(2501.6e3 + 1863·t)
const cp_da = 1003.0        # specific heat of dry air [J/kgK]
const cp_wv = 1863.0        # specific heat of water vapour [J/kgK]
const h_fg0 = 2501.6e3      # latent heat of evaporation of water at 0 °C [J/kg]

degC(t)     = t + T_0C      # °C -> K   (use this to type inputs)
in_degC(T)  = T - T_0C      # K -> °C   (use this to print)

"""
    bisection(f, a, b)

Find the root of `f` between `a` and `b` by halving the interval again and again.
`f(a)` and `f(b)` must have opposite signs. Slow but it can not fail, which is what we want.
"""
function bisection(f, a, b; iterations = 100)
    fa = f(a)
    fa * f(b) <= 0 || error("bisection: the root is not between $a and $b")
    for _ in 1:iterations
        m  = (a + b) / 2
        fm = f(m)
        if fa * fm <= 0     # the root is in the left half
            b = m
        else                # the root is in the right half
            a, fa = m, fm
        end
    end
    return (a + b) / 2
end

# ------------------------------------------------------------------------------------------
# 2. Moist air (psychrometrics)
# ------------------------------------------------------------------------------------------
# x = humidity ratio [kg water vapour / kg dry air], RH = relative humidity [0..1].
# CoolProp's HAPropsSI is used for the saturation data (it also works below 0 °C).

humidity_ratio(T, RH, p) = HAPropsSI("W", "T", T, "P", p, "R", RH)       # x from T and RH
dew_point(T, RH, p)      = HAPropsSI("D", "T", T, "P", p, "R", RH)       # dew point [K]
x_saturation(T, p)       = HAPropsSI("W", "T", T, "P", p, "R", 1.0)      # x of saturated air
volume_per_kg_dry_air(T, x, p) = HAPropsSI("Vda", "T", T, "P", p, "W", x) # [m³/kg dry air]

"Relative humidity [0..1] of air at temperature T with humidity ratio x (1.0 if saturated)."
function relative_humidity(T, x, p)
    x >= x_saturation(T, p) && return 1.0
    return HAPropsSI("R", "T", T, "P", p, "W", x)
end

"Enthalpy of moist air per kg DRY air [J/kg] (course notes), T in kelvin."
function enthalpy(T, x)
    t = in_degC(T)
    return cp_da * t + x * (h_fg0 + cp_wv * t)
end

"Inverse of `enthalpy`: temperature [K] of moist air with enthalpy h and humidity ratio x."
temperature_from_enthalpy(h, x) = degC((h - x * h_fg0) / (cp_da + x * cp_wv))

# ------------------------------------------------------------------------------------------
# 3. Data types
# ------------------------------------------------------------------------------------------

"One material layer of a wall: thickness d [m] and thermal conductivity k [W/mK]."
struct Layer
    name::String
    d::Float64
    k::Float64
end

"""
A wall construction = layers in series + the optical properties of its outside face.
For glass, `τ_sol` > 0: that fraction of the sunlight goes straight through into the cabin.
"""
Base.@kwdef struct Wall
    name::String
    layers::Vector{Layer}
    α_sol::Float64            # solar absorptance (fraction of sunlight absorbed)        [-]
    ε_ir::Float64 = 0.90      # thermal (infrared) emissivity of the outside face        [-]
    τ_sol::Float64 = 0.0      # solar transmittance at perpendicular incidence (glass)   [-]
end

"Conduction resistance of a wall per m² [m²K/W]: layers in series, R'' = Σ d/k."
R_cond(w::Wall) = sum(l.d / l.k for l in w.layers)

"""
One heat-transfer surface of the cabin (a wall panel or a window).

`tilt` is the angle between the outward normal and the vertical:
    0° = faces the sky (roof), 90° = vertical wall, 180° = faces the ground (floor).
`azimuth` is the compass direction the surface faces (0° = N, 90° = E, 180° = S, 270° = W).

A window and the panel around it are two `Surface`s that belong to the same FACE of the
cabin. The convection correlations need the size of the whole face, not of the part:
    `face`   = (length along the driving direction, other length) of the whole face [m]
    `L_flow` = length the outside air travels along the face (forced convection)    [m]
    `L_vert` = vertical height of the face (natural convection on a wall)           [m]
"""
Base.@kwdef struct Surface
    name::String
    wall::Wall
    A::Float64                       # area of this surface [m²]
    tilt::Float64                    # [°]
    azimuth::Float64 = 0.0           # [°]
    face::Tuple{Float64,Float64}     # [m]
    L_flow::Float64 = face[1]        # [m]
    L_vert::Float64 = face[2]        # [m]
    sheltered::Bool = false          # true: no sun, no sky; it "sees" objects at air temperature
    wind::Float64 = 1.0              # fraction of the vehicle speed felt by this surface (0..1)
end

"One occupant (PW5). Defaults: fire-fighter in insulated coveralls, seated."
Base.@kwdef struct Person
    M::Float64                 # metabolic rate per m² of body surface [W/m²]
    I_cl::Float64 = 1.37       # clothing insulation, standing value [clo]
    f_cl::Float64 = 1.26       # clothing area factor [-]
    seated::Bool = true        # seated: +0.15 clo (the seat insulates)
    mass::Float64 = 75.0       # [kg]
    height::Float64 = 1.75     # [m]
    h::Float64 = 7.8           # convection + radiation coefficient around the body [W/m²K]
end

"""
One design case = the ambient conditions + what we want in the cabin.
Everything that describes the "operating conditions" of the report is in here.
"""
Base.@kwdef struct Case
    name::String
    # --- ambient -------------------------------------------------------------------------
    T_amb::Float64                    # outside air temperature [K]
    RH_amb::Float64                   # outside relative humidity [0..1]
    p::Float64 = 101325.0             # atmospheric pressure [Pa]
    v::Float64 = 0.0                  # vehicle speed = air speed over the cab [m/s]
    hour::Float64 = 12.0              # solar time [h], only used for the sky temperature
    DNI::Float64 = 0.0                # direct solar irradiance, perpendicular to the rays [W/m²]
    G_dif::Float64 = 0.0              # diffuse irradiance on a horizontal surface [W/m²]
    sun_elevation::Float64 = 0.0      # height of the sun above the horizon [°]
    sun_azimuth::Float64 = 180.0      # compass direction of the sun [°] (180° = south)
    T_ground::Float64 = T_amb         # temperature of the ground around the vehicle [K]
    T_sky::Union{Nothing,Float64} = nothing   # nothing = clear sky (Martin-Berdahl); or impose [K]
    ground_radiation::Bool = true     # false = neglect radiation exchange with the ground
    # --- cabin ---------------------------------------------------------------------------
    T_cab::Float64                    # required cabin air temperature [K]
    T_vent_min::Float64 = degC(10.0)  # coldest allowed supply air [K]
    T_vent_max::Float64 = degC(40.0)  # hottest allowed supply air [K]
    V_fresh_per_person::Float64 = 0.57 / 60   # fresh outside air per person [m³/s]
    crew::Vector{Person} = Person[]
    P_equipment::Float64 = 0.0        # heat released by equipment in the cabin [W]
    h_inside::Union{Nothing,Float64} = nothing  # nothing = natural-convection correlations;
                                                # or impose one inside coefficient [W/m²K]
end

# ------------------------------------------------------------------------------------------
# 4. Ambient conditions
# ------------------------------------------------------------------------------------------

"""
Clear-sky temperature [K], Martin-Berdahl model (course notes / PW3):

    ε_sky = 0.711 + 0.0056·t_dp + 7.3e-5·t_dp² + 0.013·cos(2π·hour/24)     (t_dp in °C)
    T_sky = T_amb · ε_sky^(1/4)

Humid air (high dew point) -> emissive sky -> warm sky. Dry air -> cold sky.
"""
function sky_temperature(T_amb, T_dp, hour)
    t_dp  = in_degC(T_dp)
    ε_sky = 0.711 + 0.0056 * t_dp + 7.3e-5 * t_dp^2 + 0.013 * cos(2π * hour / 24)
    return T_amb * ε_sky^(1 / 4)
end

"Everything about the outside air that is derived from the case inputs."
function ambient(c::Case)
    x    = humidity_ratio(c.T_amb, c.RH_amb, c.p)
    T_dp = dew_point(c.T_amb, c.RH_amb, c.p)
    T_sky = c.T_sky === nothing ? sky_temperature(c.T_amb, T_dp, c.hour) : c.T_sky
    v_da = volume_per_kg_dry_air(c.T_amb, x, c.p)
    return (x = x, T_dp = T_dp, T_sky = T_sky, v_da = v_da, h = enthalpy(c.T_amb, x))
end

# ------------------------------------------------------------------------------------------
# 5. Convection coefficients
# ------------------------------------------------------------------------------------------
# Every function returns a NamedTuple with the coefficient h AND the numbers needed to
# justify it in the report (which correlation, Re, Ra, Ri, Nu, inside validity range?).

const Re_transition    = 5e5    # laminar -> turbulent boundary layer on a flat plate
const HORIZONTAL_LIMIT = 5.0    # a surface within 5° of horizontal is treated as horizontal

"Properties of dry air at temperature T [K] and pressure p [Pa] (CoolProp)."
function air_properties(T, p)
    ρ  = PropsSI("D", "T", T, "P", p, "Air")              # density [kg/m³]
    μ  = PropsSI("viscosity", "T", T, "P", p, "Air")      # dynamic viscosity [Pa s]
    k  = PropsSI("conductivity", "T", T, "P", p, "Air")   # thermal conductivity [W/mK]
    cp = PropsSI("Cp0mass", "T", T, "P", p, "Air")        # specific heat [J/kgK]
    ν  = μ / ρ                                            # kinematic viscosity [m²/s]
    return (ρ = ρ, μ = μ, k = k, cp = cp, ν = ν, Pr = μ * cp / k, β = 1 / T)
end

"""
    natural_convection(s, tilt, T_s, T_air, p)

Natural convection on one side of surface `s`. `tilt` is the tilt of THAT side:
the outside of the roof faces up (0°), but the inside of the roof faces down (180°).

Three configurations (A Heat Transfer Textbook, chapter 8):
  * vertical or inclined plate : Churchill-Chu, eq. 8.13b, L = height,
                                 gravity replaced by its component along the plate
  * horizontal, UNSTABLE       : hot side up or cold side down -> air can rise/fall freely
                                 eq. 8.35 (Ra ≥ 1e7) or 8.37a (Ra < 1e7), L = area/perimeter
  * horizontal, STABLE         : hot side down or cold side up -> air is trapped, low h
                                 eq. 8.38, L = short side of the plate
Air properties at the film temperature (average of surface and air), β = 1/T_film.
"""
function natural_convection(s::Surface, tilt, T_s, T_air, p)
    air = air_properties((T_s + T_air) / 2, p)
    ΔT  = max(abs(T_s - T_air), 0.01)                     # avoid h = 0 when ΔT = 0
    Rayleigh(L, g) = g * air.β * ΔT * L^3 / air.ν^2 * air.Pr     # Ra = Gr·Pr
    a, b = s.face

    if HORIZONTAL_LIMIT <= tilt <= 180 - HORIZONTAL_LIMIT
        L    = s.L_vert
        Ra   = Rayleigh(L, g0 * sind(tilt))               # sind(90°) = 1 for a vertical wall
        Nu   = (0.825 + 0.387 * Ra^(1 / 6) / (1 + (0.492 / air.Pr)^(9 / 16))^(8 / 27))^2
        corr = "vertical plate, Churchill-Chu (AHTT 8.13b)"
        valid = 0.1 <= Ra <= 1e12
    else
        faces_up = tilt < 90
        hot      = T_s > T_air
        if faces_up == hot                                # hot & up, or cold & down
            L  = a * b / (2 * (a + b))                    # area / perimeter
            Ra = Rayleigh(L, g0)
            if Ra >= 1e7
                Nu    = 0.14 * Ra^(1 / 3) * (1 + 0.0107 * air.Pr) / (1 + 0.01 * air.Pr)
                corr  = "horizontal unstable (AHTT 8.35)"
                valid = Ra <= 2e11
            else
                Nu    = 0.56 * Ra^(1 / 4) / (1 + (0.492 / air.Pr)^(9 / 16))^(4 / 9)
                corr  = "horizontal unstable (AHTT 8.37a)"
                valid = Ra > 1
            end
        else                                              # hot & down, or cold & up
            L     = min(a, b)
            Ra    = Rayleigh(L, g0)
            Nu    = 0.58 * Ra^(1 / 5)
            corr  = "horizontal stable (AHTT 8.38)"
            valid = 1e6 <= Ra <= 1e11
        end
    end
    return (h = Nu * air.k / L, mode = "natural", corr = corr, L = L,
            Re = NaN, Ra = Ra, Ri = Inf, Nu = Nu, valid = valid)
end

"""
    forced_convection(L, v, T_s, T_air, p)

Forced convection, air flowing at speed `v` along a flat plate of length `L` (PW4 sheet):
    laminar  (Re < 5e5) :  St = 0.664·Re^(-1/2)·Pr^(-2/3)
    mixed    (Re ≥ 5e5) :  St = (0.037·Re^(-1/5) - 871/Re)·Pr^(-2/3)
with the Stanton number St = h/(ρ·cp·v) = Nu/(Re·Pr), so Nu = St·Re·Pr.
"Mixed" = the front part of the plate is laminar, the rest turbulent.
"""
function forced_convection(L, v, T_s, T_air, p)
    air = air_properties((T_s + T_air) / 2, p)
    Re  = v * L / air.ν
    if Re < Re_transition
        St   = 0.664 * Re^(-1 / 2) * air.Pr^(-2 / 3)
        corr = "flat plate, laminar"
    else
        St   = (0.037 * Re^(-1 / 5) - 871 / Re) * air.Pr^(-2 / 3)
        corr = "flat plate, laminar+turbulent"
    end
    Nu    = St * Re * air.Pr
    valid = Re <= 1e8 && 0.6 <= air.Pr <= 60
    return (h = Nu * air.k / L, mode = "forced", corr = corr, L = L,
            Re = Re, Ra = NaN, Ri = 0.0, Nu = Nu, valid = valid)
end

"""
    outside_convection(s, T_e, c)

Convection between the outside air and the outside face (temperature T_e) of surface `s`.
Natural or forced? The Richardson number compares buoyancy with the wind:

    Ri = Gr/Re² = g·β·ΔT·L / v²
    Ri > 10   -> natural convection only
    Ri < 0.1  -> forced convection only
    between   -> mixed:  h³ = h_forced³ + h_natural³
"""
function outside_convection(s::Surface, T_e, c::Case)
    natural = natural_convection(s, s.tilt, T_e, c.T_amb, c.p)
    v = c.v * s.wind
    v < 0.1 && return natural                             # standing still

    forced = forced_convection(s.L_flow, v, T_e, c.T_amb, c.p)
    T_film = (T_e + c.T_amb) / 2
    Ri = g0 * (1 / T_film) * abs(T_e - c.T_amb) * s.L_flow / v^2
    if Ri < 0.1
        return merge(forced, (Ri = Ri,))
    elseif Ri > 10
        return merge(natural, (Ri = Ri, Re = forced.Re))
    else
        h = cbrt(forced.h^3 + natural.h^3)
        return (h = h, mode = "mixed", corr = forced.corr * " + natural", L = s.L_flow,
                Re = forced.Re, Ra = natural.Ra, Ri = Ri, Nu = h * s.L_flow / air_properties(T_film, c.p).k,
                valid = forced.valid && natural.valid)
    end
end

"""
    inside_convection(s, T_si, c)

Convection between the inside face (temperature T_si) and the cabin air.
Default = natural convection (as in PW4; the ventilation jets are not modelled).
The inside of a surface faces the opposite way of its outside: tilt_inside = 180° - tilt.
"""
function inside_convection(s::Surface, T_si, c::Case)
    if c.h_inside !== nothing
        return (h = c.h_inside, mode = "imposed", corr = "imposed value", L = NaN,
                Re = NaN, Ra = NaN, Ri = NaN, Nu = NaN, valid = true)
    end
    return natural_convection(s, 180 - s.tilt, T_si, c.T_cab, c.p)
end

# ------------------------------------------------------------------------------------------
# 6. Energy balance of ONE surface
# ------------------------------------------------------------------------------------------

# Angular dependence of the solar transmission of glass: SHGC table of the course notes,
# divided by its value at 0°. (The 90° point is added: grazing light is fully reflected.)
const SHGC_ANGLE = [0.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0]     # incidence angle θ [°]
const SHGC_VALUE = [0.86, 0.84, 0.82, 0.78, 0.67, 0.42, 0.0]

"Transmission at incidence angle θ relative to perpendicular incidence (linear interpolation)."
function angular_factor(θ)
    θ >= 90 && return 0.0
    i = findlast(a -> a <= θ, SHGC_ANGLE)
    f = (θ - SHGC_ANGLE[i]) / (SHGC_ANGLE[i+1] - SHGC_ANGLE[i])
    return (SHGC_VALUE[i] + f * (SHGC_VALUE[i+1] - SHGC_VALUE[i])) / SHGC_VALUE[1]
end

"""
    solve_surface(s, c, amb)

Two unknowns: T_e (outside face) and T_si (inside face). Two equations (PW4):

    outside node :  Q_sun + Q_conv_o + Q_rad = Q_k
    inside node  :  Q_k = Q_i

    Q_sun    = α·A·G                               absorbed sunlight
    Q_conv_o = h_o·A·(T_amb - T_e)                 outside convection
    Q_rad    = ε·A·σ·[F_sky·(T_sky⁴ - T_e⁴) + F_sur·(T_sur⁴ - T_e⁴)]   infrared radiation
    Q_k      = (T_e - T_si)/R_k                    conduction through the layers
    Q_i      = h_i·A·(T_si - T_cab)                inside convection

For glass there is one more flow, straight into the cabin: Q_trans = τ(θ)·A·G.

HOW IT IS SOLVED. The difficulty: h_o and h_i depend on the unknown temperatures.
    1. guess T_e and T_si
    2. compute h_o and h_i with these temperatures               (h is now "frozen")
    3. with frozen h the inside node is linear, so T_si follows from T_e, and only ONE
       equation in T_e is left (it contains T_e⁴) -> solved by bisection
    4. move the guess half-way towards the new temperatures and go back to 2,
       until the temperatures do not change any more.
"""
function solve_surface(s::Surface, c::Case, amb)
    w, A = s.wall, s.A
    R_k  = R_cond(w) / A                                  # conduction resistance [K/W]

    # --- what does the outside face "see"? (view factors) and how much sun does it get? ---
    if s.sheltered
        F_sky, F_sur, T_sur = 0.0, 1.0, c.T_amb           # sees only objects at air temperature
        θ, G_dir, G_dif = NaN, 0.0, 0.0
    else
        F_sky = (1 + cosd(s.tilt)) / 2                    # roof 1, vertical wall 0.5, floor 0
        F_sur = 1 - F_sky                                 # the rest of the view is ground
        T_sur = c.T_ground
        # angle θ between the sun rays and the normal of the surface
        cosθ  = sind(c.sun_elevation) * cosd(s.tilt) +
                cosd(c.sun_elevation) * sind(s.tilt) * cosd(c.sun_azimuth - s.azimuth)
        θ     = c.DNI > 0 ? acosd(clamp(cosθ, -1, 1)) : NaN      # no meaning at night
        G_dir = c.DNI * max(cosθ, 0)                      # direct sun (0 if the sun is behind)
        G_dif = c.G_dif * F_sky                           # diffuse light comes from the sky dome
    end
    c.ground_radiation || (F_sur = 0.0)
    G = G_dir + G_dif                                     # irradiance on the surface [W/m²]

    Q_sun   = w.α_sol * A * G
    τ_direct = G_dir > 0 ? w.τ_sol * angular_factor(θ) : 0.0      # glass only (τ_sol = 0 for walls)
    τ_diffuse = w.τ_sol * angular_factor(60.0)                    # diffuse light: average angle ≈ 60°
    Q_trans = A * (τ_direct * G_dir + τ_diffuse * G_dif)

    # --- the heat flows as functions of the temperatures ---
    Q_conv_o(T_e, h_o) = h_o * A * (c.T_amb - T_e)
    Q_rad(T_e) = w.ε_ir * A * σ_SB * (F_sky * (amb.T_sky^4 - T_e^4) + F_sur * (T_sur^4 - T_e^4))
    Q_k(T_e, T_si)     = (T_e - T_si) / R_k
    Q_i(T_si, h_i)     = h_i * A * (T_si - c.T_cab)

    # --- iteration (steps 1-4 of the docstring) ---
    T_e, T_si = c.T_amb + 1.0, c.T_cab + 1.0              # 1. first guess
    converged = false
    for _ in 1:500
        h_o = outside_convection(s, T_e, c).h             # 2. freeze h
        h_i = inside_convection(s, T_si, c).h
        R_i  = 1 / (h_i * A)                              #    inside convection resistance [K/W]
        # 3. inside node: (T_e - T_si)/R_k = (T_si - T_cab)/R_i  ->  Q_k = (T_e - T_cab)/(R_k + R_i)
        balance(T) = Q_sun + Q_conv_o(T, h_o) + Q_rad(T) - (T - c.T_cab) / (R_k + R_i)
        T_e_new  = bisection(balance, 100.0, 700.0)
        T_si_new = c.T_cab + (T_e_new - c.T_cab) * R_i / (R_k + R_i)
        change   = max(abs(T_e_new - T_e), abs(T_si_new - T_si))
        T_e  += 0.5 * (T_e_new - T_e)                     # 4. move half-way
        T_si += 0.5 * (T_si_new - T_si)
        if change < 1e-9
            converged = true
            break
        end
    end
    converged || @warn "surface $(s.name): iteration did not converge"

    # --- final values, and the proof that both balances close (residuals ≈ 0) ---
    out = outside_convection(s, T_e, c)
    ins = inside_convection(s, T_si, c)
    q_conv, q_rad = Q_conv_o(T_e, out.h), Q_rad(T_e)
    q_k, q_i      = Q_k(T_e, T_si), Q_i(T_si, ins.h)

    return (name = s.name, A = A, θ = θ, G = G, F_sky = F_sky, F_sur = F_sur,
            T_e = T_e, T_si = T_si, outside = out, inside = ins,
            Q_sun = Q_sun, Q_conv_o = q_conv, Q_rad = q_rad, Q_k = q_k, Q_i = q_i,
            Q_trans = Q_trans,
            Q_cabin = q_i + Q_trans,                      # total heat this surface gives to the cabin
            residual_outside = Q_sun + q_conv + q_rad - q_k,
            residual_inside  = q_k - q_i,
            # resistances for the electrical analogue [K/W]
            R_conv_o = 1 / (out.h * A), R_k = R_k, R_conv_i = 1 / (ins.h * A),
            # overall coefficient air-to-air, convection + conduction only [W/m²K]
            U = 1 / (1 / out.h + R_cond(w) + 1 / ins.h))
end

# ------------------------------------------------------------------------------------------
# 7. Heat released by the occupants (PW5)
# ------------------------------------------------------------------------------------------

"""
    occupant_heat(person, T_air, T_rad)

Sensible and latent heat [W] released by one person (method of PW5):

    A_D     = 0.202·m^0.425·l^0.725                    Dubois body area [m²]
    t_sk    = 35.7 - 0.0275·M                          skin temperature [°C]   (work W = 0)
    t_o     = (t_air + t_rad)/2                        operative temperature [°C]
    q_skin  = (t_sk - t_o)/(R_cl + 1/(f_cl·h))         sensible heat through the clothing
    q_resp  = 0.0014·M·(34 - t_air)                    sensible heat by breathing
    q_lat   = M - q_skin - q_resp                      the rest leaves as water vapour
"""
function occupant_heat(person::Person, T_air, T_rad)
    t_a, t_r = in_degC(T_air), in_degC(T_rad)
    A_D   = 0.202 * person.mass^0.425 * person.height^0.725
    t_sk  = 35.7 - 0.0275 * person.M
    t_o   = (t_a + t_r) / 2
    R_cl  = 0.155 * (person.I_cl + (person.seated ? 0.15 : 0.0))     # [m²K/W]
    q_skin = (t_sk - t_o) / (R_cl + 1 / (person.f_cl * person.h))    # [W/m²]
    q_resp = 0.0014 * person.M * (34 - t_a)                          # [W/m²]
    q_sens = q_skin + q_resp
    q_lat  = max(person.M - q_sens, 0.0)
    return (sensible = q_sens * A_D, latent = q_lat * A_D, A_D = A_D, t_sk = t_sk, t_o = t_o)
end

# ------------------------------------------------------------------------------------------
# 8. Cabin balance, ventilation, recirculation, heat exchanger
# ------------------------------------------------------------------------------------------

"""
    air_handling(c, amb, ṁ_fresh, ṁ_total, T_vent, ṁ_w)

State of the air around the loop  outside -> mixing box -> heat exchanger -> cabin.
ṁ_total is blown into the cabin at T_vent; ṁ_fresh of it is outside air, the rest
(ṁ_total - ṁ_fresh) is cabin air that is recirculated. ṁ_w = water vapour from the crew.

Water balance. Without condensation all the water of the crew leaves with the fresh air:
    x_cab = x_out + ṁ_w/ṁ_fresh
If the mixed air holds more water than saturated air at T_vent can hold, the surplus
condenses on the cooling coil (ideal coil: air leaves saturated at T_vent). Then:
    x_supply = x_sat(T_vent)   and   x_cab = x_supply + ṁ_w/ṁ_total

Heat exchanger power (positive = heating, negative = cooling):
    Q_hx = ṁ_total·[h(T_vent, x_supply) - h_mix]
"""
function air_handling(c::Case, amb, ṁ_fresh, ṁ_total, T_vent, ṁ_w)
    ṁ_rec = ṁ_total - ṁ_fresh
    x_sat = x_saturation(T_vent, c.p)

    x_cab = amb.x + ṁ_w / ṁ_fresh                              # try: no condensation
    x_mix = (ṁ_fresh * amb.x + ṁ_rec * x_cab) / ṁ_total
    condensing = x_mix > x_sat
    if condensing
        x_sup = x_sat
        x_cab = x_sup + ṁ_w / ṁ_total
        x_mix = (ṁ_fresh * amb.x + ṁ_rec * x_cab) / ṁ_total
    else
        x_sup = x_mix
    end

    h_cab = enthalpy(c.T_cab, x_cab)
    h_mix = (ṁ_fresh * amb.h + ṁ_rec * h_cab) / ṁ_total       # mixing of two air streams
    h_sup = enthalpy(T_vent, x_sup)
    return (ṁ_rec = ṁ_rec, x_mix = x_mix, x_sup = x_sup, x_cab = x_cab,
            T_mix = temperature_from_enthalpy(h_mix, x_mix),
            h_mix = h_mix, h_sup = h_sup, h_cab = h_cab,
            condensing = condensing, ṁ_condensate = ṁ_total * (x_mix - x_sup),
            RH_cab = relative_humidity(c.T_cab, x_cab, c.p),
            Q_hx = ṁ_total * (h_sup - h_mix))
end

"""
    ventilation(c, amb, Q_walls, P_sens, P_lat, V_fresh)

Cabin energy balance (PW4):   ṁ·cp·(T_vent - T_cab) + P_sens + Q_walls = 0
(the latent heat does not appear: it leaves as water vapour, see the water balance).

 g) With the hygienic minimum of fresh air `V_fresh` [m³/s]: which T_vent is needed?
 h) If that T_vent is outside the allowed limits: fix T_vent at the limit, compute the
    (larger) flow that is then needed, and recirculate everything above the fresh-air minimum.
 6) Power of the heat exchanger.
"""
function ventilation(c::Case, amb, Q_walls, P_sens, P_lat, V_fresh)
    Q_load  = Q_walls + P_sens              # sensible heat entering the cabin air [W]
    ṁ_fresh = V_fresh / amb.v_da            # [kg dry air/s], volume measured at outside conditions
    ṁ_w     = P_lat / h_fg0                 # water vapour released by the crew [kg/s]

    # g) hygienic flow only
    cp_out     = cp_da + cp_wv * amb.x      # heat capacity of moist air per kg dry air [J/kgK]
    T_vent_hyg = c.T_cab - Q_load / (ṁ_fresh * cp_out)
    limited    = !(c.T_vent_min <= T_vent_hyg <= c.T_vent_max)

    # h) limit the supply temperature and increase the flow if necessary.
    #    cp depends (a little) on the humidity of the supply air, which is only known once
    #    the flows are known -> repeat a few times until x_supply does not change any more.
    T_vent  = clamp(T_vent_hyg, c.T_vent_min, c.T_vent_max)
    ṁ_total = ṁ_fresh
    state   = air_handling(c, amb, ṁ_fresh, ṁ_total, T_vent, ṁ_w)
    for _ in 1:50
        cp = cp_da + cp_wv * state.x_sup
        if limited
            ṁ_total = max(ṁ_fresh, -Q_load / (cp * (T_vent - c.T_cab)))
        else
            T_vent = c.T_cab - Q_load / (ṁ_fresh * cp)
        end
        new_state = air_handling(c, amb, ṁ_fresh, ṁ_total, T_vent, ṁ_w)
        done  = abs(new_state.x_sup - state.x_sup) < 1e-13
        state = new_state
        done && break
    end

    # for comparison: same flow and supply temperature but WITHOUT recirculation (all outside air)
    x_sup_fresh    = min(amb.x, x_saturation(T_vent, c.p))
    Q_hx_all_fresh = ṁ_total * (enthalpy(T_vent, x_sup_fresh) - amb.h)

    # proof: enthalpy balance of the cabin air (must be ≈ 0)
    residual = ṁ_total * (state.h_sup - state.h_cab) + Q_walls + P_sens +
               ṁ_w * (h_fg0 + cp_wv * in_degC(c.T_cab))

    return merge(state, (Q_load = Q_load, ṁ_fresh = ṁ_fresh, ṁ_w = ṁ_w, V_fresh = V_fresh,
                         T_vent_hyg = T_vent_hyg, limited = limited,
                         T_vent = T_vent, ṁ_total = ṁ_total,
                         V_supply = ṁ_total * volume_per_kg_dry_air(T_vent, state.x_sup, c.p),
                         recirculation = state.ṁ_rec / ṁ_total,
                         Q_hx_all_fresh = Q_hx_all_fresh, residual = residual))
end

# ------------------------------------------------------------------------------------------
# 9. Solve a complete case
# ------------------------------------------------------------------------------------------

"""
    solve_case(surfaces, c)

1. ambient conditions -> 2. every surface on its own (T_cab is imposed, so the surfaces do
not influence each other) -> 3. occupants -> 4. cabin balance and ventilation.
"""
function solve_case(surfaces::Vector{Surface}, c::Case)
    amb   = ambient(c)
    walls = [solve_surface(s, c, amb) for s in surfaces]
    Q_walls = sum(w.Q_cabin for w in walls)

    # mean radiant temperature felt by the crew = area-weighted inside surface temperature
    T_rad = sum(w.A * w.T_si for w in walls) / sum(w.A for w in walls)
    crew  = [occupant_heat(person, c.T_cab, T_rad) for person in c.crew]
    P_sens_crew = sum(o.sensible for o in crew; init = 0.0)
    P_lat       = sum(o.latent for o in crew; init = 0.0)
    P_sens      = P_sens_crew + c.P_equipment

    V_fresh = length(c.crew) * c.V_fresh_per_person
    vent    = ventilation(c, amb, Q_walls, P_sens, P_lat, V_fresh)

    return (case = c, amb = amb, walls = walls, Q_walls = Q_walls, T_rad = T_rad, crew = crew,
            P_sens_crew = P_sens_crew, P_sens = P_sens, P_lat = P_lat, vent = vent)
end

# ------------------------------------------------------------------------------------------
# 10. Printing
# ------------------------------------------------------------------------------------------

function print_surfaces(walls)
    println("\n  Temperatures and coefficients per surface")
    @printf("  %-22s %6s %6s %7s %7s %7s %7s %7s %7s\n",
            "surface", "A", "theta", "G_inc", "T_e", "T_si", "h_o", "h_i", "U")
    @printf("  %-22s %6s %6s %7s %7s %7s %7s %7s %7s\n",
            "", "m2", "deg", "W/m2", "degC", "degC", "W/m2K", "W/m2K", "W/m2K")
    for w in walls
        @printf("  %-22s %6.2f %6.1f %7.1f %7.2f %7.2f %7.2f %7.2f %7.2f\n",
                w.name, w.A, w.θ, w.G, in_degC(w.T_e), in_degC(w.T_si), w.outside.h, w.inside.h, w.U)
    end

    println("\n  Heat flows per surface [W]   (positive = towards the cabin)")
    @printf("  %-22s %8s %8s %8s %8s %8s %8s %8s %9s\n",
            "surface", "Q_sun", "Q_conv_o", "Q_rad", "Q_k", "Q_i", "Q_trans", "Q_cabin", "residual")
    for w in walls
        @printf("  %-22s %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %9.1e\n",
                w.name, w.Q_sun, w.Q_conv_o, w.Q_rad, w.Q_k, w.Q_i, w.Q_trans, w.Q_cabin,
                max(abs(w.residual_outside), abs(w.residual_inside)))
    end
    @printf("  %-22s %8.1f %8s %8s %8s %8.1f %8.1f %8.1f\n", "TOTAL",
            sum(w.Q_sun for w in walls), "", "", "",
            sum(w.Q_i for w in walls), sum(w.Q_trans for w in walls), sum(w.Q_cabin for w in walls))

    println("\n  Justification of the convection coefficients   (ok = inside validity range)")
    @printf("  %-22s %-4s %-8s %-44s %5s %9s %9s %8s %7s %3s\n",
            "surface", "side", "regime", "correlation", "L[m]", "Re", "Ra", "Ri", "Nu", "ok")
    for w in walls, (side, cv) in (("out", w.outside), ("in", w.inside))
        @printf("  %-22s %-4s %-8s %-44s %5.2f %9.2e %9.2e %8.1e %7.1f %3s\n",
                w.name, side, cv.mode, cv.corr, cv.L, cv.Re, cv.Ra, cv.Ri, cv.Nu, cv.valid ? "yes" : "NO")
    end

    println("\n  Resistances for the electrical analogue [K/W]")
    @printf("  %-22s %10s %10s %10s\n", "surface", "R_conv_o", "R_cond", "R_conv_i")
    for w in walls
        @printf("  %-22s %10.4f %10.4f %10.4f\n", w.name, w.R_conv_o, w.R_k, w.R_conv_i)
    end
end

function print_ventilation(c::Case, amb, v)
    println("\n  Cabin balance and ventilation")
    @printf("    sensible load on the cabin air          %9.1f W\n", v.Q_load)
    @printf("    fresh air (hygienic minimum)            %9.4f kg/s  = %.2f m3/min\n", v.ṁ_fresh, v.V_fresh * 60)
    @printf("  g) T_vent needed with the fresh air only  %9.1f degC  (allowed: %.0f to %.0f degC) -> %s\n",
            in_degC(v.T_vent_hyg), in_degC(c.T_vent_min), in_degC(c.T_vent_max),
            v.limited ? "NOT acceptable" : "acceptable")
    @printf("  h) supply temperature used                %9.1f degC\n", in_degC(v.T_vent))
    @printf("    total supply flow                       %9.4f kg/s  = %.2f m3/min\n", v.ṁ_total, v.V_supply * 60)
    @printf("    recirculated flow                       %9.4f kg/s  = %.0f %% of the supply\n", v.ṁ_rec, 100 * v.recirculation)
    @printf("    mixed air before the heat exchanger     %9.1f degC, x = %.2f g/kg\n", in_degC(v.T_mix), 1000 * v.x_mix)
    @printf("    supply air                              %9.1f degC, x = %.2f g/kg\n", in_degC(v.T_vent), 1000 * v.x_sup)
    @printf("    condensation on the coil                %9.3f kg/h\n", 3600 * v.ṁ_condensate)
    @printf("    cabin air                               %9.1f degC, x = %.2f g/kg, RH = %.0f %%\n",
            in_degC(c.T_cab), 1000 * v.x_cab, 100 * v.RH_cab)
    @printf("  6) HEAT EXCHANGER POWER                   %9.1f W  (%s)\n", v.Q_hx, v.Q_hx >= 0 ? "heating" : "cooling")
    @printf("    same supply but without recirculation   %9.1f W\n", v.Q_hx_all_fresh)
    @printf("    check: residual of the cabin balance    %9.1e W\n", v.residual)
end

function print_report(r)
    c, amb = r.case, r.amb
    println("\n", "="^120)
    println(" ", c.name)
    println("="^120)
    @printf("  Outside : %.1f degC, RH %.0f %%, x = %.2f g/kg, dew point %.1f degC, speed %.0f km/h\n",
            in_degC(c.T_amb), 100 * c.RH_amb, 1000 * amb.x, in_degC(amb.T_dp), 3.6 * c.v)
    @printf("            sky %.1f degC, ground %.1f degC, direct sun %.0f W/m2 at %.1f deg elevation, diffuse %.0f W/m2\n",
            in_degC(amb.T_sky), in_degC(c.T_ground), c.DNI, c.sun_elevation, c.G_dif)
    @printf("  Cabin   : %.1f degC required, %d occupants, equipment %.0f W\n",
            in_degC(c.T_cab), length(c.crew), c.P_equipment)

    print_surfaces(r.walls)

    println("\n  Occupants (PW5 method), mean radiant temperature ", @sprintf("%.1f degC", in_degC(r.T_rad)))
    for (i, o) in enumerate(r.crew)
        @printf("    person %d: M = %5.0f W/m2, sensible %6.1f W, latent %6.1f W\n",
                i, c.crew[i].M, o.sensible, o.latent)
    end
    @printf("    crew sensible %.1f W + equipment %.1f W = %.1f W ; crew latent %.1f W\n",
            r.P_sens_crew, c.P_equipment, r.P_sens, r.P_lat)

    print_ventilation(c, amb, r.vent)
    return nothing
end
