## Transfer Matrix Method implemented in Julia (TMM.jl)
# written by Olof Sjödin <olsj@kth.se> 2023

# This file includes a couple of functions in order to calculate the
# reflectance and transmittance of multi layer dielectric stack.

"""
    r_ip(n0, n1, φ0, φ1)

Computes the interface matrix for a p-polarised wave.
"""
function r_ip(n0, n1, φ0, φ1)
    Rp = ( n1*cosd(φ0) - n0*cosd(φ1) ) / ( n1*cosd(φ0) + n0*cosd(φ1) );
    Tp = 2*n0*cosd(φ0) / (n1*cosd(φ0) + n0*cosd(φ1));

    1/Tp.*[ 1 Rp; Rp 1 ]
end

"""
    r_is(n0, n1, φ0, φ1)

Computes the interface matrix for a s-polarised wave.
"""
function r_is(n0, n1, φ0, φ1)
    Rs = ( n0*cosd(φ0) - n1*cosd(φ1) ) / ( n0*cosd(φ0) + n1*cosd(φ1) );
    Ts = 2*n0*cosd(φ0) / (n0*cosd(φ0) + n1*cosd(φ1));

    1/Ts.*[ 1 Rs; Rs 1 ]
end

"""
    r_layer(n, d, φ, λ)

Computes the 2x2 propagation matrix in a dielectric with refractive index n and
thickness d. The wave is propagating the medium with incident angle φ and
wavelength λ.

"""
function r_layer(n, d, φ, λ)
    a = (2*pi*d*n/λ) * cosd(φ);

    [exp(im*a) 0 ; 0 exp(-im*a)]
end

@enum Polarisation begin
    Spol
    Ppol
end

"""
    calc_TR(n::Vector, d::Vector, φ0::Real, λ::Real)

Computes the transmission (reflection) coefficient through (from) a multi layer
dielectric stack.

# Arguments

- n::Vector: the refractive indices (real numbers) of a multi-layer stack. n[1]
  corresponds to the first layer of the stack.
- d::Vector: the thicknesses of each layer of the multi-layer stack. d[1]
  corresponds to the first layer and has a unit of nanometers.
- φ0::Real: the incident angle in degrees of the impinging wave.
- λ::Real: the wavelength of the wave in nanometers.
- pol::Polarisation: the impinging wave's polarisation.
"""

function calc_TR(n::Vector, d::Vector, φ0::Real, λ::Real, pol::Polarisation)
    if pol == Spol
        r_i = r_is
    else
        r_i = r_ip
    end

    φ = Vector{Float64}(undef, length(n))
    φ[1] = φ0

    # does this work even though n is complex?
    for k=2:length(n)
        φ[k] = asind(n[k-1]*sind(φ[k-1])/n[k])
    end

    LxI = Vector{Matrix{ComplexF64}}(undef, (length(n)-1))
    LxI[1] = r_i(n[1], n[2], φ[1], φ[2])

    for k=2:(length(n)-1)
        LxI[k] = r_layer.(n[k], d[k], φ[k], λ).*r_i.(n[k], n[k+1], φ[k], φ[k+1])
    end

    S = foldl(*, LxI)

    T, R = abs2(1 / S[1,1]), abs2(S[2,1] / S[1,1])
end

using ArgParse

s = ArgParseSettings(description = "Calculates the transmittance and reflectance
of a dielectric slab.")

@add_arg_table! s begin
    "--wavelength", "-w"
    help = "The wavelength (in nanometers) if doing an angle sweep. Wavelength span (start:step:end) in nanometers. Example: 1550, 1500:0.1:1600"
    required = true
    "--angle", "-a"
    help = "The incident angle range if doing an angle sweep. Or incident angle if doing a wavelength sweep. Example: 0:0.1:89.9, 10"
    required = true
    "--polarisation", "-p"
    help = "The polarisation of the incident wave (S/P or TE/TM). Possible values: TE, TM, S, P"
    required = true
    "-d"
    help = "Slab thickness in terms of nanometers."
    nargs = '+'
    required = true
    "-n"
    help = "Slab refractive index."
    nargs = '+'
    required = true
end

parsed_args = parse_args(ARGS, s)
println("Parsed args:")

for (key, val) in parsed_args
    println(" $key => $(repr(val))")
end

# Check each argument to detect whether they are ranges or not.
splitted_wvl = split(parsed_args["wavelength"], ":")
splitted_ang = split(parsed_args["angle"], ":")

# Sanity check
if length(splitted_wvl)*length(splitted_ang) != 3
    error("One of the arguments must be an interval! Such as 1500:1:1600, 1300:2:1800, etc.")
end

@enum PlotVS AngleScan WavelengthScan

if length(splitted_wvl) == 3
    scan_type = WavelengthScan

    wvl_min = parse(Float64, splitted_wvl[1])
    wvl_step = parse(Float64, splitted_wvl[2])
    wvl_max = parse(Float64, splitted_wvl[3])

    angle = parse(Float64, splitted_ang[1])

    if !(wvl_min < wvl_max && (wvl_max-wvl_min) >= wvl_step)
        error("The wavelength interval specified is not correct.")
    end

elseif length(splitted_ang) == 3
    scan_type = AngleScan

    ang_min = parse(Float64, splitted_ang[1])
    ang_step = parse(Float64, splitted_ang[2])
    ang_max = parse(Float64, splitted_ang[3])

    wavelength = parse(Float64, splitted_wvl[1])

    if !(ang_min < ang_max && (ang_max-ang_min) >= ang_step)
        error("The angle interval specified is not correct.")
    end
else
    error("Unknown error.")
end

if uppercase(parsed_args["polarisation"]) == "TE" || uppercase(parsed_args["polarisation"]) == "S"
    pol_type = Spol
elseif uppercase(parsed_args["polarisation"]) == "TM" || uppercase(parsed_args["polarisation"]) == "P"
    pol_type = Ppol
end

using Plots

# Install with
# using Pkg
# Pkg.add("Plots")

## It is neccessary to do this in order to install matplotlib automatically.
# using Pkg
# ENV["PYTHON"]=""
# Pkg.build("PyCall")
#
# Then you can install PyPlot.jl with
# Pkg.add("PyPlot")

# To use the PyPlot backend in Plots
pyplot()

nArr = parse.(Float64, parsed_args["n"])
dArr = parse.(Float64, parsed_args["d"])

# Create the structure
n = vcat([1], nArr, [1])
d = vcat([0], dArr, [0]);

if scan_type == WavelengthScan
    N = trunc(Int, (wvl_max-wvl_min)/wvl_step)

    wavelength = LinRange(wvl_min, wvl_max, N)
    T = zeros(N)
    R = zeros(N)

    @time begin
        Threads.@threads for k = 1:N
            Tc, Rc = calc_TR(n, d, angle, wavelength[k], pol_type)

            T[k] = Tc
            R[k] = Rc
        end
    end

    plot()
    plot!(wavelength, T, label="Transmittance")
    plot!(wavelength, R, label="Reflectance")
    gui()
    read(stdin, Char)

elseif scan_type == AngleScan
    N = trunc(Int, (ang_max-ang_min)/ang_step)

    angle = LinRange(ang_min, ang_max, N)
    T = zeros(N)
    R = zeros(N)

    Threads.@threads for k = 1:N
        Tc, Rc = calc_TR(n, d, angle[k], wavelength, pol_type)

        T[k] = Tc
        R[k] = Rc
    end

    plot()
    plot!(angle, T, label="Transmittance")
    plot!(angle, R, label="Reflectance")
    gui()
    read(stdin, Char)

end
