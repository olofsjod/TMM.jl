include("TMM.jl")

using ArgParse

s = ArgParseSettings(description = "Calculates the transmittance and reflectance
of a dielectric slab.")

@add_arg_table! s begin
    "wavelength"
    help = "The wavelength (in nanometers) if doing an angle sweep. Wavelength span (start:step:end) in nanometers. Example: 1550, 1500:0.1:1600"
    required = true
    "angle"
    help = "The incident angle range if doing an angle sweep. Or incident angle if doing a wavelength sweep. Example: 0:0.1:89.9, 10"
    required = true
    "polarisation"
    help = "The polarisation of the incident wave (S/P or TE/TM). Possible values: TE, TM, S, P"
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
elseif uppercase(parsed_args["polarisation"]) == "TE" || uppercase(parsed_args["polarisation"]) == "S"
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

# Create the structure
n = Complex{Float64}[];
d = Float64[];
wavelength = 1000;
nSiO2 = 1.75;

push!(d, 0)
push!(n, 1)

push!(d, 2000)
push!(n, nSiO2)

push!(d, 0)
push!(n, 1)

if scan_type == WavelengthScan
    N = trunc(Int, (wvl_max-wvl_min)/wvl_step)

    wavelength = LinRange(wvl_min, wvl_max, N)
    T = zeros(N)
    R = zeros(N)

    Threads.@threads for k = 1:N
        Tc, Rc = calc_TR(n, d, angle, wavelength[k], pol_type)

        T[k] = Tc
        R[k] = Rc
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
