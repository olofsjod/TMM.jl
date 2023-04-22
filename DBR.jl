include("TMM.jl")

function calc_TR_vs_periods(p, λ)
    n = Complex{Float64}[];
    d = Float64[];
    wavelength = 1550;
    nSi = 3.30151480384;

    for m = 1:p
        push!(d, wavelength/4, wavelength/(4*nSi))
        push!(n, 1, nSi)
    end
    calc_TR(n, d, 0, λ)
end

# We will calculate the stopband for each period and plot it in a colormap
λmin = 500
λstep = 0.1
λmax = 4000
N = trunc(Int64, (λmax-λmin)/λstep) + 1
M = 100
Λ = λmin:λstep:λmax
Rmx = zeros(Float64, M, N);
Tmx = zeros(Float64, M, N);

for p in 1:M
    TR = [calc_TR_vs_periods(p, λ) for λ in λmin:λstep:λmax];
    T = first.(TR);
    R = last.(TR);

    Tmx[p, :] = T;
    Rmx[p, :] = R;
end

using PyPlot

pcolormesh(Λ, 1:M, Rmx)
colorbar()

