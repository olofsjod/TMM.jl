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
