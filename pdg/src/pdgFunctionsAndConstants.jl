# Pulse Data Generation

using Trapz
using Interpolations
using Statistics

#* constants
const femto = 1e-15
const nano = 1e-9
const mikro = 1e-6
const milli = 1e-3
const Kilo = 1e3
const Mega = 1e6
const Giga = 1e9
const Tera = 1e12
const angFregTHz = 2*pi * 1e12

const clight = 299_792_458    # Lichtgeschwindigkeit
const epsilon0 = 8.8541878128*1e-12
const c2p = 2*pi*clight
const sc2p = sqrt(c2p)

export femto, nano, mikro, milli, Kilo, Mega, Giga, Tera, angFregTHz, clight, epsilon0, c2p, sc2p



rd1(x) = round(x, digits=1)
rd2(x) = round(x, digits=2)

#* Small helper functions

# kleine Summary-Funktion für die Achsen-Vektoren
function mys(v)
    ex = extrema(v)
    println("$(length(v)) Elements from $(ex[1]) to $(ex[2])")
    m = mean(diff(v))
    s = std(diff(v))
    println("Increment: $(m) with cv = $(s/m)")
end


"""
  zeroCrossings: berechnet die Nulldurchgänge eines Vektors y
"""
function zeroCrossings(x::AbstractVector,
                       y::AbstractVector;
                       TOL = 1e-12)::Vector
    N = length(x)
    if length(y) != N
        error("Argumente müssen die gleiche Länge haben")
    end

    xzero = []

    if abs(y[1]) < TOL
        push!(xzero, x[1])
    end

    for n=2:N
        if abs(y[n]) < TOL
            push!(xzero, x[n])
            continue
        end
        if y[n]*y[n-1] < 0 && abs(y[n-1]) > TOL
            # fit a straight line bewteen the points and find the zero of this line
            Dx = x[n] - x[n-1]
            if !(Dx>0)
                error("Abstand der x-Werte muss > 0 sein")
            end
            Dy = y[n] - y[n-1]
            m = Dy / Dx
            @assert abs(m) > 0
            b1 = y[n-1] - m * x[n-1]
            b2 = y[n] - m * x[n]
            b = (b1 + b2) / 2
            push!(xzero, -b/m)
        end
    end
    return xzero
end


function berechneFWHM(yd::Vector{<:Real}, tt::AbstractVector)::Real
    (xm, im) = findmax(yd)
    xz = zeroCrossings(tt, yd .- xm/2)
    if length(xz) < 2
        println("FHWM kann nicht berechnet werden, da <2 HW-Punkte")
        # error("FHWM kann nicht berechnet werden, da <2 HW-Punkte")
        return -1.0
    else # mehr als 2
        ext = extrema(xz)
        return ext[2] - ext[1]
    end
end


function berechneRMSBreite(yd::Vector{<:Real},
                           tt::Union{Vector{<:Number}, AbstractRange{<:Number}})::Real
    N = length(yd)
    if length(tt) != N
        error("yd und tt müssen die gleiche Länge haben")
    end

    normFac = sum(abs.(yd).^2)

    mom1 = 1 / normFac * sum(tt .* abs.(yd).^2)
    mom2 = 1 / normFac * sum(tt.^2 .* abs.(yd).^2)

    md = mom2 - mom1^2
    if md < 0.0
        return -1.0
    else
        return sqrt(md)
    end
end


function berechneCOM(yd::Vector{<:Real},
                     tt::Union{Vector{<:Number}, AbstractRange{<:Number}})::Real
    N = length(yd)
    if length(tt) != N
        error("yd und tt müssen die gleiche Länge haben")
    end

    normFac = sum(abs.(yd).^2)

    mom1 = 1 / normFac * sum(tt .* abs.(yd).^2)
    return mom1
end




function mydiff5(y, h)
    # MYDIFF5 Numerische Ableitung per 5-Punkt-Formel
    N = length(y)
    k5 = 3:(N-2)
    dynum5 = (y[k5.-2] .- 8*y[k5.-1] .+ 8*y[k5.+1] .- y[k5.+2])/(12*h)
    return dynum5
end
# kleiner Test
# xx = range(-3, 3, 100)
# yy = 3*xx .* 4*xx.^2
# maximum(abs(mydiff5(yy, xx[2] - xx[1]) .- (3 .+ 8*xx[3:end-2])))  # e-13



# bestimme alle Indizes zwischen den äußersten Werten des Vektors 'sig',
# die über thres*maximum(sig) liehen
function idxRangeAboveThres(sig, thres)
    maxSig = maximum(abs.(sig))
    idxLims = extrema(findall(abs.(sig) .> thres * maxSig))
    return idxLims[1]:idxLims[2]
end

# bestimme den Index-Bereich, für den die Werte von sig innerhalb der limits liegen
function idxRangeWithinLimits(sig, limits)
    @assert issorted(sig)
    @assert limits[2] > limits[1]
    idx1 = findfirst(sig .> limits[1])
    idx2 = findlast(sig .< limits[2])
    return (idx1, idx2)
end


function resample(x::AbstractVector,
                  f::AbstractVector,
                  newx::AbstractVector)::AbstractVector
    itp = linear_interpolation(x, f)   # cubic doesn't work for irregular grids
    newf = itp.(newx)
    @assert length(newf) == length(newx)
    return newf
end


function regGridInterp(xIreg, funcVals; N=length(xIreg))
    @assert issorted(xIreg)
    @assert length(xIreg) == length(funcVals)
    x = range(xIreg[1], xIreg[end], N)
    # itp = linear_interpolation(xIreg, funcVals)
    # return (x, itp.(x))
    return (x, resample(xIreg, funcVals, x))
end



# Test von regGridInterp
# ff = 2.0:0.2:4.0
# funcf(f) = (f-2.6).^2
# ll = sort(1.0 ./ ff)
# funcl(l) = funcf(1/l)
# yy = funcl.(ll)
# f = Figure()
# ax1 = Axis(f[1, 1], title = "f-Achse", xlabel = "f", ylabel="func")
# scatterlines!(ax1, ff, funcf.(ff), label="", linewidth=2, linestyle=:solid) #, color=:blue;)
# ax2 = Axis(f[2, 1], title = "l-Achse", xlabel = "l", ylabel="func")
# scatterlines!(ax2, ll, yy, label="Orig", linewidth=2, linestyle=:solid)
# # Interpolation der Funktion auf lGrid-Werten, nur unter Nutzung der (ll, yy)-Punkte
# (lG, lf) = regGridInterp(ll, yy)
# (lG2, lf2) = regGridInterp(ll, yy; N=25)
# scatterlines!(ax2, lG, lf, label="Interp", linewidth=2, linestyle=:solid)
# scatterlines!(ax2, lG2, lf2, label="Interp2", linewidth=2, linestyle=:solid)
# axislegend(ax2, position=:lb)
# display(f)

function delayToSpecRes(delay, nSpec, wCenter)
    @assert iseven(nSpec)
    idxVec = (-nSpec ÷ 2) : (nSpec ÷ 2 - 1)

    ws = 2*pi / delay   # Abtastrate
    w0 = ws / nSpec            # Aufloesung im Frequenzbereich
    wAxis = idxVec * w0 .+ 2 * wCenter

    lamdaMin = c2p / maximum(wAxis)
    lamdaMax = c2p / minimum(wAxis)
    lambdaAxis = LinRange(lamdaMin, lamdaMax, nSpec)

    return (wAxis[2] - wAxis[1], lambdaAxis[2] - lambdaAxis[1])

end




function intensityMatFreq2Wavelength(ww, freqIntMat)
    @assert length(ww) == size(freqIntMat, 2)
    @assert issorted(ww)
    wvIntMat = zeros(eltype(freqIntMat), size(freqIntMat))


    ll =  c2p ./ reverse(ww)
    llGrid = range(ll[1], ll[end], length(ll))
    @assert issorted(ll)

    for (i, Sw) in enumerate(eachrow(freqIntMat))
        Sl = reverse(Sw .* ww.^2 / c2p)
        itp = linear_interpolation(ll, Sl)
        wvIntMat[i,:] = itp.(llGrid)
    end

    return (llGrid, wvIntMat)
end


function specMatCorrectlySampled(mat::Matrix{T}; nPix=3, thres=1e-4) where {T <: Real}
    low = thres * maximum(mat)
    #@show maximum(mat)
    #@show low
    max1 = maximum(mat[1:nPix,:])
    max2 = maximum(mat[:, 1:nPix])
    max3 = maximum(mat[end-nPix+1,:])
    max4 = maximum(mat[:, end-nPix+1])

    if max1 > low
        #println("Problem oben")
        return false
    end
    if max2 > low
        #println("Problem links")
        return false
    end
    if max3 > low
        #println("Problem unten")
        return false
    end
    if max4 > low
        #println("Problem rechts")
        return false
    end
    return true


    # return any(mat[1:nPix,:] .> low) ||
    #     any(mat[:, 1:nPix] .> low)  ||
    #     any(mat[end-nPix+1,:] .> low) ||
    #     any(mat[:, end-nPix+1] .> low)
end


function autocorr(x::AbstractVector; Ts =1)
    N = length(x)
    x1 = [zeros(eltype(x), N); x]
    x2 = [x; zeros(eltype(x), N)]
    ret = zeros(eltype(x), 2*N)
    t = 1:(2*N)
    for i in 1:(2*N)
        xs = circshift(x2, i)
        ret[i] = trapz(t, x1.*xs) * Ts
    end
    return (ret, (-N:(N-1))*Ts)
end
# kleiner Test


# tt = range(-5, 5, 100)
# Ts = tt[2] - tt[1]
# yy = exp.(-tt.^2)
# (ac, aci) = autocorr(yy; Ts=Ts)
# n = length(ac)

# f = lines(aci, ac , label="ac", linewidth=2, linestyle=:solid, #color=:blue;
#            axis = (; title = "Ein Plot", xlabel = "x", ylabel="y"))
# lines!(tt, yy, label="y", linewidth=2, linestyle=:solid)
# # lines!(1:n, ys, label="ys", linewidth=2, linestyle=:solid)
# axislegend(position=:lt)
# display(f)



#* Fehlerfunktionen

function mseUpToScale(A::AbstractArray, B::AbstractArray)
    mu = sum(A .* B) / sum(B.*B)
    return sum( (A .- mu * B).^2 )
end

"""
    trace error zweier Frog-Traces, siehe Gl. (11) und (12) im pypret-Paper.
    Das sollte das gleiche sein wie
"""
function traceError(Tmeas::Matrix{T}, Tref::Matrix{T})::T where {T <: Real}
    r = mseUpToScale(Tmeas, Tref)
    normFactor = prod(size(Tmeas)) * maximum(Tmeas)^2
    return sqrt(r/normFactor)
end



function mseUpToSignAndLin(v1, v2)
    N = length(v1)
    @assert length(v2) == N
    x = 1:N
    D = hcat(ones(N), x)

    y1 = v1 .- v2
    beta1 = D'*D \ D'*y1
    mse1 = 1/N*sum((D*beta1 .- y1).^2)
    y2 = v1 .+ v2
    beta2 = D'*D \ D'*y2
    mse2 = 1/N*sum((D*beta2 .- y2).^2)
    if mse1 < mse2
        return (mse1, D*beta1 .- y1)
    else
        return (mse2, D*beta2 .- y2)
    end
end


export berechneFWHM, berechneRMSBreite, berechneCOM, intensityMatFreq2Wavelength, specMatCorrectlySampled, autocorr, mseUpToScale, traceError, mseUpToSignAndLin, idxRangeAboveThres, idxRangeWithinLimits
