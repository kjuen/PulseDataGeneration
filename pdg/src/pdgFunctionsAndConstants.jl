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


# Trebino Gl. (2.17)
function Sw2Sl(ww, Sw)
    @assert length(ww) == length(Sw)
    @assert issorted(ww)
    ll =  c2p ./ reverse(ww)
    llGrid = range(ll[1], ll[end], length(ll))   # äquidistant
    Sl = reverse(Sw .* ww.^2 / c2p)
    itp = linear_interpolation(ll, Sl)
    SlGrid = itp.(llGrid)

    return (llGrid, SlGrid)
end

# Trebino Gl. (2.17)
function Sl2Sw(ll, Sl)
    @assert length(ll) == length(Sl)
    @assert issorted(ll)

    ww =  c2p ./ reverse(ll)
    @assert issorted(ww)
    wwGrid = range(ww[1], ww[end], length(ww))

    Sw = reverse(Sl .* ll.^2 / c2p)
    itp = linear_interpolation(ww, Sw)
    SwGrid = itp.(wwGrid)

    return (wwGrid, SwGrid)
end

function intensityMatFreq2Wavelength(ww, freqIntMat)
    @assert length(ww) == size(freqIntMat, 2)
    @assert issorted(ww)
    wvIntMat = zeros(eltype(freqIntMat), size(freqIntMat))


    ll =  c2p ./ reverse(ww)   # nicht äquidistant!
    # @assert issorted(ll)
    llGrid = range(ll[1], ll[end], length(ll))   # äquidistant

    for (i, Sw) in enumerate(eachrow(freqIntMat))
        Sl = reverse(Sw .* ww.^2 / c2p)
        itp = linear_interpolation(ll, Sl)
        wvIntMat[i,:] = itp.(llGrid)
    end

    return (llGrid, wvIntMat)
end

function intensityMatWavelength2Freq(ll, wvIntMat)
    @assert length(ll) == size(wvIntMat, 2)
    @assert issorted(ll)
    freqIntMat = zeros(eltype(wvIntMat), size(wvIntMat))

    ww =  c2p ./ reverse(ll)
    @assert issorted(ww)
    wwGrid = range(ww[1], ww[end], length(ww))

    for (i, Sl) in enumerate(eachrow(wvIntMat))
        Sw = reverse(Sl .* ll.^2 / c2p)
        itp = linear_interpolation(ww, Sw)
        freqIntMat[i,:] = itp.(wwGrid)
    end

    return (wwGrid, freqIntMat)
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
    return (ret, ((-N+1):(N))*Ts)
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

#* Rana

"""
 Rana-Algorithm from https://arxiv.org/abs/1811.01470
 # Arguments:
 - `ScSw`: Auto convolution of the spectrum S(w)
 - `w0`: frequency resolution
 - `alpha`: parameter controlling continuity
 - `beta`: parameter controlling first derivative
 - `gamma`: parameter controlling second derivative
"""
function rana(ScSw, w0; alpha=0.09, beta=0.425, gamma=1.0)
    N = length(ScSw)
    # We assume that ScSw lives in the centered baseband.
    # st2c is centered around t=0
    st2c = fftshift(ifft(fftshift(ScSw)))  # Equation (4) in Rana-Paper
    st2c ./= maximum(abs.(st2c))

    # Take square root of st2 accordung to RANA paper
    st = similar(st2c)
    st[1] = sqrt(st2c[1])
    for i = 2:N
        sp = sqrt(st2c[i])
        sm = - sp
        # Equation (6) in Rana-Paper:
        Deltap = alpha * abs(sp - st[i-1])
        Deltam = alpha * abs(sm - st[i-1])

        # Equation (7) in Rana-Paper:
        if i > 2
            Deltap += beta * abs(sp - 2*st[i-1] + st[i-2])
            Deltam += beta * abs(sm - 2*st[i-1] + st[i-2])
        end

        # Equation (8) in Rana-Paper:
        if i > 3
            Deltap += gamma * abs(sp - 3*st[i-1] + 3 * st[i-2] - st[i-3])
            Deltam += gamma * abs(sm - 3*st[i-1] + 3 * st[i-2] - st[i-3])
        end

        # Take closer root, see bottom of page 4 in Rana-Paper
        if Deltap < Deltam
            st[i] = sp
        else
            st[i] = sm
        end
    end

    # Go back to frequency domain
    Sw = fftshift(fft(fftshift(st)))   # this shoud be real !!
    Sw ./= maximum(abs.(Sw))

    @show maximum(imag.(Sw)) / maximum(real.(Sw))
    ws = N*w0
    wAxis = (0:(N-1))*w0 .- ws/2
    return (Sw, wAxis)
end


function ranaDetails(ScSw, w0; alpha=0.09, beta=0.425, gamma=1.0)
    N = length(ScSw)
    # We assume that ScSw lives in the centered baseband.
    # st2c is centered around t=0
    st2c = fftshift(ifft(fftshift(ScSw)))  # Equation (4) in Rana-Paper
    st2c ./= maximum(abs.(st2c))

    # Take square root of st2 accordung to RANA paper
    stc = similar(st2c)
    stc[1] = sqrt(st2c[1])
    for i = 2:N
        sp = sqrt(st2c[i])
        sm = - sp
        # Equation (6) in Rana-Paper:
        Deltap = alpha * abs(sp - stc[i-1])
        Deltam = alpha * abs(sm - stc[i-1])

        # Equation (7) in Rana-Paper:
        if i > 2
            Deltap += beta * abs(sp - 2*stc[i-1] + stc[i-2])
            Deltam += beta * abs(sm - 2*stc[i-1] + stc[i-2])
        end

        # Equation (8) in Rana-Paper:
        if i > 3
            Deltap += gamma * abs(sp - 3*stc[i-1] + 3 * stc[i-2] - stc[i-3])
            Deltam += gamma * abs(sm - 3*stc[i-1] + 3 * stc[i-2] - stc[i-3])
        end

        # Take closer root, see bottom of page 4 in Rana-Paper
        if Deltap < Deltam
            stc[i] = sp
        else
            stc[i] = sm
        end
    end

    # Go back to frequency domain
    Sw = fftshift(fft(fftshift(stc)))   # this shoud be real !!
    Sw ./= maximum(abs.(Sw))

    @show maximum(imag.(Sw)) / maximum(real.(Sw))
    ws = N*w0
    Ts = 2*pi/ws
    T0 = N*Ts
    ttc = (0:(N-1))*Ts .- T0/2
    wAxis = (0:(N-1))*w0 .- ws/2
    return (Sw, wAxis, st2c, stc, ttc)
end

"""
   Faltung so dass das Ausgangssignal die gleiche Länge wie das Eingangssignal hat.
   Außerdem Normierung auf Max = 1
"""
function myAutoConv(s)
    sf = fft(s)
    c =  abs.(fftshift(ifft(sf.^2)))
    c ./= maximum(abs.(c))
    return c
end


function subtractBackground(sig; thres=0.02)
    absSig = abs.(sig)
    thresVal = maximum(absSig)*thres
    bgVals = filter(s -> s<thresVal, absSig)
    bg = isempty(bgVals) ? 0.0 : sum(bgVals) / length(bgVals)
    return absSig .- bg
end

function rmsDiff(s, sRef)
    @assert length(s) == length(sRef)
    Z = sqrt(sum(sRef.^2))
    return 1/Z * sqrt(sum((s .- sRef).^2))
end


function ranaImproved(ScSw, w0; alpha=0.09, beta=0.425, gamma=1.0)
    N = length(ScSw)
    @assert iseven(N)
    # We assume that ScSw lives in the centered baseband.
    # st2c is centered around t=0
    st2c = fftshift(ifft(fftshift(ScSw)))  # Equation (4) in Rana-Paper
    st2c ./= maximum(abs.(st2c))

    # Take square root of st2 accordung to RANA paper
    stc = similar(st2c)
    stc[1] = sqrt(st2c[1])
    for i = 2:N
        sp = sqrt(st2c[i])
        sm = - sp
        # Equation (6) in Rana-Paper:
        Deltap = alpha * abs(sp - stc[i-1])
        Deltam = alpha * abs(sm - stc[i-1])

        # Equation (7) in Rana-Paper:
        if i > 2
            Deltap += beta * abs(sp - 2*stc[i-1] + stc[i-2])
            Deltam += beta * abs(sm - 2*stc[i-1] + stc[i-2])
        end

        # Equation (8) in Rana-Paper:
        if i > 3
            Deltap += gamma * abs(sp - 3*stc[i-1] + 3 * stc[i-2] - stc[i-3])
            Deltam += gamma * abs(sm - 3*stc[i-1] + 3 * stc[i-2] - stc[i-3])
        end

        # Take closer root, see bottom of page 4 in Rana-Paper
        if Deltap < Deltam
            stc[i] = sp
        else
            stc[i] = sm
        end
    end
    # Go back to frequency domain
    SwSimple = fftshift(fft(fftshift(stc)))
    SwSimple ./= maximum(abs.(SwSimple))
    # Take magnitude und subtract backgound
    SwSimple = subtractBackground(SwSimple)
    ScSwSimple = myAutoConv(SwSimple)
    rmsDiffSimple = rmsDiff(ScSwSimple, ScSw)

    # Symmetrisiing of stc:
    middleIndex = N÷2 + 1
    stcSym = copy(stc)
    stcSym[(middleIndex+1):end] = conj.(reverse(stcSym[2:(middleIndex-1)]))
    SwSym = fftshift(fft(fftshift(stcSym)))
    SwSym ./= maximum(abs.(SwSym))
    # Take magnitude und subtract backgound
    SwSym = subtractBackground(SwSym)
    ScSwSym = myAutoConv(SwSym)
    rmsDiffSym = rmsDiff(ScSwSym, ScSw)

    # Loop just over the right hand side,  starting in the middle
    stcRight = similar(st2c)
    stcRight[1] = sqrt(st2c[1])  # this one remains undetermined
    stcRight[middleIndex] = sqrt(st2c[middleIndex])
    for ii = 1:(N÷2 - 1)
        i = ii + middleIndex
        sp = sqrt(st2c[i])
        sm = - sp
        # Equation (6) in Rana-Paper:
        Deltap = alpha * abs(sp - stcRight[i-1])
        Deltam = alpha * abs(sm - stcRight[i-1])

        # Equation (7) in Rana-Paper:
        if ii > 1
            Deltap += beta * abs(sp - 2*stcRight[i-1] + stcRight[i-2])
            Deltam += beta * abs(sm - 2*stcRight[i-1] + stcRight[i-2])
        end

        # Equation (8) in Rana-Paper:
        if ii > 2
            Deltap += gamma * abs(sp - 3*stcRight[i-1] + 3 * stcRight[i-2] - stcRight[i-3])
            Deltam += gamma * abs(sm - 3*stcRight[i-1] + 3 * stcRight[i-2] - stcRight[i-3])
        end

        # Take closer root, see bottom of page 4 in Rana-Paper
        if Deltap < Deltam
            stcRight[i] = sp
        else
            stcRight[i] = sm
        end
    end
    # fill left part by mirroring the right part
    stcRight[2:(middleIndex-1)] = conj.(reverse(stcRight[(middleIndex+1):end]))

    # Go back to frequency domain
    SwRight = fftshift(fft(fftshift(stcRight)))
    SwRight ./= maximum(abs.(SwRight))
    # Take magnitude und subtract backgound
    SwRight = subtractBackground(SwRight)
    ScSwRight = myAutoConv(SwRight)
    rmsDiffRight = rmsDiff(ScSwRight, ScSw)


     # Loop just over the left hand side,  starting in the middle
    stcLeft = similar(st2c)
    stcLeft[middleIndex] = sqrt(st2c[middleIndex])
    for ii = 1:(N÷2)
        i = middleIndex - ii
        sp = sqrt(st2c[i])
        sm = - sp
        # Equation (6) in Rana-Paper:
        Deltap = alpha * abs(sp - stcLeft[i+1])
        Deltam = alpha * abs(sm - stcLeft[i+1])

        # Equation (7) in Rana-Paper:
        if ii > 1
            Deltap += beta * abs(sp - 2*stcLeft[i+1] + stcLeft[i+2])
            Deltam += beta * abs(sm - 2*stcLeft[i+1] + stcLeft[i+2])
        end

        # Equation (8) in Rana-Paper:
        if ii > 2
            Deltap += gamma * abs(sp - 3*stcLeft[i+1] + 3 * stcLeft[i+2] - stcLeft[i+3])
            Deltam += gamma * abs(sm - 3*stcLeft[i+1] + 3 * stcLeft[i+2] - stcLeft[i+3])
        end

        # Take closer root, see bottom of page 4 in Rana-Paper
        if Deltap < Deltam
            stcLeft[i] = sp
        else
            stcLeft[i] = sm
        end
    end
    # fill right hand part by mirroring the left part
    stcLeft[(middleIndex+1):end] = conj.(reverse(stcLeft[2:(middleIndex-1)]))

    # Go back to frequency domain
    SwLeft = fftshift(fft(fftshift(stcLeft)))
    SwLeft ./= maximum(abs.(SwLeft))
    # Take magnitude und subtract backgound
    SwLeft = subtractBackground(SwLeft)# Take magnitude und subtract backgound
    ScSwLeft = myAutoConv(SwLeft)
    rmsDiffLeft = rmsDiff(ScSwLeft, ScSw)

    # So now we have four candicates: pick the best
    rmsVec = [rmsDiffSimple, rmsDiffSym, rmsDiffRight, rmsDiffLeft]
    #@show rmsVec
    (minRms, minIdx) = findmin(rmsVec)
    Sw = SwLeft
    if minIdx == 1
        Sw = SwSimple
    elseif minIdx == 2
        Sw = SwSym
    elseif minIdx == 3
        Sw = SwRight
    end

    # construct frequency axis the returned signal lives on
    ws = N*w0
    wAxis = (0:(N-1))*w0 .- ws/2

    return (Sw, wAxis, minIdx, minRms)
end




export berechneFWHM, berechneRMSBreite, berechneCOM, intensityMatFreq2Wavelength, intensityMatWavelength2Freq, Sw2Sl, Sl2Sw, specMatCorrectlySampled, autocorr, mseUpToScale, traceError, mseUpToSignAndLin, idxRangeAboveThres, idxRangeWithinLimits, rana, ranaDetails, ranaImproved
