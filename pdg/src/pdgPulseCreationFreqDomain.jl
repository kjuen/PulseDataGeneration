# cd(@__DIR__)
# Pkg.activate("jlenvPulseCreation")


# module PulseCreation

import Trapz
using FFTW

# import Distributions
rdVal(lo, hi) = lo + (hi-lo)*rand()
rdSign() = (-1)^convert(Int, round(rand()))


# function cubicPhaseFunc(om, omPeak, timeShift, gdd, tod)
#     ph = om * timeShift .+ gdd/2 * (om .- omPeak).^2 .+ tod/6*(om .- omPeak).^3
#     gd = timeShift .+ gdd*(om .- omPeak) .+ tod/2 * (om .- omPeak).^2
#     return (ph, gd)
# end

function cubicPhaseFuncAngFreq(omPeak, timeShift, gdd, tod)
    ph(om) = (om .- omPeak) * timeShift .+ gdd/2 * (om .- omPeak).^2 .+ tod/6*(om .- omPeak).^3
    gd(om) = timeShift .+ gdd*(om .- omPeak) .+ tod/2 * (om .- omPeak).^2
    return (ph, gd)
end




function gaussShape(x, xPeak, rmsw, A)
    # gaussSahpe create gauss shape with
    # - x: Vector
    # - maximum at 'xPeak'
    # - given root mean square width
    # - given area 'A'

    s2 = 2*rmsw^2
    y = A/(sqrt(2*pi*s2)) .* exp.(-1/(2*s2) * (x .- xPeak).^2)
    return y
end

function symExpShape(x, xPeak, rmsw, A)
    #SYMEXP create symmetric decreasing exponential with
    # - maximum at 'xPeak'
    # - given root mean square width
    # - given area 'A'

    a = 1/(rmsw*sqrt(2))
    b = a*A/2
    y = b*exp.(-a*abs.(x .- xPeak))
    return y

end

function sech2Shape(x, xPeak, rmsw, A)
    #SECH2SHAPE create sech^2 shape with
    # - maximum at 'xPeak'
    # - given root mean square width
    # - given area 'A'

    lambda = 0.5679 / rmsw;
    a = lambda * A / 2;
    y = a*(sech.(lambda*(x .- xPeak))).^2;

end


function peakWithAbsorption(singlePeak, x, xPeak, xAbs, rmsw, rmswAbs, A)
    # create a peak with an absorption line
    #  - maximum at 'xPeak'
    # - given root mean square width
    # - given area 'A'

    Aabs = A / rdVal(5, 10)
    y1 = singlePeak(x, xPeak, rmsw, A)
    y1xabs = singlePeak(xAbs, xPeak, rmsw, A)
    # Das hier stellt sicher, dass der Puls im Absorptionsbereich nicht ganz auf 0 geht.
    # Das bringt sonst Probleme mit der an der Stelle unbestimmten Phase
    yabs = gaussShape(x, xAbs, rmswAbs, Aabs)
    yabs .*= rdVal(0.7, 0.95) * y1xabs / maximum(yabs)


    y = max.(y1 .- yabs, 0.0) # 1*maximum(y1))  # evtl. negative Werte auf 0 setzen
    Atmp = Trapz.trapz(x,y)
    return y * A/Atmp

end

function doublePeak(singlePeak, x, peak1, peak2, rmsw1, rmsw2, a, A)
    #doublePeak: create double peak from a single peak function
    # - peakDiff: differnce of peak functions
    # - maximum at 'xPeak'
    # - given root mean square width
    # - given area 'A'

    y1 = singlePeak(x, peak1, rmsw1, a*A)
    y2 = singlePeak(x, peak2, rmsw2, (1-a)*A)
    y = y1 + y2
    Atmp = Trapz.trapz(x,y)
    return y * A/Atmp

end

function triplePeak(singlePeak, x, peak1, peak2, peak3, rmsw1, rmsw2, rmsw3, a1, a2, A)
    #doublePeak: create triple peak from a single peak function
    # - peakDiff: differnce of peak functions
    # - maximum at 'xPeak'
    # - given root mean square width
    # - given area 'A'


    y1 = singlePeak(x, peak1, rmsw1, a1*A)
    y2 = singlePeak(x, peak2, rmsw2, a2*A)
    y3 = singlePeak(x, peak3, rmsw3, (1-a1-a2)*A)
    y = y1 + y2 + y3
    Atmp = Trapz.trapz(x,y)
    y = y * A/Atmp
end

function createRandomSignal(params)
    #createRandomSignal: erzeuge Zufallssignal mit gleicher Laenge wie ff
    # mit den Parametern params

    i = rand(1:3)
    # i = 1
    if i == 1
        shapeFunc = gaussShape
    elseif i == 2
        shapeFunc = sech2Shape
    elseif i == 3
        shapeFunc = symExpShape
    end

    i = rand(1:4)
    # i = 4
    # peakDiff = rdVal(0.2, 0.3)
    if i == 1    # absorption
        wAbs = params.wPeak + rdVal(-2, 2) * params.rmsw
        rmswAbs = params.rmsw / rdVal(10, 25)
        abssig = ww -> peakWithAbsorption(shapeFunc, ww, params.wPeak, wAbs, params.rmsw, rmswAbs, params.A)
    elseif i == 2   # double peak
        peak1 = params.wPeak + rdVal(-1.0, -0.5) * params.rmsw
        peak2 = params.wPeak + rdVal(0.5, 1.0) * params.rmsw
        rmsw1 = rdVal(0.2, 0.6) * params.rmsw
        rmsw2 = rdVal(0.2, 0.6) * params.rmsw
        a = rdVal(0.4, 0.6)
        abssig = ww -> doublePeak(shapeFunc, ww, peak1, peak2, rmsw1, rmsw2, a, params.A)
    elseif i == 3  # triple peak
        peak2 = params.wPeak + rdVal(-2.0, -1.0) * params.rmsw
        peak3 = params.wPeak + rdVal(1.0, 2.0) * params.rmsw
        rmsw1 = rdVal(0.4, 0.6) * params.rmsw
        rmsw2 = rdVal(0.4, 0.6) * params.rmsw
        rmsw3 = rdVal(0.4, 0.6) * params.rmsw
        a1 = rdVal(0.2, 0.4)
        a2 = rdVal(0.2, 0.4)
        abssig = ww -> triplePeak(shapeFunc, ww, params.wPeak, peak2, peak3,
                                  rmsw1, rmsw2, rmsw3, a1, a2, params.A)
    elseif i == 4   # single peak
        abssig = ww -> shapeFunc(ww, params.wPeak, params.rmsw, params.A)
    end
    # cubicPhaseFunc(ff, params.fPeak, params.timeShift, params.gdd, params.tod);

    # Rückgabe: Betrag, Phase, Group-Delay
    return (abssig, cubicPhaseFuncAngFreq(params.wPeak, params.timeShift, params.gdd, params.tod)...)
end

# end




"""
  createSHGAmplitude: berechnet die SHG-Frog-Matrix aus dem analytischen Signal im Zeitbereich.
  Dabei wird angenommen, dass das Zeitsignal im Bereich [-T0/2, T0/2 - Ts] definiert ist.
  Wegen der Frequenzverdopplung wird das Produktsignal um 2*wCenter im Frequenzbereich verschoben.
  # Arguments
  `yta::Vector{Complex{T}}`: analytisches Signal im Zeitbereich
  `Ts::T`: Abtastzeit = Delay
  `wCenter::T`: Center-Frequenz
  # Return value: Eine quadratische, komplexe NxN-Matrix, wobei N die Laenge des Zeitvektors yta ist.
"""
function createSHGAmplitude(yta::Vector{Complex{T}}, Ts::T,
                            wCenter::T)::Matrix{Complex{T}} where{T<:Real}
    N = length(yta)
    @assert iseven(N)
    delayIdxVec = (-N ÷ 2):(N ÷ 2 - 1)
    shiftFactor = exp.(- 1im * 2*wCenter * Ts * delayIdxVec)

    shgMat = zeros(Complex{T}, N, N)
    for (matIdx, delayIdx) in enumerate(delayIdxVec)
        ytaShifted = circshift(yta, delayIdx)
        shgMat[matIdx, :] = Ts * fftshift(fft(fftshift(yta .* ytaShifted .* shiftFactor)))
    end
    return shgMat
end


export rdVal, rdSign, createRandomSignal, createSHGAmplitude
