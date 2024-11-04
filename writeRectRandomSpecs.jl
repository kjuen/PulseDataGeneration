# const hst = gethostname()

# # pdgDir = ""
# if hst == "Manni"
#     pdgDir = "/home/kjuen/dvlp/kai/Julia/PulseDataGeneration/pdg"
# elseif hst == "hpc-fat"
#     pdgDir = "/home/klaus/kai/julia/pdg"
# end
# cd(pdgDir)
#

using Pkg
Pkg.activate(".")
Pkg.develop(path="./pdg")
# Pkg.activate("./pdg")
using pdg
using FFTW
using DSP: unwrap
using Statistics
using Interpolations
using Logging

logger = SimpleLogger(stdout, Logging.Info)
disable_logging(Logging.Info)

using CairoMakie
CairoMakie.activate!(visible=false)
# using GLMakie
# GLMakie.activate!()




"""
  generate pulse that is closer to the experimental frog traces than generatePulseSquareSpec
  # Example data (1951):
  # 613: Number of delays (rows)
  # 1064: Number of points in spectrum (columns)
  # 3.2678: Distance between delays [fs]
  # 0.0207: Resolution of spectrum [nm]
  # 516.005: Center Wavelenght [nm]
"""
function generatePulseRectSpec(dirPath::AbstractString,
                               si::SpectrogramInfo;
                               simple::Bool=false,
                               skipWrite::Bool=false)::Tuple


    # TODO: Alles auf single-Datentyp umstellen
    # Ansatz: krasses zero-Padding, so dass am Ende aber w0 und Ts <= die Werte aus dem si-Objekt sind.
    # Dann vermeidet man upsampling beim Zurechtschneiden auf die Zielgröße

    wLaser = 1/2 * c2p / si.CenterWavelen  # hier lebt der Laserpuls
    @info "Laser" wLaser / angFregTHz
    Ts = si.DelayDist
    (ttTarget, llTarget) = generateAxes(si)
    wTargetRange = (c2p/llTarget[end], c2p/llTarget[1])

    w0Prelim = (wTargetRange[2] - wTargetRange[1]) / (si.NumSpecPoints - 1)

    T0Prelim = 2*pi/w0Prelim
    # N = nextpow(2, T0Prelim / si.DelayDist) # gibt eine Float32-SpecMat von etwa 128MB
    # N = nextpow(2, T0Prelim / si.DelayDist) ÷ 2
    N = round(Int, T0Prelim / Ts)
    @info "N" N w0Prelim/angFregTHz T0Prelim/femto
    # TODO Ausprobieren: wie ändert sich die Laufzeit, wenn man
    # - N/2 nimmt
    # - Nicht die nextpow von N nimmt

    idxVec = ((-N ÷ 2):(N ÷ 2 - 1))
    tt = idxVec * Ts
    ws = 2*pi / Ts
    w0 = ws / N
    wwBaseBand = idxVec * w0     # hier lebt der Laserpuls
    wwLaser = wwBaseBand .+ wLaser
    wwFrogTrace = wwBaseBand .+ 2*wLaser   # hier lebt das Spektrogramm

    @info "Zeitachse" Ts extrema(ttTarget)./femto extrema(tt)./femto
    # @info "Target Wellenlängenachse" si.ResSpec/nano extrema(llTarget)./nano
    @info "Große Wellenlängenachse" wTargetRange ./ angFregTHz extrema(wwFrogTrace)./angFregTHz (wwFrogTrace[2] - wwFrogTrace[1]) / angFregTHz


    # Erzeugung eines Pulses auf der dichten Achse
    gdd = 200000  / wLaser^2
    tod = 300 * gdd/wLaser

    if simple   # das hier macht keinen eigenen scope!
        # fixed parameters
        paramsPulse = (wPeak = wLaser*rdVal(0.999, 1.001),
                       A = 1e14,
                       rmsw = 0.6*angFregTHz* rdVal(0.8, 1.2),
                       timeShift = 0,
                       gdd = gdd/50* rdVal(0.0, 1.0),   # gdd >= 0.0 wegen Eindeutigkeit
                       tod = 0)

        # Das hier muesste Fourier-limitiert sein:
        # paramsPulse = (wPeak = wLaser,
        #                A = 1e14,
        #                rmsw = 0.6*angFregTHz,
        #                timeShift = 0,
        #                gdd = 0,   # gdd >= 0.0 wegen Eindeutigkeit
        #                tod = 0)

        Ew = w -> gaussShape(w, paramsPulse.wPeak, paramsPulse.rmsw, paramsPulse.A)
        (Phw, Gdw) = cubicPhaseFuncAngFreq(paramsPulse.wPeak, paramsPulse.timeShift, paramsPulse.gdd, paramsPulse.tod)
    else
        paramsPulse = (wPeak = wLaser*rdVal(0.999, 1.001),
                       A = 1e14,
                       rmsw = 0.75*angFregTHz* rdVal(0.8, 1.2),
                       timeShift = 0,
                       gdd = gdd * rdVal(0.0, 1.0),      # gdd >= 0.0 wegen Eindeutigkeit
                       tod = tod * rdVal(-1.0, 1.0))

        (Ew, Phw, Gdw) = createRandomSignal(paramsPulse)
    end

    Phl(l) = Phw(c2p ./ l)              # Eq. (2.13) im Trebino-Buch
    El(l) = Ew(c2p ./ l) * sc2p ./ l     # Eq. (2.17) im Trebino-Buch

    #FIXME: wenn die Signale eh alle auf max=1 normiert werden, kann man sich das mit der Flaeche im Signalgenerator sparen...
    Ywp = Ew(wwLaser) .* exp.(1im * Phw(wwLaser))
    Ywp ./= maximum(abs.(Ywp))
    # FIXME: Das Zeitsignal scheint symmetrisch zu sein! Das liegt an der nahezu konstanten Phase im Frequenzbereich

    fwhmW = berechneFWHM(abs.(Ywp).^2, wwLaser)
    if fwhmW < 0
        @warn "fwhmW < 0"
        return (3, nothing)
    end
    rmsW = berechneRMSBreite(abs.(Ywp).^2, wwLaser)
    if rmsW < 0
        @warn "rmsW < 0"
        return (3, nothing)
    end
    @info "Laserpuls-Freqbereich" fwhmW/angFregTHz rmsW/angFregTHz

    # Da die Signale auf max=1 normiert werden, braucht man kein 1/delay
    ytaBaseBand = fftshift(ifft(fftshift(Ywp))) # .* exp.(1im * wCenter * delayAxis)
    ytaBaseBand ./= maximum(abs.(ytaBaseBand))
    ytaBaseBand .+= rdVal(0, 1/50) *maximum(abs.(ytaBaseBand))

    # Das hier ist tricky: unwrap(angle.(yta)) liefert ein falsches Ergebnis, denn im Band um
    # wCenter ist ja das Abtasttheorem nicht erfuellt. Daher wird die Phase im Basisband berechnet
    # und einfach wCenter*t addiert.
    ytaPhase = unwrap(angle.(ytaBaseBand)) .+ wLaser * tt
    yta = ytaBaseBand .* exp.(1im * wLaser * tt)

    fwhmT = berechneFWHM(abs.(yta).^2, tt)
    if fwhmT < 0
        @warn "fwhmT < 0"
        return (3, nothing)
    end
    rmsT = berechneRMSBreite(abs.(yta).^2, tt)
    if rmsT < 0
        @warn "rmsT < 0"
        return (3, nothing)
    end
    @info "Laserpuls-Zeitbereich" fwhmT/femto rmsT/femto
    @info "Laser-TBP" fwhmT*fwhmW rmsT*rmsW


    if rmsT < 3*Ts || rmsT*rmsW < 0.2499
        # In diesem Fall kann rmsT nicht mehr gut berechnet werden
        # println("Puls zu schmal")
        @warn "rmsT < 3*Ts || rmsT*rmsW < 0.2499"
        @info rmsT/femto Ts/femto rmsT*rmsW
        return (4, nothing)
    end
    # @show rmsT*rmsW
    # @show fwhmT / femto
    # @show rmsT / femto


    # checke den Abfall des Zeitsignals
    # FIXME: das gleich muesste man doch auch im Frequenzbereich machen
    idxW = idxRangeAboveThres(abs.(Ywp), 1/1000)
    if (idxW[1] == 1) || (idxW[end] == N)
        @warn "Signal fällt nicht auf 0 im Frequenzbereich ab"
        return (5, nothing)
    end
    idxT = idxRangeAboveThres(abs.(yta), 1/1000)
    # if (idxT[1] == 1) || (idxT[end] == N)
    #     @warn "Signal fällt nicht auf 0 im Zeitbereich ab"
    #     return (1, nothing)
    # end
    # @show rmsW / angFregTHz

    # Zeitbereichsdaten für die Datei: Da muss die Zeitachse wird auf ttTarget runtergesampled
    ipctAbs = cubic_spline_interpolation(tt, abs.(yta))
    ipctPhase = cubic_spline_interpolation(tt, ytaPhase)

    # FIXME: Warum sind die Zeitdaten symmetrisch um 0????
    EsData = ComplexSignal(ttTarget/femto, ipctAbs.(ttTarget), ipctPhase.(ttTarget))



    # Uebergang zum Wellenlaengen-Bereich: in der Datei Speck.dat nicht äquidistant!
    # llLaser = range(c2p/wwLaser[end], c2p/wwLaser[1], si.NumSpecPoints)
    llLaser = llTarget .+ si.CenterWavelen
    @info "llLaser" extrema(llLaser)./ nano (llLaser[2] - llLaser[1]) / nano si.ResSpec / nano
    YlpAmpl = El(llLaser)
    YlpAmpl ./= maximum(abs.(YlpAmpl))
    SpecsData = ComplexSignal(llLaser/nano, YlpAmpl, Phl(llLaser))
    fwhmWv = berechneFWHM(reverse(abs.(YlpAmpl).^2), llLaser)
    if fwhmWv < 0
        @warn "fwhmWv < 0"
        return (3, nothing)
    end
    rmsWv = berechneRMSBreite(reverse(abs.(YlpAmpl).^2), llLaser)
    if rmsWv < 0
        @warn "rmsWv < 0"
        return (3, nothing)
    end

    # Spec Mat: lebt in delayAxis x wAxisShifted2
    # wAxisShifted2 = wAxisBaseBand .+ 2*wCenter   # hier lebt das Produkt-Signal
    # wAxisShifted2 = ww .+ wCenter   # hier lebt das Produkt-Signal
    @time shgAmplMat = createSHGAmplitude(yta, Ts, wLaser)

    # Zu der Groesse reduzieren, die fuer das Target-Spektrogramm ausreicht
    (tidx1, tidx2) = idxRangeWithinLimits(tt, extrema(ttTarget))
    tIdxRange = (tidx1-1):tidx2
    ttOut = tt[tIdxRange]
    @info "Zeitauschnitt" extrema(ttOut)./femto  length(ttOut) extrema(ttTarget)./femto length(ttTarget)

    (widx1, widx2) = idxRangeWithinLimits(wwFrogTrace, wTargetRange)
    wIdxRange =  (widx1-1):widx2
    wwOut = wwFrogTrace[wIdxRange]
    @info "Frequenzausschnitt" extrema(wwOut)./angFregTHz length(wwOut) wTargetRange./angFregTHz si.NumSpecPoints

    shgAmplMatOut = shgAmplMat[tIdxRange, wIdxRange]



    # 2D-Interpolation: ist das wirklich noetig??
    # itpc = cubic_spline_interpolation((ww[wIdxRange], tt[tIdxRange]), shgAmplMatOut);
    #shgAmplMatInterp =


    shgIntMatOut = abs.(shgAmplMatOut).^2
    # Uebergang zum Wellenlaengenbereich
    (llOut, shgIntMatWvOut) = intensityMatFreq2Wavelength(wwOut, shgIntMatOut)
    @info size(shgIntMatWvOut)

    # if !specMatCorrectlySampled(shgIntMatWvOut; nPix=12, thres=1e-3)
    #     @warn "SpecMat not an island in a sea of zeros (2)"
    #     return (2, nothing)
    # end

    # Gauss-Noise zufügen
    shgIntMatWvMax = maximum(shgIntMatWvOut)
    # shgIntMatWv_gn01 = copy(shgIntMatWv) + 0.001 * shgIntMatWvMax * randn(size(shgIntMatWv))
    shgIntMatWv_gn10 = copy(shgIntMatWvOut) + 0.01 * shgIntMatWvMax * randn(size(shgIntMatWvOut))
    # shgIntMatWv_gn30 = copy(shgIntMatWv) + 0.03 * shgIntMatWvMax * randn(size(shgIntMatWv))
    # @show traceError(shgIntMatWv, shgIntMatWv_gn01)
    # @show traceError(shgIntMatWv, shgIntMatWv_gn10)
    # @show traceError(shgIntMatWv, shgIntMatWv_gn30)  Komisch: der Wert hier passt nicht


    asData = (ttOut, llOut,
              (;gn00 = shgIntMatWvOut,
               gn10 = shgIntMatWv_gn10))
    # In Verzeichnis schreiben
    if !skipWrite
        pd = (fwhmT=fwhmT, rmsT=rmsT,
              fwhmW = fwhmW, rmsW = rmsW,
              fwhmWv=fwhmWv, rmsWv=rmsWv,
              wCenter=2*wLaser,
              wp=paramsPulse.wPeak)
        if length(dirPath) > 0
            writePulseDataToDir(dirPath, EsData, SpecsData, asData)
            writePulseFeatures(dirPath, pd)
        end
        plotPulseData(dirPath, EsData, SpecsData, asData, pd)
    end


    return (0, rmsT * rmsW)
end




#* Loop zum Erzeugen der Daten
# p = joinpath(pdgDir, "simulatedTraces_3", "abc"x1)
# retCode = generatePulse(p)
# @show retCode)./nano

targetDir = "/home/kjuen/dvlp/kai/data/simulatedTraces_6"
si = SpectrogramInfo(613, 1064, 3.2678*femto, 0.0207*nano, 516.005*nano)

count = 0
success = 0
simpleCount = 0
errCodeCounter = zeros(Int, 5)
tbpVec = []
nFiles = 10
@time begin
    while success < nFiles
        global count += 1
        if mod(count, 5) == 0
            @show count
        end
        local p = joinpath(targetDir, "s"*string(success+1))
        simple= rand() < 1/300
        local (retCode, tbp) = generatePulseRectSpec(p, si; simple=simple,
                                                     skipWrite=false)
        if retCode == 0
            global success += 1
            if simple
                global simpleCount += 1
            end
            println(success)
            push!(tbpVec, tbp)
        else
            global errCodeCounter[retCode] += 1
        end

    end
end
@assert success == nFiles
# println("Erolgsrate: $(round(100 * nFiles/count, digits=1))%")
# @show count
# @show simpleCount
# @show keinAbfallZeitbereich
# @show keinAbfallFreqbereich
# @show noIsland
# @show noFWHM
# @show rmstTooSmall


# f = Figure()
# hist(f[1,1], tbpVec, bins=50)
# # save("tbp.pdf", f)
