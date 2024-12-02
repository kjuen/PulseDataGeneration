

# const hst = gethostname()

# # pdgDir = ""
# if hst == "Manni"
#     pdgDir = "/home/kjuen/dvlp/kai/Julia/PulseDataGeneration/pdg"
# elseif hst == "hpc-fat"
#     pdgDir = "/home/klaus/kai/julia/pdg"
# end
# cd(pdgDir)


using Pkg
Pkg.activate("..")
Pkg.develop(path="../pdg")
# Pkg.activate("./pdg")
using pdg
using FFTW
using DSP: unwrap, conv
using Statistics

using CairoMakie
CairoMakie.activate!(visible=false)
#using GLMakie
#GLMakie.activate!()


function relDiff(s1, s2)
    N = length(s1)
    @assert length(s2) == N
    Z = sum(abs.(s1).^2)
    return sum(abs.(s1 - s2).^2) / Z
end




function plotRana(It, delayAxis, Ew, wAxisShifted,
                  shgIntMatWv, newWvAxis,
                  wAxisShifted2, shgIntMat,
                  pulseData, dirPath)
        # Plot der Marginals
    colors = Makie.wong_colors()

    meanDelayProfile = map(mean, eachrow(shgIntMatWv))
    meanDelayProfile ./= maximum(meanDelayProfile)
    @assert length(meanDelayProfile) == length(delayAxis)
    fwhmT2 = berechneFWHM(meanDelayProfile, delayAxis) / sqrt(2)

    f2= Figure()
    tPlotLims = extrema(delayAxis)
    ax = Axis(f2[1,1];
              ylabel="Intensity (a.u.)",
              xlabel="Time in fs",
              title="Delay marginal: FWHM=$(round(fwhmT2/femto, digits=2))fs",
              limits=(tPlotLims ./ femto, nothing))
    lines!(ax, delayAxis/femto, meanDelayProfile, label="Delay marginal",
           color = colors[1], linestyle=:solid, linewidth=2)
    Ts = delayAxis[2] - delayAxis[1]
    (iac, tiac) = autocorr(It; Ts)
    iac ./= maximum(iac)
    lines!(ax, tiac/femto, iac, label="Intensity-AC", color = colors[2],
                  linestyle=:dash, linewidth=2)
    axislegend(ax, position=:rt; framevisible = false, labelsize=8)

    meanSpecProfile = map(mean, eachcol(shgIntMat))
    meanSpecProfile ./= maximum(meanSpecProfile)
    @assert length(wAxisShifted2) == length(meanSpecProfile)
    fwhmW2 = berechneFWHM(meanSpecProfile, wAxisShifted2)

    Sw = abs.(Ew).^2
    ScSw = conv(Sw, Sw)
    ScSw ./= maximum(ScSw)
    # Rumspielen mit Rauschen
    # ScSw .*= (0.02 * randn(size(ScSw)) .+ 1.0)
    # Frequenzachse fuer Faltung:
    Nc = length(ScSw)
    w0 = wAxisShifted[2] - wAxisShifted[1]
    wwConv = (0:(Nc-1))*w0 .+ 2*wAxisShifted[1]

    @show wAxisShifted[1] / angFregTHz

    ax2 = Axis(f2[2,1];
               ylabel="Intensity (a.u.)",
               xlabel="Frequency in THz",
               title="Frequency marginal: FWHM=$(round(fwhmW2/angFregTHz, digits=2))THz",
               limits=(extrema(wAxisShifted2) ./ angFregTHz, nothing))
    lines!(ax2, wAxisShifted2/angFregTHz, meanSpecProfile, label="Freq marginal",
           color = colors[1], linestyle=:solid, linewidth=2)
    lines!(ax2, wwConv/angFregTHz, ScSw, label="S*S(w)",
           color = colors[2], linestyle=:dash, linewidth=2)
    axislegend(ax2, position=:lt; framevisible = false, labelsize=8)


    (SwRana, wAxisRana) = rana(ScSw, w0; alpha=1.0, beta=10.0, gamma = 10.0)
    ax3 = Axis(f2[1,2];
               xlabel="Frequenz in THz",
               ylabel="S(w)",
               title="Frequenzbereich",
               limits=(extrema(wAxisShifted/angFregTHz), nothing))
    lines!(ax3, wAxisShifted/angFregTHz, Sw, label="Sw",
           color = colors[1], linestyle=:solid, linewidth=2)
    lines!(ax3, (wAxisRana .+ pulseData.wCenter)/angFregTHz, abs.(SwRana), color = colors[2], linestyle=:dash, linewidth=2)
    # lines!(ax3, (wAxisRana)/angFregTHz, abs.(fftshift(SwRana)), color = colors[4], linestyle=:dash, linewidth=2)


    wAxisRana2 = (0:(Nc-1))*w0 #  .+ wAxisShifted[1]
    @assert length(SwRana) == Nc
    ax4 = Axis(f2[2,2];
               xlabel="Frequency in THz",
               ylabel="S*S",
               title="Frequenzbereich",
               limits=(extrema(wAxisRana) ./ angFregTHz, nothing))
    lines!(ax4, (wAxisRana .+ wAxisShifted[1])/angFregTHz, abs.(SwRana), color = colors[1], linestyle=:solid, linewidth=2)
    # lines!(ax4, wwConv/angFregTHz, ScSw, label="S*S(w)", color = colors[1], linestyle=:solid, linewidth=2)


    save(joinpath(dirPath, "SimulatedMarginals.pdf"), f2)
end


#* Random pulse generation

# Erzeuge ein quadratisches Spektrogramm der Groesse nGrid x nGrid
# nGrid: Groesse der Frog-Trace
# Strategie
# Erst mit vielen Sample-Points samplen, dann F.T. in den Zeitbereich
# Damit wird dann das Spektrogramm berechnet
# Dann wird im Frequenzbereich runtergesampelt,
# und damit noch mal in den Zeitbereich transformiert, um den runtergesampelten Zeitbereich zu bekommen
function generatePulseSquareSpec(dirPath::AbstractString;
                                 nGrid::Int=256,
                                 delay::Float64 = 1.5*femto,
                                 wCenter::Float64 = 500*angFregTHz,
                                 simple::Bool=false,
                                 skipWrite::Bool=false)::Tuple

    #* Erzeugung eines Pulses auf der dichten Achse
    gdd = 800  / wCenter^2
    tod = 25 * gdd/wCenter

    if simple   # das hier macht keinen eigenen scope!
        # fixed parameters
        paramsPulse = (wPeak = wCenter*rdVal(0.95, 1.05),
                       A = 1e14,
                       rmsw = 50*angFregTHz* rdVal(0.1, 1.5),
                       timeShift = 0,
                       gdd = gdd* rdVal(0.0, 1.0),   # gdd >= 0.0 wegen Eindeutigkeit
                       tod = 0)

        # Das hier muesste Fourier-limitiert sein:
        # paramsPulse = (wPeak = wCenter,
        #                A = 1e14,
        #                rmsw = 50*angFregTHz,
        #                timeShift = 0,
        #                gdd = 0,   # gdd >= 0.0 wegen Eindeutigkeit
        #                tod = 0)

        Ew = w -> gaussShape(w, paramsPulse.wPeak, paramsPulse.rmsw, paramsPulse.A)
        (Phw, Gdw) = cubicPhaseFuncAngFreq(paramsPulse.wPeak, paramsPulse.timeShift, paramsPulse.gdd, paramsPulse.tod)
    else
        paramsPulse = (wPeak = wCenter*rdVal(0.95, 1.05),
                       A = 1e14*rdVal(0.9, 1.1),
                       rmsw = 50.0*angFregTHz* rdVal(0.5, 1.8),
                       timeShift = 0,
                       gdd = gdd/3 * rdVal(0.0, 1.0),      # gdd >= 0.0 wegen Eindeutigkeit
                       tod = tod * rdVal(-1.0, 1.0))

        (Ew, Phw, Gdw) = createRandomSignal(paramsPulse)
    end

    # Theta-foermige Phase
    # Phw(om) = (om .> paramsPulse.wPeak) * pi
    # Phw = om -> zeros(Float64, size(om))


    Phl(l) = Phw(c2p ./ l)              # Eq. (2.13) im Trebino-Buch
    El(l) = Ew(c2p ./ l) * sc2p ./ l     # Eq. (2.17) im Trebino-Buch


    # Konstruktion von Zeit- und Frequenzachse
    idxVec = ((-nGrid ÷ 2):(nGrid ÷ 2 - 1))
    delayAxis = idxVec * delay
    ws = 2*pi/delay
    w0 = ws/nGrid
    wAxisBaseBand = idxVec * w0
    wAxisShifted = wAxisBaseBand .+ wCenter    # hier lebt der Laser-Puls


    #FIXME: wenn die Signale eh alle auf max=1 normiert werden, kann man sich das mit der Flaeche im Signalgenerator sparen...
    Ywp = Ew(wAxisShifted) .* exp.(1im * Phw(wAxisShifted))
    Ywp ./= maximum(abs.(Ywp))
    fwhmW = berechneFWHM(abs.(Ywp).^2, wAxisShifted)
    if fwhmW < 0
        println("no fwhmW")
        return (3, nothing)
    end
    rmsW = berechneRMSBreite(abs.(Ywp).^2, wAxisShifted)
    if rmsW < 0
        println("no rmsW")
        return (3, nothing)
    end

    # Da die Signale auf max=1 normiert werden, braucht man kein 1/delay
    ytaBaseBand = fftshift(ifft(fftshift(Ywp))) # .* exp.(1im * wCenter * delayAxis)
    ytaBaseBand ./= maximum(abs.(ytaBaseBand))

    # Das hier ist tricky: unwrap(angle.(yta)) liefert ein falsches Ergebnis, denn im Band um
    # wCenter ist ja das Abtasttheorem nicht erfuellt. Daher wird die Phase im Basisband berechnet
    # und einfach wCenter*t addiert.
    ytaPhase = unwrap(angle.(ytaBaseBand)) .+ wCenter * delayAxis
    yta = ytaBaseBand .* exp.(1im * wCenter * delayAxis)

    fwhmT = berechneFWHM(abs.(yta).^2, delayAxis)
    if fwhmT < 0
        println("no fwhmT")
        return (3, nothing)
    end
    rmsT = berechneRMSBreite(abs.(yta).^2, delayAxis)
    if rmsT < 0
        println("no rmsT")
        return (3, nothing)
    end


    if rmsT < 3*delay || rmsT*rmsW < 0.25
        # In diesem Fall kann rmsT nicht mehr gut berechnet werden
        return (4, nothing)
    end

    # checke den Abfall des Zeitsignals
    # FIXME: das gleich muesste man doch auch im Frequenzbereich machen
    idxW = idxRangeAboveThres(abs.(Ywp), 1/1000)
    if (idxW[1] == 1) || (idxW[end] == nGrid)
        # println("Signal fällt nicht auf 0 im Frequenzbereich ab")
        return (5, nothing)
    end
    idxT = idxRangeAboveThres(abs.(yta), 1/1000)
    if (idxT[1] == 1) || (idxT[end] == nGrid)
        # println("Signal fällt nicht auf 0 im Zeitbereich ab")
        return (1, nothing)
    end
    # @show rmsW / angFregTHz

    EsData = ComplexSignal(delayAxis/femto, abs.(yta), ytaPhase)

    # Uebergang zum Wellenlaengen-Bereich: in der Datei Speck.dat nicht äquidistant!
    lAxisShifted = c2p ./ wAxisShifted
    YlpAmpl = El(lAxisShifted)
    YlpAmpl ./= maximum(abs.(YlpAmpl))
    SpecsData = ComplexSignal(lAxisShifted/nano, YlpAmpl, Phl(lAxisShifted))
    fwhmWv = berechneFWHM(reverse(abs.(YlpAmpl).^2), reverse(lAxisShifted))
    if fwhmWv < 0
        # println("no fwhmWv")
        return (3, nothing)
    end
    rmsWv = berechneRMSBreite(reverse(abs.(YlpAmpl).^2), reverse(lAxisShifted))
    if rmsWv < 0
        # println("no rmsWv")
        return (3, nothing)
    end


    # Spec Mat: lebt in delayAxis x wAxisShifted2
    wAxisShifted2 = wAxisBaseBand .+ 2*wCenter   # hier lebt das Produkt-Signal
    shgAmplMat = createSHGAmplitude(yta, delay, wCenter)
    shgIntMat = abs.(shgAmplMat).^2


    # Uebergang zum Wellenlaengenbereich
    (lAxisShifted2, shgIntMatWv) = intensityMatFreq2Wavelength(wAxisShifted2, shgIntMat)
    (wAxisShifted2_2, shgIntMat_2) = intensityMatWavelength2Freq(lAxisShifted2, shgIntMatWv)
    @show maximum(abs.(wAxisShifted2_2 .- wAxisShifted2)) ./ maximum(wAxisShifted2)
    @show maximum(abs.(shgIntMat_2 .- shgIntMat)) ./ maximum(shgIntMat)

    if !specMatCorrectlySampled(shgIntMatWv; nPix=12, thres=1e-5)
        # println("SpecMat not an island in a sea of zeros (2)")
        return (2, nothing)
    end

    # Gauss-Noise zufügen
    # shgIntMatWvMax = maximum(shgIntMatWv)
    # shgIntMatWv_gn01 = copy(shgIntMatWv) + 0.001 * shgIntMatWvMax * randn(size(shgIntMatWv))
    # shgIntMatWv_gn10 = copy(shgIntMatWv) + 0.01 * shgIntMatWvMax * randn(size(shgIntMatWv))
    # shgIntMatWv_gn30 = copy(shgIntMatWv) + 0.03 * shgIntMatWvMax * randn(size(shgIntMatWv))
    # @show traceError(shgIntMatWv, shgIntMatWv_gn01)
    # @show traceError(shgIntMatWv, shgIntMatWv_gn10)
    # @show traceError(shgIntMatWv, shgIntMatWv_gn30)  Komisch: der Wert hier passt nicht


    # asData = (delayAxis, lAxisShifted2,
    #           (gn00 = shgIntMatWv,
    #            gn01 = shgIntMatWv_gn01,
    #            gn10 = shgIntMatWv_gn10,
    #            gn30 = shgIntMatWv_gn30))

    asData = (delayAxis, lAxisShifted2,
              (;gn00 = shgIntMatWv))


    println("*********")
    # @show fwhmWv / nano
    # @show fwhmT / femto

    meanDelayProfile = map(mean, eachrow(shgIntMatWv))
    fwhmT2 = berechneFWHM(meanDelayProfile, delayAxis) / sqrt(2)
    @show round(fwhmT2 / fwhmT, digits=3)

    meanSpecProfile = map(mean, eachcol(shgIntMatWv))
    fwhmWv2 = berechneFWHM(meanSpecProfile, lAxisShifted2) * 2 * sqrt(2)
    @show round(fwhmWv2 / fwhmWv, digits=3)
    fwhmTbd = fwhmWv * fwhmT
    fwhmTbd2 = fwhmWv2 * fwhmT2
    @show fwhmTbd2 / fwhmTbd

    # In Verzeichnis schreiben
    if !skipWrite
        pd = (fwhmT=fwhmT, rmsT=rmsT,
              fwhmW = fwhmW, rmsW = rmsW,
              fwhmWv=fwhmWv, rmsWv=rmsWv,
              wCenter=wCenter,
              wp=paramsPulse.wPeak)
        if length(dirPath) > 0
            writePulseDataToDir(dirPath, EsData, SpecsData, asData)
            writePulseFeatures(dirPath, pd)
        end
        plotPulseData(dirPath, EsData, SpecsData, asData, pd)
        plotRana(abs.(yta).^2, delayAxis, Ywp, wAxisShifted, shgIntMatWv, lAxisShifted2,
                 wAxisShifted2_2, shgIntMat_2, pd, dirPath)
    end
    # @show paramsPulse.rmsw  /angFregTHz
    # @show paramsPulse.gdd
    # @show paramsPulse.tod

    # Check, ob man den Imaginärteil per Hilbert-Trafo berechnen kann
    # yta2 = hilbert(real.(yta))
    # @show relDiff(real.(yta), real.(yta2))
    # Warum hier ein -? Irgendwie kann das doch nicht richtig sein!
    # @show relDiff(imag.(yta), -imag.(yta2))

    return (0, rmsT * rmsW)
end








#* Loop zum Erzeugen der Daten
# p = joinpath(pdgDir, "simulatedTraces_3", "abc"x1)
# retCode = generatePulse(p)
# @show retCode)./nano

# targetDir = joinpath(pdgDir, "testWriteRandomSignal")
targetDir = "/home/kjuen/dvlp/kai/data/simulatedTraces_3"
# targetDir = joinpath(pdgDir, "VergleichSimulationFrogEinfachePulse_120924")
# targetDir = joinpath(pdgDir, "VergleichSimulationFrogMehrPulse_120924")
# targetDir = "/home/klaus/kai/data/grid_512_v1"

count = 0
success = 0
simpleCount = 0
keinAbfallZeitbereich = 0
keinAbfallFreqbereich = 0
noIsland = 0
noFWHM = 0
rmstTooSmall = 0
tbpVec = []
nFiles = 10
@time begin
    while success < nFiles
        global count += 1
        if mod(count, 10000) == 0
            @show count
        end
        local p = joinpath(targetDir, "s"*string(success+1))
        simple= rand() < 1/300
        local (retCode, tbp) = generatePulseSquareSpec(p; delay= 1.5*femto, simple=false,
                                                       nGrid=256, skipWrite=false)
        @show retCode
        if retCode == 0
            global success += 1
            if simple
                global simpleCount += 1
            end
            println(success)

            push!(tbpVec, tbp)
        elseif retCode == 1
            global keinAbfallZeitbereich += 1
        elseif retCode == 2
            global noIsland +=1
        elseif retCode == 3
            global noFWHM += 1
        elseif retCode == 4
            global rmstTooSmall += 1
        elseif retCode == 5
            global keinAbfallFreqbereich += 1
        else
            error("Unknown ret-code")
        end

    end
end
@assert success == nFiles
println("Erolgsrate: $(round(100 * nFiles/count, digits=1))%")
@show count
@show simpleCount
@show keinAbfallZeitbereich
@show keinAbfallFreqbereich
@show noIsland
@show noFWHM
@show rmstTooSmall


# f = Figure()
# hist(f[1,1], tbpVec, bins=50)
# # save("tbp.pdf", f)