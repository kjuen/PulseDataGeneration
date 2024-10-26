

# const hst = gethostname()

# # pdgDir = ""
# if hst == "Manni"
#     pdgDir = "/home/kjuen/dvlp/kai/Julia/PulseDataGeneration/pdg"
# elseif hst == "hpc-fat"
#     pdgDir = "/home/klaus/kai/julia/pdg"
# end
# cd(pdgDir)

# # using GLMakie
# using DSP.Unwrap: unwrap
# using DSP: conv
# using DSP.Util: hilbert
# using FFTW
# using Interpolations
# using Statistics
# using CSV
# using CairoMakie

using Pkg
Pkg.activate(".")
Pkg.develop(path="./pdg")
# Pkg.activate("./pdg")
using pdg
using CairoMakie
using FFTW
using DSP: unwrap
using Statistics

CairoMakie.activate!(visible=false)
# using GLMakie
# GLMakie.activate!()

# include(joinpath(pdgDir, "PDG.jl"))

function relDiff(s1, s2)
    N = length(s1)
    @assert length(s2) == N
    Z = sum(abs.(s1).^2)
    return sum(abs.(s1 - s2).^2) / Z
end


# include("pdgFunctionsAndConstants.jl")
# include("pdgPulseCreationFreqDomain.jl")
# include("pdgPlots.jl")
# include("pdgFileIO.jl")

# Änderungen am Pulsgenerator
# - RMS-Breiten werden jetzt auch mit den Intensitäten berechnet (2.9.24)



"""
 write all available data to a directory. This consists of
 Es.dat: signal in time domain (five column format)
 Specs.dat: signal in wavelength domain (five column format)
 as.dat: shg matrix
 # Arguments:
 - `fullPath`: full path of directory to write to
 - `EsData`: ComplexSignal object (see pdgFileIO.jl) for file Es.dat.
 - `SpecsData`: ComplexSignal object (see pdgFileIO.jl) for file Specs.dat.
 - `asData`: tuple of (time vector, wavelength vector, shg-Mat) of consistent size
 - `delim="  "`: delimiter

"""
function writePulseDataToDir(fullPath::String,
                             EsData::ComplexSignal{T},
                             SpecsData::ComplexSignal{T},
                             asData::Tuple{AbstractVector{T},AbstractVector{T}, NamedTuple},
                             delim="  ") where{T<:Real}
    if isfile(fullPath)
        error("$(fullPath) is an existing file")
    end
    if !isdir(fullPath)
        mkpath(fullPath)
    end
    @assert length(asData[1]) == size(asData[3][1], 1)
    @assert length(asData[2]) == size(asData[3][1], 2)


    writeFiveColumns(joinpath(fullPath, "Es.dat"),
                     EsData;
                     delim=delim)

    writeFiveColumns(joinpath(fullPath, "Specs.dat"),
                     SpecsData;
                     delim=delim)


    for matName in keys(asData[3])
        fn = "as_" * string(matName) * ".dat"
        # Skaliert auf max=2^16
        writeFrogTraceToFile(joinpath(fullPath, fn),
                             asData[1], asData[2], asData[3][matName])
    end

end


"""
write some pulse features (width and peak frequency) to a text file
"""
function writePulseFeatures(fullPath::String, pulseData=Tuple)
    if length(fullPath) != 0
        @assert isdir(fullPath)

        open(joinpath(fullPath, "SimulatedPulseData.txt"), "w") do io
            for p in pairs(pulseData)
                writedlm(io, [p.first p.second], ";")
            end
        end
    end
end


"""
    Plot pulse data and spectrogram to pdf file
"""
function plotPulseData(dirPath::String,
                       EsData::ComplexSignal{T},
                       SpecsData::ComplexSignal{T},
                       asData::Tuple{AbstractVector{T},AbstractVector{T},NamedTuple},
                       pulseData=Tuple) where {T<:Real}
    colors = Makie.wong_colors()
    commonAxAttr = (titlesize=12, xlabelsize=10, ylabelsize=10,
                    xticklabelsize=8, yticklabelsize=8, ylabelvisible=false, xlabelvisible=false,
                    xminorgridvisible=true, yminorgridvisible=true,
                    topspinevisible=false, leftspinecolor = colors[1], rightspinecolor = colors[2])
    commonLeftAttr = (ylabelcolor = colors[1], yticklabelcolor = colors[1], ytickcolor = colors[1])
    commonRightAttr = (ylabelcolor = colors[2], yticklabelcolor = colors[2], ytickcolor = colors[2],
                       yaxisposition = :right)
    legendAttr = (framevisible = false, labelsize=8)

    delayAxis = EsData.xVals * femto
    Es = getSignal(EsData)
    ll = SpecsData.xVals * nano
    Ywl = getSignal(SpecsData)

    newWvAxis = asData[2]
    shgIntMatWv = asData[3][1]

    # Plotten
    f1= Figure()
    # tPlotLims = (-3, 3) .* pulseData.fwhmT
    tPlotLims = extrema(delayAxis)

    # Time Domain
    tempIntensity = abs.(Es).^2
    idxPhase = idxRangeAboveThres(tempIntensity, 1/100)
    omt = mydiff5(EsData.phase[idxPhase], delayAxis[2] - delayAxis[1])# .+ pulseData.wp
    instFreqLims =  (0.90, 1.1) .* pulseData.wCenter

    axl = Axis(f1[1,1]; commonAxAttr..., commonLeftAttr...,
               ylabel="Intensity (a.u.)",
               xlabel="Time in fs",
               title="Pulse on time axis (in fs): FWHM=$(round(pulseData.fwhmT/femto, digits=2))fs",
               limits=(tPlotLims ./ femto, nothing))
    scatterlines!(axl, delayAxis/femto, tempIntensity, label="Intensity",
                  color = colors[1], linestyle=:solid, linewidth=1, marker=:x, markersize=4)
    axislegend(axl, position=:lt; legendAttr...)

    axr = Axis(f1[1,1]; commonAxAttr..., commonRightAttr...,
               ylabel="Inst. freq. in THz",
               limits=(tPlotLims ./ femto,  instFreqLims ./ angFregTHz))
    lines!(axr, delayAxis[idxPhase[3:end-2]]/femto,
           omt / angFregTHz, label="Inst. Freq (THz)",
           color = colors[2], linestyle=:solid, linewidth=1)

    # scatterlines!(axr, delayAxis[idxPhase]/femto,
    #        EsData.phase[idxPhase], label="Phase",
    #        color = colors[2], linestyle=:solid, linewidth=1)


    axislegend(axr, position=:rt; legendAttr...)

    hidespines!(axr)
    hidexdecorations!(axr)
    hideydecorations!(axr, ticks=false, ticklabels = false, label=false)

    # Wavelength Domain
    #llLims = (((-3, 3) .* pulseData.fwhmWv .+ c2p /pulseData.wCenter)./nano, nothing)
    llLims = (extrema(ll)./nano, nothing)
    wvInt = abs.(Ywl).^2
    wvInt ./= maximum(wvInt)
    idxPhase = idxRangeAboveThres(wvInt, 1/100)
    phl = unwrap(angle.(Ywl))
    axl = Axis(f1[2,1]; commonAxAttr..., commonLeftAttr...,
               title="Pulse on wavelength axis (in nm): FWHM=$(round(pulseData.fwhmWv/nano, digits=2))nm",
               limits=llLims)
    scatterlines!(axl, ll/nano, wvInt , label="Spectral Intensity",
                  color = colors[1], linestyle=:solid, linewidth=1, marker=:x, markersize=4)
    axislegend(axl, position=:lt; legendAttr...)

    axr = Axis(f1[2,1]; commonAxAttr..., commonRightAttr...,
               ylabel="Phase",
               limits=llLims)
    lines!(axr, ll[idxPhase]/nano, phl[idxPhase], label="Phase in rad",
           color = colors[2], linestyle=:solid, linewidth=1)
    axislegend(axr, position=:rt; legendAttr...)

    hidespines!(axr)
    hidexdecorations!(axr)
    hideydecorations!(axr, ticks=false, ticklabels = false, label=false)

    axSHG1 = plotWavelengthSpectrogram(f1[1,2], delayAxis, newWvAxis,
                                       shgIntMatWv)
    tbp = round(pulseData.rmsW * pulseData.rmsT, digits=2)
    axSHG1.title="Wavelength-SHG: TBP=$(tbp)"


    # MAN BEACHTE: durch das Interpolieren kann die Matrix shgIntMatWv vom Betrag
    # her sehr kleine, aber negative Werte haben. Dann erzeugt die log-Skala einen Fehler.
    # Daher wird hier noch mal der Betrag der Matrix genommen.
    axSHG2 = plotWavelengthSpectrogram(f1[2,2], delayAxis, newWvAxis,
                                       abs.(shgIntMatWv);
                                       logScale=true);
    axSHG2.title="Wavelength-SHG: Log Scale"


    if length(dirPath) > 0
        @assert isdir(dirPath)
        save(joinpath(dirPath, "SimulatedData.pdf"), f1)
    # else
    #     display(f1)
    end

    # Plot der 4 Noise-Levels
    # f2= Figure()
    # matNames = keys(asData[3])
    # axSHG11 = plotWavelengthSpectrogram(f2[1,1], delayAxis, newWvAxis,
    #                                     asData[3][matNames[1]])
    # axSHG11.title=string(matNames[1])
    # axSHG12 = plotWavelengthSpectrogram(f2[1,2], delayAxis, newWvAxis,
    #                                     asData[3][matNames[2]])
    # axSHG12.title=string(matNames[2])
    # axSHG21 = plotWavelengthSpectrogram(f2[2,1], delayAxis, newWvAxis,
    #                                     asData[3][matNames[3]])
    # axSHG21.title=string(matNames[3])
    # axSHG22 = plotWavelengthSpectrogram(f2[2,2], delayAxis, newWvAxis,
    #                                    asData[3][matNames[4]])
    # axSHG22.title=string(matNames[4])


    # if length(dirPath) > 0
    #     @assert isdir(dirPath)
    #     save(joinpath(dirPath, "SimulatedData_noisyTraces.pdf"), f2)
    # end

end


#* Random pulse generation


# nGrid: Groesse der Frog-Trace
# Strategie
# Erst mit vielen Sample-Points samplen, dann F.T. in den Zeitbereich
# Damit wird dann das Spektrogramm berechnet
# Dann wird im Frequenzbereich runtergesampelt,
# und damit noch mal in den Zeitbereich transformiert, um den runtergesampelten Zeitbereich zu bekommen
function generatePulse(dirPath::AbstractString;
                       nGrid::Int=256,
                       delay::Float64 = 1.5*femto,
                       wCenter::Float64 = 500*angFregTHz,
                       simple::Bool=false,
                       skipWrite::Bool=false)::Tuple

    #* Erzeugung eines Pulses auf der dichten Achse
    gdd = 4000  / wCenter^2
    tod = 75 * gdd/wCenter

    if simple   # das hier macht keinen eigenen scope!
        # fixed parameters
        paramsPulse = (wPeak = wCenter*rdVal(0.95, 1.05),
                       A = 1e14,
                       rmsw = 50*angFregTHz* rdVal(0.1, 1.5),
                       timeShift = 0,
                       gdd = gdd/5* rdVal(0.0, 1.0),   # gdd >= 0.0 wegen Eindeutigkeit
                       tod = 0)

        # Das hier muesste Fourier-limitiert sein:
        # paramsPulse = (wPeak = wCenter,
        #                A = 1e14,
        #                rmsw = 10*angFregTHz,
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
                       gdd = gdd * rdVal(0.0, 1.0),      # gdd >= 0.0 wegen Eindeutigkeit
                       tod = tod * rdVal(-1.0, 1.0))

        (Ew, Phw, Gdw) = createRandomSignal(paramsPulse)
    end

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
targetDir = "/home/kjuen/dvlp/kai/data/simulatedTraces_4"
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
        local (retCode, tbp) = generatePulse(p; delay= 1.5*femto, simple=false,
                                             nGrid=256, skipWrite=true)
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
