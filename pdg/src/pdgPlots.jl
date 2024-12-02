# Plots of the generated pulse profiles

# using GLMakie

# include("pdgFunctionsAndConstants.jl")

using GLMakie
using CairoMakie
colors = Makie.wong_colors()


commonAxAttr = (topspinevisible=false, leftspinecolor = colors[1], rightspinecolor = colors[2])
commonLeftAttr = (ylabelcolor = colors[1], yticklabelcolor = colors[1], ytickcolor = colors[1])
commonRightAttr = (ylabelcolor = colors[2], yticklabelcolor = colors[2], ytickcolor = colors[2],
                   yaxisposition = :right)




"""
 Plot of amplitude and group delay in angular frequency domain on THz axis
 # Arguments
 - `f`: the grid position of the plot
 - `ww`: angular frequency vector
 - `Ywp`: complex signal, positive frequence part
 - `specFreqRange`: frequency range which is displayed in plot
 # Optional Arguments
 - `thres (default=0.01)`: Only where amplitude exceeds this fraction of maximum, phasi ingo is displayed
"""
function plotFreqDomain(f, ww, Ywp, specFreqRange; thres=0.01)

    @assert length(ww) == length(Ywp)
    w0 = ww[2] - ww[1]
    Sw = abs2.(Ywp)
    axl = Axis(f; commonAxAttr..., commonLeftAttr...,
               ylabel="Intensity (a.u.)",
               xlabel="Angular Frequency in THz",
               title="Pulse on frequency axis",
               limits=(specFreqRange./angFregTHz, nothing))
    lines!(axl, ww/angFregTHz, Sw./maximum(Sw), label="Intensity",
           color = colors[1], linestyle=:solid, linewidth=2)
    # axislegend(axl, position=:lt, framevisible = true)


    axr = Axis(f; commonAxAttr..., commonRightAttr...,
               ylabel="Group delay in fs",
               limits=(specFreqRange./angFregTHz, nothing))

    idxPhase = idxRangeAboveThres(Ywp[3:end-2], thres)

    ph = unwrap(angle.(Ywp))
    gd = mydiff5(ph, w0)
    lines!(axr, ww[3:end-2][idxPhase]/angFregTHz, gd[idxPhase] / femto,
           label="Group Delay", color = colors[2], linestyle=:solid, linewidth=2)
    # axislegend(axr, position=:rt, framevisible = false)


    hidespines!(axr)
    hideydecorations!(axr, ticks=false, ticklabels = false, label=false)
    return (axl, axr)
end


function plotWavelengthDomain(f, ww, Ywp, specFreqRange; thres=0.01)
    @assert length(ww) == length(Ywp)

    wMask = ww .>= specFreqRange[1] .&& ww .<= specFreqRange[2]
    ww = ww[wMask]
    @assert(ww[1] > 0.0)
    Ywp = Ywp[wMask]

    # w0 = ww[2] - ww[1]
    ll = c2p ./ ww
    Sw = abs2.(Ywp)
    Sl = Sw .* ww.^2 / c2p

    (llGrid, SlGrid) = regGridInterp(reverse(ll), reverse(Sl))

    axl = Axis(f; commonAxAttr..., commonLeftAttr...,
               ylabel="Intensity (a.u.)",
               xlabel="Wavelength in nm",
               title="Pulse on wavelength axis",
               limits = ( (llGrid[1], llGrid[end] )./nano, nothing))
    lines!(axl, llGrid./nano, SlGrid./maximum(SlGrid), label="Intensity",
           color = colors[1], linestyle=:solid, linewidth=2)
    # axislegend(axl, position=:lt, framevisible = true)

    ph = unwrap(angle.(Ywp))
    idxPhase = idxRangeAboveThres(Ywp, thres)
    (llGrid, phGrid) = regGridInterp(reverse(ll[idxPhase]), reverse(ph[idxPhase]))
    axr = Axis(f; commonAxAttr..., commonRightAttr...,
               ylabel="Phase in rad",
               limits=((llGrid[1], llGrid[end] )./nano, extrema(phGrid)))
    lines!(axr, llGrid./nano, phGrid,
           label="Phase", color = colors[2], linestyle=:solid, linewidth=2)
    # axislegend(axr, position=:rt, framevisible = true)


    hidespines!(axr)
    hideydecorations!(axr, ticks=false, ticklabels = false, label=false)

    return (axl, axr)
end





"""
 Plot of amplitude and instantenous frequency in time domain on fs axis
 # Arguments
 - `f`: the grid position of the plot
 - `tt`: time vector, centered around zero
 - `ytana`: complex analytic signal
 - `timeLims`: time range which is displayed in ampl plot
 - `freqLims`: freq range which is displayed in inst freq plot
 # Optional Arguments
 - `thres (default=0.01)`: Only where amplitude exceeds this fraction of maximum, phasi ingo is displayed
"""
function plotTimeDomain(f, tt, ytana, timeLims, freqLims; thres=0.01, showRealField=false, wp=0)
    @assert length(tt) == length(ytana)
    Ts = tt[2] - tt[1]

    pht = unwrap(angle.(ytana))
    omt = mydiff5(pht, Ts) .+ wp
    idxPhase = idxRangeAboveThres(ytana[3:end-2], thres)

    axl = Axis(f[1,1]; commonAxAttr..., commonLeftAttr...,
               ylabel="Intensity (a.u.)",
               xlabel="Time in fs", title="Pulse on time axis",
               limits=(timeLims./femto, nothing))
    lines!(axl, tt/femto, abs.(ytana), label="Amplitude",
           color = colors[1], linestyle=:solid, linewidth=2)
    if showRealField
        lines!(axl, tt/femto, real.(ytana), label="Field",
               color = colors[3], linestyle=:solid, linewidth=1)
    end
    # axislegend(axl, position=:lt, framevisible = false)

    axr = Axis(f[1,1]; commonAxAttr..., commonRightAttr...,
               ylabel="Inst freq in THz",
               xlabel="Time in fs") #,
               # limits=(timeLims./femto, freqLims./angFregTHz))
    lines!(axr, tt[3:end-2][idxPhase]/femto, omt[idxPhase] / angFregTHz, label="Inst. freq",
           color = colors[2], linestyle=:solid, linewidth=2)
    # axislegend(axr, position=:rt, framevisible = false)

    hidespines!(axr)
    hidexdecorations!(axr)
    hideydecorations!(axr, ticks=false, ticklabels = false, label=false)

    return(axl, axr)
end


function plotFreqSpectrogram(f, delayVec, wVec, intensityMat;
                             delayRange=nothing, wRange=nothing, unitMax = true)
    if !isequal(size(intensityMat), (length(delayVec), length(wVec)))
        error("Size of input arguments inconsistent")
    end
    if unitMax
        intensityMat ./= maximum(intensityMat)
    end

    dIdx = 1:size(intensityMat, 1)
    if !isnothing(delayRange)
        dIdx = findall(delayVec .>= delayRange[1] .&& delayVec .<= delayRange[2])
    end
    wIdx = 1:size(intensityMat, 2)
    if !isnothing(wRange)
        wIdx = findall(wVec .>= wRange[1] .&& wVec .<= wRange[2])
    end
    smallSpecMat = intensityMat[dIdx, wIdx]
    glHM = f[1,1] = GridLayout()
    axHM = Axis(glHM[1,1], title="Spec-Mat",
                xlabel="Delay in fs", ylabel="Frequency in THz")
    hm = heatmap!(axHM, delayVec[dIdx]/femto, wVec[wIdx]/angFregTHz, smallSpecMat,
                  colormap = :nipy_spectral, interpolate=true)
    Colorbar(glHM[1,2], hm)
    colgap!(glHM, 8)
    return axHM
end


function plotWavelengthSpectrogram(f, delayVec, lVec, intensityMat;
                                   delayRange=nothing, lRange=nothing, unitMax = true,
                                   logScale=false, logThres=-5, interpolate=true)
    if !isequal(size(intensityMat), (length(delayVec), length(lVec)))
        error("Size of input arguments inconsistent")
    end
    if unitMax
        # wenn man hier kein copy macht, wird die Matrix im Scope des callers geändert!!
        intensityMat = copy(intensityMat) ./ maximum(intensityMat)
    end

    dIdx = 1:size(intensityMat, 1)
    if !isnothing(delayRange)
        dIdx = findall(delayVec .>= delayRange[1] .&& delayVec .<= delayRange[2])
    end

    lIdx = 1:size(intensityMat, 2)
    if !isnothing(lRange)
        lIdx = findall(lVec .>= lRange[1] .&& lVec .<= lRange[2])
    end
    smallSpecMat = intensityMat[dIdx, lIdx]
    if logScale
        smallSpecMat = log10.(smallSpecMat)
        smallSpecMat[smallSpecMat.<logThres] .= logThres
    end

    glHM = f[1,1] = GridLayout()
    axHM = Axis(glHM[1,1], title="Spec-Mat",
                xlabel="Delay in fs", ylabel="Wavelength in nm",
                titlesize=12, xlabelsize=10, ylabelsize=10,
                xticklabelsize=8, yticklabelsize=8, ylabelvisible=false, xlabelvisible=false)
    hm = heatmap!(axHM, delayVec[dIdx]/femto, lVec[lIdx]/nano, smallSpecMat,
                  colormap = :nipy_spectral, interpolate=interpolate)
    Colorbar(glHM[1,2], hm; ticklabelsize=8)
    colgap!(glHM, 8)
    return axHM
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
    commonAxAttr = (titlesize=10, xlabelsize=10, ylabelsize=10,
                    xticklabelsize=8, yticklabelsize=8, ylabelvisible=false, xlabelvisible=true,
                    xminorgridvisible=true, yminorgridvisible=true,
                    topspinevisible=false, leftspinecolor = colors[1], rightspinecolor = colors[2])
    commonLeftAttr = (ylabelcolor = colors[1], yticklabelcolor = colors[1], ytickcolor = colors[1])
    commonRightAttr = (ylabelcolor = colors[2], yticklabelcolor = colors[2], ytickcolor = colors[2],
                       yaxisposition = :right)
    legendAttr = (framevisible = false, labelsize=8)

    wLaser = pulseData.wCenter / 2
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
    instFreqLims =  (0.975, 1.025) .* wLaser

    axl = Axis(f1[1,1]; commonAxAttr..., commonLeftAttr...,
               ylabel="Intensity (a.u.)",
               xlabel="Time in fs",
               title="Time axis: FWHM=$(round(pulseData.fwhmT/femto, digits=2))fs",
               limits=(tPlotLims ./ femto, nothing))
    lines!(axl, delayAxis/femto, tempIntensity, label="Intensity",
           color = colors[1], linestyle=:solid, linewidth=2) # , marker=:x, markersize=4)
    #lines!(axl, delayAxis/femto, reverse(tempIntensity), label="Intensity-Rev", color = colors[4],
    #       linestyle=:solid, linewidth=2) # , marker=:x, markersize=4)
    axislegend(axl, position=:lt; legendAttr...)

    axr = Axis(f1[1,1]; commonAxAttr..., commonRightAttr...,
               ylabel="Inst. freq. in THz",
               limits=(tPlotLims ./ femto,  instFreqLims ./ angFregTHz))
               # limits=(tPlotLims ./ femto,  nothing))
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
    # llLims = (((-5, 5) .* pulseData.fwhmWv .+ 2*c2p /pulseData.wCenter)./nano, nothing)
    llLims = (extrema(ll)./nano, nothing)
    wvInt = abs.(Ywl).^2
    wvInt ./= maximum(wvInt)
    idxPhase = idxRangeAboveThres(wvInt, 1/100)
    phl = unwrap(angle.(Ywl))
    axl = Axis(f1[2,1]; commonAxAttr..., commonLeftAttr...,
               title="Pulse on wavelength axis (in nm): FWHM=$(round(pulseData.fwhmWv/nano, digits=2))nm",
               limits=llLims)
    lines!(axl, ll/nano, wvInt , label="Spectral Intensity",
           color = colors[1], linestyle=:solid, linewidth=2) # , marker=:x, markersize=4)
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



    if length(dirPath) > 0
        @assert isdir(dirPath)
        save(joinpath(dirPath, "SimulatedData.pdf"), f1)
        # else
        #     display(f1)
    end

end



export plotFreqDomain, plotWavelengthDomain, plotTimeDomain, plotFreqSpectrogram, plotWavelengthSpectrogram, plotPulseData
