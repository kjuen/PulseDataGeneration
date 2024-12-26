using FreqTables
using DelimitedFiles
using DSP: unwrap
using DataFrames
using CSV

# include("/home/kjuen/dvlp/kai/Julia/PulseDataGeneration/pdgFunctionsAndConstants.jl")

# Matrix einlesen mit
# mat = readdlm(fn, '\t', Float32; skipstart=1, header=false);

# Spalten fuer die Dateien mit 5 Spalten (werden aber nicht in die Dateien geschrieben)
# const timeColNames=["Delay", "Intensity", "Phase", "Real", "Imag"]
# const wavColNames=["Wavelength", "Intensity", "Phase", "Real", "Imag"]


#* Spectrogram Header Info Datenstruktur
struct SpectrogramInfo
    NumDelays::Int64          # Number of delays (rows)
    NumSpecPoints::Int64      # Number of points in spectrum (columns)
    DelayDist::Float64        # Distance between delays
    ResSpec::Float64          # Resolution of spectrum
    CenterWavelen::Float64    # Center Wavelength
    function SpectrogramInfo(headerString::AbstractString)
        items = collect(eachsplit(headerString))
        if length(items) != 5
            error("Headerstring besteht nicht aus 5 Werten")
        end
        new(parse(Int64, items[1]),
            parse(Int64, items[2]),
            parse(Float64, items[3])*femto,  # given in fs
            parse(Float64, items[4])*nano,   # given in nm
            parse(Float64, items[5])*nano)   # given in nm
    end

    function SpectrogramInfo(NumDelays, NumSpecPoints, DelayDist, ResSpec, CenterWavelen)
        new(NumDelays, NumSpecPoints, DelayDist, ResSpec, CenterWavelen)
    end


    function SpectrogramInfo(delayAxis::AbstractVector,
                             wavAxis::AbstractVector)
        de = extrema(diff(delayAxis))
        if ( (de[2] - de[1]) / de[2] ) > 1e-8
            error("delayAxis not equidistant")
        end
        we = extrema(diff(wavAxis))
        if ( (we[2] - we[1]) / we[2] ) > 1e-8
            error("wavAxis not equidistant")
        end
        # calculate center wavelength
        lenWav = length(wavAxis)
        cw = 0.0
        if isodd(lenWav)
            cw = wavAxis[(lenWav+1)÷2]
        else
            cw = wavAxis[(lenWav÷2)+1]
        end

        new(length(delayAxis),
            length(wavAxis),
            1/2*(de[1] + de[2]),
            1/2*(we[1] + we[2]),
            cw)
    end
end

# Liest Header-Info aus Datei ein. Die Header-Info kann über eine oder über zwei Zeilen verteilt sein.
function readSpectrogramInfoFromFile(fileName::AbstractString)::Tuple{SpectrogramInfo, Int}
    headerString = ""
    for (idx, val) in enumerate(eachline(fileName))
        if idx > 2
            error("Unknown header format: more than 2 lines")
        end
        headerString *= (" " * replace(val, ","=>"."))
        items = collect(eachsplit(headerString))
        if length(items) == 5  # complete header line
            return (SpectrogramInfo(headerString), idx)
        end
    end
end


showFieldFreqs(symbol, SpecInfo) = freqtable(map(hi->getproperty(hi, symbol), values(SpecInfo)))



function generateAxis(N, res, center=0.0)
    idx = ceil(Int, -N/2) : 1 : floor(Int, (N-1)/2)
    @assert length(idx) == N
    return idx * res .+ center
end

function generateAxes(hi::SpectrogramInfo)
    delayAxis = generateAxis(hi.NumDelays, hi.DelayDist)
    @assert delayAxis[1] < delayAxis[end]
    wavAxis = generateAxis(hi.NumSpecPoints, hi.ResSpec, hi.CenterWavelen)
     #@assert wavAxis[1] < wavAxis[end]
    return (delayAxis, wavAxis)
end

function getExpectedFreqResolutionFromDelay(hi::SpectrogramInfo)
    return delayToSpecRes(hi.DelayDist, hi.NumSpecPoints, c2p / hi.CenterWavelen)

end

function print(hi::SpectrogramInfo; digits=2)
    ws = 2*pi / hi.DelayDist   # DelayDist ist in SI-Einheit angegeben
    w0 = ws / hi.NumSpecPoints
    wCenter = c2p / hi.CenterWavelen
    N = hi.NumSpecPoints
    idxVec = (-N ÷ 2):(N ÷ 2 - 1)
    fftAxis = idxVec * w0 .+ wCenter
    fftAxisRange = extrema(fftAxis)

    println("Distance between delays: $(round(hi.DelayDist/femto, digits=digits))fs")

    (delayAxis, wavAxis) = generateAxes(hi)
    println("Number of delays (rows): $(hi.NumDelays) from $(round(delayAxis[1]/femto, digits=digits))fs to $(round(delayAxis[end]/femto, digits=digits))fs")
    println("FFT-Frequency axis: $(round(fftAxisRange[1]/angFregTHz, digits=2))THz to $(round(fftAxisRange[2]/angFregTHz, digits=2))THz with resolution $(round(w0/angFregTHz, digits=2))THz")

    (wres, lres) = getExpectedFreqResolutionFromDelay(hi)
    @show wres / angFregTHz
    @show lres / nano
    println("Resolution of spectrum: $(round(hi.ResSpec / nano, digits=digits))nm")
    println("Center Wavelength : $(round(hi.CenterWavelen / nano, digits=digits ))nm")
    println("Center angular frequency : $(round( wCenter / angFregTHz, digits=digits ))THz")
    println("Number of points in spectrum (columns): $(hi.NumSpecPoints)")
    println("Wavelength range: $(round(wavAxis[1]/nano, digits=digits))nm to $(round(wavAxis[end]/nano, digits=digits))nm")
    println("Angular frequence range: $(round(c2p / wavAxis[end] / angFregTHz, digits=digits))THz to $(round(c2p / wavAxis[1]/angFregTHz, digits=digits))THz")
end


function loadSpecMat(fn::String)
    (hi, offset) = readSpectrogramInfoFromFile(fn)

    mat = readdlm(fn, '\t', Float64; skipstart=offset, header=false);
    return (hi, mat)
end


function createHeaderline(hi::SpectrogramInfo; digits=4)
    return "$(hi.NumDelays) $(hi.NumSpecPoints) $(round(hi.DelayDist / femto, digits=4)) $(round(hi.ResSpec / nano, digits=4)) $(round(hi.CenterWavelen / nano, digits=3))\n"

end



function writeFrogTraceToFile(fn, delayAxis, wavAxis, intensityMat;
                              maxVal=maximum(intensityMat))
    hi = SpectrogramInfo(delayAxis, wavAxis)
    io = open(fn, "w")
    write(io, createHeaderline(hi))
    maxI = maximum(intensityMat)
    # intensityMat .*= maxVal/maxI
    intensityMat = round.(intensityMat*maxVal/maxI, digits=4)

    writedlm(io, intensityMat, '\t')
    close(io)

end


struct ComplexSignal{T <: Real}
    xVals::AbstractVector{T}
    ampl::Vector{T}
    phase::Vector{T}
    function ComplexSignal{T}(xv::AbstractVector{T},
                              csig::Vector{Complex{T}}) where {T<:Real}
        if length(xv) != length(csig)
            error("arguments must have equal length")
        end
        new(xv, abs.(csig), unwrap(angle.(csig)))
    end
    function ComplexSignal{T}(xv::AbstractVector{T},
                              a::Vector{T}, p::Vector{T}) where {T<:Real}
        if (length(xv) != length(a)) || (length(xv) != length(p))
            error("all arguments must have equal length")
        end
        new(xv, a, p)
    end

end
# Dies hier ist noetig, sonst klappt es nicht,
# siehe https://docs.julialang.org/en/v1/manual/constructors/#Parametric-Constructors
ComplexSignal(xv::AbstractVector{T}, csig::Vector{Complex{T}}) where {T<:Real} = ComplexSignal{T}(xv, csig)
ComplexSignal(xv::AbstractVector{T}, a::Vector{T}, p::Vector{T}) where {T<:Real} = ComplexSignal{T}(xv, a, p)

import Base: length
length(cs::ComplexSignal{T}) where{T<:Real} = length(cs.ampl)
function getSignal(cs::ComplexSignal{T})::Vector{Complex{T}} where{T<:Real}
    return cs.ampl .* exp.(1im * cs.phase)
end


function writeFiveColumns(fileName::String, csig::ComplexSignal{T}; delim="  ") where{T<:Real}
    c = csig.ampl .* exp.(1im * csig.phase)
    CSV.write(fileName,
              (x= csig.xVals,
               i= csig.ampl.^2,
               p= csig.phase,
               r= real.(c),
               im= imag.(c));
              delim=delim,
              writeheader=false,
              decimal=',')
end


function readFrogDat(frogFileName)

    df = DataFrame(MinFrogErr=Float64[],
                   TempFWHM=Float64[],  # in fs
                   SpecFWHM=Float64[],  # in nm
                   AutocorrFWHM=Float64[], # in fs
                   TBP_FWHM = Float64[],
                   TBP_RMS = Float64[],
                   TBP_TempLap = Float64[],
                   TPB_SpecLap = Float64[])

    MinFrogErr = NaN
    TemporalFWHM = NaN
    SpectralFWHM = NaN
    AutocorrelationFWHM = NaN
    TimeBandwidthProductFWHM = NaN
    TimeBandwidthProductRMS = NaN
    TimeBandwidthProductTemporalLaplacian = NaN
    TimeBandwidthProductSpectralLaplacian = NaN

    for l in eachline(frogFileName)

        if match(r"minimum\s+FROG\s+error", l) !== nothing
            items = collect(eachsplit(l))
            MinFrogErr=parse(Float64, replace(items[end], ","=>"."))
        end


        if match(r"^Temporal\s+FWHM", l) !== nothing
            items = collect(eachsplit(l))
            @assert length(items) == 4
            @assert occursin("fs", items[4])
            TemporalFWHM = parse(Float64, replace(items[3], ","=>"."))
        end

        if match(r"^Spectral\s+FWHM", l) !== nothing
            items = collect(eachsplit(l))
            @assert length(items) == 4
            @assert occursin("nm", items[4])
            SpectralFWHM = parse(Float64, replace(items[3], ","=>"."))
        end

        if match(r"^Autocorrelation\s+FWHM", l) !== nothing
            items = collect(eachsplit(l))
            @assert length(items) == 4
            @assert occursin("fs", items[4])
            AutocorrelationFWHM = parse(Float64, replace(items[3], ","=>"."))
        end

        if match(r"^Time-Bandwidth\s+Product\s+-\s+FWHM", l) !== nothing
            items = collect(eachsplit(l))
            TimeBandwidthProductFWHM = parse(Float64, replace(items[end], ","=>"."))
        end

        if match(r"^Time-Bandwidth\s+Product\s+-\s+RMS", l) !== nothing
            items = collect(eachsplit(l))
            TimeBandwidthProductRMS = parse(Float64, replace(items[end], ","=>"."))
        end

        if match(r"^Time-Bandwidth\s+Product\s+-\s+Temporal\s+Laplacian", l) !== nothing
            items = collect(eachsplit(l))
            TimeBandwidthProductTemporalLaplacian = parse(Float64, replace(items[end], ","=>"."))
        end

        if match(r"^Time-Bandwidth\s+Product\s+-\s+Spectral\s+Laplacian", l) !== nothing
            items = collect(eachsplit(l))
            TimeBandwidthProductSpectralLaplacian = parse(Float64, replace(items[end], ","=>"."))
        end

    end
    push!(df, [MinFrogErr,
               TemporalFWHM,
               SpectralFWHM,
               AutocorrelationFWHM,
               TimeBandwidthProductFWHM,
               TimeBandwidthProductRMS,
               TimeBandwidthProductTemporalLaplacian,
               TimeBandwidthProductSpectralLaplacian])
    return df

end

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
                             asData[1], asData[2], asData[3][matName];
                             maxVal=1.0)
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




export SpectrogramInfo, readSpectrogramInfoFromFile, showFieldFreqs, generateAxes, print, loadSpecMat, writeFrogTraceToFile, ComplexSignal, writeFiveColumns, readFrogDat, writePulseDataToDir, writePulseFeatures
