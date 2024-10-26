# Tests
using Test
using FFTW
using DSP: unwrap
using GLMakie

using pdg



function runVariousTests()
    tolSmall = 1e-14
    TOL = 1e-10

    N = 4*128
    Ts = 0.1 / 2
    T0 = N*Ts
    tt = (0:(N-1))*Ts
    tShift = T0/5
    yt = exp.(-4*(tt .- tShift).^2)


    w0 = 2*pi/T0
    ws = N * w0
    wwc = (0:(N-1))*w0 .- ws/2
    Ywc = fftshift(fft(yt))
    amplMask = abs.(Ywc) .> maximum(abs.(Ywc)) / 1000
    # f = Figure()
    # ax1 = Axis(f[1, 1], title = "Zeitbereich", xlabel = "t", ylabel="y(t)")
    # lines!(ax1, tt, yt, label="y(t)", linewidth=2, linestyle=:solid) #, color=:blue;)
    # ax2 = Axis(f[2, 1], title = "Frequenzbereich", xlabel = "w", ylabel="Y(w)")
    # lines!(ax2, wwc , abs.(Ywc), label="Y(w)", linewidth=2, linestyle=:solid)
    # display(f)

    # f = lines(wwc[amplMask],
    #           unwrap(angle.(Ywc))[amplMask], label="", linewidth=2, linestyle=:solid, # color=:blue;
    #           axis = (; title = "Ein Plot", xlabel = "x", ylabel="y"))
    # display(f)


    @testset "Breiten-Berechnungen" begin

        fwhmt = berechneFWHM(yt, tt)
        fwhmt2 = berechneFWHM(yt.^2, tt)
        fwhmw = berechneFWHM(abs.(Ywc), wwc)
        fwhmw2 = berechneFWHM(abs.(Ywc).^2, wwc)

        rmst = berechneRMSBreite(yt, tt)
        rmst2 = berechneRMSBreite(yt.^2, tt)
        rmsw = berechneRMSBreite(abs.(Ywc), wwc)
        rmsw2 = berechneRMSBreite(abs.(Ywc).^2, wwc)


        # fwhm ist sehr sensitiv bzgl. Diskretisierung!
        @test isapprox(fwhmw * fwhmt, 8*log(2); atol=2/100)
        @test isapprox(fwhmw2 * fwhmt2, 4*log(2); atol= 2/100)
        @test isapprox(rmsw * rmst, 1/2; atol= TOL)
        @test isapprox(rmsw2 * rmst2, 1/4; atol = TOL)

    end

    @testset "mseUpTo" begin

        msePhase = mseUpToSignAndLin(unwrap(angle.(Ywc[amplMask])), zeros(sum(amplMask)))[1]
        @test abs(msePhase) < TOL

        v1 = 1:10
        v2 = 2.2 * v1
        @test mseUpToScale(v1, v2) < TOL

        A = rand(100);
        B = copy(A);
        @test mseUpToScale(A, 0.24 * B) < TOL


        v1 = rand(10)
        v2 = -copy(v1) .+ 2.1 * (1:10) .+ 7.4
        @test mseUpToSignAndLin(v1, v2)[1] < TOL


    end

    @testset "berechneCOM" begin
        @test abs(berechneCOM(yt, tt) - tShift) < TOL
    end


    @testset "traceError" begin
        A = rand(100, 100)
        B = 2.34 * copy(A)
        @test traceError(A, B) < TOL

        s = 1/100
        B = 2.34 * copy(A) + s * randn(size(A))
        #  Der traceError ist ungefaehr s.
        # Durch die Skalierung 2.34 sollte der traceError hier aber sicher < s sein.
        @test traceError(A, B) < s
    end
end


function runSHGTests()

    # Amplitude im Frequenzbereich: ein einfacher Gauss-Puls ohne Phase
    A(w, wp, sw) = exp.(-1/(2*sw^2) * (w .- wp).^2)
    ytaTheo(t, wp, sw) = sw / sqrt(2*pi) * exp.(-1/2*sw^2 * t.^2 + 1im * wp * t)

    # SHG-Frog-Amplitude theoretisch:

    # exakte Loesung der Frog-Amplitude
    # Siehe @. Macro, um das auf einem Grid auszuwerten,
# z.B. https://docs.juliaplots.org/stable/series_types/contour/
    Afrog(tau, w, wp, sw) = sw/(2*sqrt(pi)) * exp(-1/4*sw^2 * tau^2 - 1/(4*sw^2) * (w - 2*wp)^2  - 1/2 * 1im *w * tau)





    @testset "createSHG: Gauss 1" begin
        wCenter = 110.0
        N = 128  # angestrebte Grid-Größe
        wp = 100.0 # Traeger-Frequenz
        sw = 4.0   # Breite im Frequenzbereich
        delay = 0.05  # = Abtastzeit
        idxVecCentered = ((-N ÷ 2):(N ÷ 2 - 1))
        delayAxis = idxVecCentered * delay

        ws = 2*pi/delay
        w0 = ws/N
        wAxis = idxVecCentered * w0 .+ wCenter
        wAxisShifted = wAxis .+ wCenter

        yta = ytaTheo(delayAxis, wp, sw)
        Afft = createSHGAmplitude(yta, delay, wCenter)
        # f = Figure()
        # plotFreqSpectrogram(f[1,1], delayAxis, wAxisShifted, abs.(Afft).^2)
        # display(f)

        AfrogMat = @. Afrog(delayAxis, wAxisShifted', wp, sw)

        AdiffAbs = abs.(Afft) .- abs.(AfrogMat)
        Adiff = Afft .- AfrogMat
        @test maximum(abs.(Adiff)) < 1e-12

    end

    @testset "createSHG: Gauss 2" begin
        wCenter = 110.0
        N = 128  # angestrebte Grid-Größe
        wp = 100.0 # Traeger-Frequenz
        sw = 8.0   # Breite im Frequenzbereich
        delay = 0.03  # = Abtastzeit
        idxVecCentered = ((-N ÷ 2):(N ÷ 2 - 1))
        delayAxis = idxVecCentered * delay

        ws = 2*pi/delay
        w0 = ws/N
        wAxis = idxVecCentered * w0 .+ wCenter
        wAxisShifted = wAxis .+ wCenter

        yta = ytaTheo(delayAxis, wp, sw)
        Afft = createSHGAmplitude(yta, delay, wCenter)
        # f = Figure()
        # plotFreqSpectrogram(f[1,1], delayAxis, wAxisShifted, abs.(Afft).^2)
        # display(f)

        AfrogMat = @. Afrog(delayAxis, wAxisShifted', wp, sw)

        AdiffAbs = abs.(Afft) .- abs.(AfrogMat)
        Adiff = Afft .- AfrogMat
        @test maximum(abs.(Adiff)) < 1e-12

    end

     @testset "createSHG: Gauss 3" begin
        wCenter = 102.0
        N = 128*2  # angestrebte Grid-Größe
        wp = 100.0 # Traeger-Frequenz
        sw = 2.0   # Breite im Frequenzbereich
        delay = 0.04  # = Abtastzeit
        idxVecCentered = ((-N ÷ 2):(N ÷ 2 - 1))
        delayAxis = idxVecCentered * delay

        ws = 2*pi/delay
        w0 = ws/N
        wAxis = idxVecCentered * w0 .+ wCenter
        wAxisShifted = wAxis .+ wCenter

        yta = ytaTheo(delayAxis, wp, sw)
        Afft = createSHGAmplitude(yta, delay, wCenter)
        # f = Figure()
        # plotFreqSpectrogram(f[1,1], delayAxis, wAxisShifted, abs.(Afft).^2)
        # display(f)

        AfrogMat = @. Afrog(delayAxis, wAxisShifted', wp, sw)

        AdiffAbs = abs.(Afft) .- abs.(AfrogMat)
        Adiff = Afft .- AfrogMat
        @test maximum(abs.(Adiff)) < 1e-10

    end


end
runSHGTests()
runVariousTests()
