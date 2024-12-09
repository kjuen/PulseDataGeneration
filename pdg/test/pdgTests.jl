# Tests
using Test
using FFTW
using DSP: unwrap
using GLMakie

#using Pkg
#Pkg.activate("..")
#Pkg.develop(path="../../pdg")
using pdg

const TolSmall = 1e-14
const TOL = 1e-10

function runVariousTests()


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
        #@show fwhmt
        fwhmt2 = berechneFWHM(yt.^2, tt)
        fwhmw = berechneFWHM(abs.(Ywc), wwc)
        fwhmw2 = berechneFWHM(abs.(Ywc).^2, wwc)

        rmst = berechneRMSBreite(yt, tt)
        rmst2 = berechneRMSBreite(yt.^2, tt)
        rmsw = berechneRMSBreite(abs.(Ywc), wwc)
        rmsw2 = berechneRMSBreite(abs.(Ywc).^2, wwc)


        # fwhm ist sehr sensitiv bzgl. Diskretisierung!
        @test isapprox(fwhmt, sqrt(log(2)); atol=2/100)
        @test isapprox(fwhmw * fwhmt, 8*log(2); atol=2/100)
        @test isapprox(fwhmw2 * fwhmt2, 4*log(2); atol= 2/100)
        @test isapprox(rmst, 1/4; atol= TOL)
        @test isapprox(rmsw, 2; atol= TOL)
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


    @testset "berechneBreite2" begin
        # Das hier basiert auf dem Buch von Cohen, Example 1.6, Seite 16
        alpha = 13.58
        beta = 122.23
        gamma = -212.23
        w0 = 134.7
        phi(t) = gamma/3*t^3 + beta/2*t^2 + w0*t
        s(t) = (alpha/pi)^(1/4) * exp(-alpha/2*t^2 + 1im*phi(t))  # Gl 1.99

        N = 1024
        Ts = 1/100
        T0 = N*Ts
        tt = (0:(N-1))*Ts
        tShift = T0/3
        yt = s.(tt .- tShift)

        # f = lines(tt, abs.(yt) , label="", linewidth=2, linestyle=:solid, # color=:blue;
        #           axis = (; title = "Ein Plot", xlabel = "x", ylabel="y"))
        # display(f)


        w0 = 2*pi/T0
        ws = N * w0
        wwc = (0:(N-1))*w0 .- ws/2
        Ywc = Ts * fftshift(fft(yt))

        rmsw = berechneRMSBreite(abs.(Ywc), wwc)
        rmswTheo = sqrt((alpha^2 + beta^2)/(2*alpha) + gamma^2 / (2*alpha^2))  # Gl (1.103)
        @test abs(rmsw - rmswTheo)/rmswTheo < 1e-7

        # f = lines(wwc, abs.(Ywc) , label="", linewidth=2, linestyle=:solid, # color=:blue;
        #            axis = (; title = "Ein Plot", xlabel = "x", ylabel="y"))
        # display(f)



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

    @testset "Sw-Sl conversion" begin

        N = 100
        a = 1/25
        ww = range(50, 100, N)
        wp = sum(ww)/N
        Sw = exp.(-a*(ww .- wp).^2)
        (ll, Sl) = Sw2Sl(ww, Sw)
        (ww2, Sw2) = Sl2Sw(ll, Sl)
        # f = lines(ww , Sw, label="", linewidth=2, linestyle=:solid, #color=:blue;
        #            axis = (; title = "Ein Plot", xlabel = "x", ylabel="y"))
        # lines!(ww2, Sw2, label="", linewidth=2, linestyle=:solid)
        # axislegend(position=:lt)
        # display(f)

        mse = sqrt(1/N * sum((Sw .- Sw2).^2))
        A = 1/N * sum(Sw)
        @test mse / A < 0.01



    end


    @testset "Intensity Matrix conversion" begin

        N = 100
        w = range(445, 455, N) * angFregTHz
        wc = sum(w)/N
        wr = w'
        y = range(64, 74, N)
        Mw = @. exp(-1/8*((wr-wc)^2+(y-70)^2))
        (l, Ml) = intensityMatFreq2Wavelength(w, Mw)
        (w2, Mw2) = intensityMatFreq2Wavelength(l, Ml)
        @test maximum(abs.(w2 .- w)) ./ maximum(w) < TolSmall
        # @show maximum(abs.(Mw2 .- Mw)) < TolSmall


    end


end
runSHGTests()
runVariousTests()
