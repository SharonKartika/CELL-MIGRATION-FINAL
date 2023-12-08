using Plots
using DelimitedFiles
theme(:dark)

function circle(center=[0,0], radius=1)
    theta = 0:0.1:2*π
    x = center[1].+radius .* cos.(theta)
    y = center[2].+radius .* sin.(theta)
    x, y
end

function generateAnimation(nskip=1)
    data = readdlm("outdata.csv", ' ')[:, begin:end-1] .|> float
    anim = @animate for i = 1:nskip:size(data)[1]
        x, y = data[i, 1:2:end], data[i, 2:2:end]
        lim = 1200 * 2
        xl = [-lim, lim]
        yl = [-lim, lim]
        plot(x, y, seriestype=:scatter,
            xlims=xl, ylims=yl,
            title="$i",
            markerstrokewidth=0,
            c=:royalblue,
            aspectratio=:equal,
            legend=nothing,
            axis=nothing,
            framestyle=:box)
    end
    gif(anim, "tplot.gif", fps=30)
end

function generateAnimationVariableLength(nskip=1)
    fall = open("outdata.csv", "r")
    fbound = open("boundaryOutData.csv", "r")
    ftbound = open("topBoundaryOutData.csv", "r")
    fsbound = open("boundarySequenceOutData.csv", "r")
    lines = readlines(fall)
    linesBound = readlines(fbound)
    linesTopBound = readlines(ftbound)
    linesSequence = readlines(fsbound)
    anim = @animate for i = 1:nskip:length(lines)
        line = lines[i]
        D = split(line, " ")[begin:end-1] .|> x -> parse(Float64, x)
        x, y = D[1:2:end], D[2:2:end]
        lim = 1200
        xl = [0, lim]
        yl = [0, lim+100]
        plot(xlims=xl, ylims=yl,
            title="$i",
            markerstrokewidth=0,
            c=:royalblue,
            aspectratio=:equal,
            legend=nothing,
            axis=nothing,
            framestyle=:box)
        plot!(x, y, seriestype=:scatter, color=:royalblue)
        lineb = linesBound[i]
        Db = split(lineb, " ")[begin:end-1] .|> x -> parse(Float64, x)
        xb, yb = Db[1:2:end], Db[2:2:end]
        plot!(xb, yb, seriestype=:scatter, c=:orange)

        linetb = linesTopBound[i]
        Dtb = split(linetb, " ")[begin:end-1] .|> x -> parse(Float64, x)
        xtb, ytb = Dtb[1:2:end], Dtb[2:2:end]
        plot!(xtb, ytb, seriestype=:scatter, c=:red)
        # hline!([lim - 50])

        linets = linesSequence[i]
        Dts = split(linets, " ")[begin:end-1] .|> x -> parse(Float64, x)
        xts, yts = Dts[1:2:end], Dts[2:2:end]
        plot!(xts, yts, c=:skyblue)

        plot!(circle([lim/2, lim+100], 185), color=:yellow)
    end
    close(fall)
    close(fbound)
    # gif(anim, "tplot.gif", fps=30)
    mp4(anim, "tplot.mp4", fps=60)
end

begin
    t1 = @elapsed run(`make`)
    println("Compilation: $(t1)s")
    t2 = @elapsed run(`./Cell_Migration`)
    println("Running: $(t2)s")
    t3 = @elapsed generateAnimationVariableLength(30)
    println("Plotting: $(t3)s")
end

(function testLinePlot()
    fsbound = open("boundarySequenceOutData.csv", "r")
    linesSequence = readlines(fsbound)

    line = linesSequence[end]
    D = split(line, " ")[begin:end-1] .|> x -> parse(Float64, x)
    x, y = D[1:2:end], D[2:2:end]
    println(length(x), length(y))
    plot(x, y)
end)()

(function testCirclePlot()
    fall = open("outdata.csv", "r")
    lines = readlines(fall)
    line = lines[end]
    D = split(line, " ")[begin:end-1] .|> x -> parse(Float64, x)
    x, y = D[1:2:end], D[2:2:end]
    lim = 1200
    xl = [0, lim]
    yl = [0, lim+400]
    plot(x, y, xlims=xl, ylims=yl,
        aspectratio=:equal,
        legend=nothing, seriestype=:scatter)

end)()


