using Plots
using DelimitedFiles
theme(:default)

function circle(center=[0, 0], radius=1)
    theta = 0:0.1:2*π
    x = center[1] .+ radius .* cos.(theta)
    y = center[2] .+ radius .* sin.(theta)
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
        yl = [0, lim + 100]
        plot(xlims=xl .* 2 .- [lim, 0], ylims=yl .* 1.8,
            # title="$i",
            markerstrokewidth=0,
            c=:royalblue,
            aspectratio=:equal,
            legend=nothing,
            axis=nothing,
            framestyle=:box)
        plot!(x, y, seriestype=:scatter, markerstrokewidth=0, color=:royalblue)
        lineb = linesBound[i]
        Db = split(lineb, " ")[begin:end-1] .|> x -> parse(Float64, x)
        xb, yb = Db[1:2:end], Db[2:2:end]
        plot!(xb, yb, seriestype=:scatter, markerstrokewidth=0, c=:orange)

        linetb = linesTopBound[i]
        Dtb = split(linetb, " ")[begin:end-1] .|> x -> parse(Float64, x)
        xtb, ytb = Dtb[1:2:end], Dtb[2:2:end]
        plot!(xtb, ytb, seriestype=:scatter, markerstrokewidth=0, c=:red)
        # hline!([lim - 50])

        linets = linesSequence[i]
        Dts = split(linets, " ")[begin:end-1] .|> x -> parse(Float64, x)
        xts, yts = Dts[1:2:end], Dts[2:2:end]
        plot!(xts, yts, c=:magenta, lw=3, alpha=0.9)

        plot!(circle([lim / 2, lim + 100], 250), color=:black, lw=2)
    end
    close(fall)
    close(fbound)
    gif(anim, "tplot.gif", fps=30)
    # mp4(anim, "tplot.mp4", fps=60)
end

begin
    t1 = @elapsed run(`make`)
    println("Compilation: $(t1)s")
    t2 = @elapsed run(`./Cell_Migration`)
    println("Running: $(t2)s")
    t3 = @elapsed generateAnimationVariableLength(30)
    println("Plotting: $(t3)s")
end

# (function testLinePlot()
#     fsbound = open("boundarySequenceOutData.csv", "r")
#     linesSequence = readlines(fsbound)

#     line = linesSequence[end]
#     D = split(line, " ")[begin:end-1] .|> x -> parse(Float64, x)
#     x, y = D[1:2:end], D[2:2:end]
#     println(length(x), length(y))
#     plot(x, y)
# end)()

# (function testCirclePlot()
#     fall = open("outdata.csv", "r")
#     lines = readlines(fall)
#     line = lines[end]
#     D = split(line, " ")[begin:end-1] .|> x -> parse(Float64, x)
#     x, y = D[1:2:end], D[2:2:end]
#     lim = 1200
#     xl = [0, lim]
#     yl = [0, lim+400]
#     plot(x, y, xlims=xl, ylims=yl,
#         aspectratio=:equal,
#         legend=nothing, seriestype=:scatter)

# end)()
#  nirgov interaction force
# Hv(x) = float(x > 0.0)
# function getforce(r)
#     rt = 70
#     if r > rt
#         return 0.0
#     end
#     U0 = 2650
#     U1 = 30
#     U2 = 2
#     U3 = 1
#     A0 = 8
#     A1 = 2
#     A2 = 25
#     A3 = 26
#     force = 0
#     force += U0 * r * exp(-(((r / A0)^2)))
#     force += U2 * exp(-r / A2)
#     force -= U3 * (r - A3)^2 * Hv(r - A3)
#     force += U1 * (r - A1) * Hv(r - A1)
#     return force
# end


# function boundaryattractiveForce(r)
#     req = 70
#     maxrep = 3000
#     maxadh = -1200
#     rmax = 70 * 2

#     if (r < req)
#         Fmag = ((maxadh - maxrep) / req) * r + maxrep
#     elseif (r < rmax)
#         Fmag = ((-maxadh) / (rmax - req)) * (r - req) + (maxadh)
#     else
#         Fmag = 0
#     end
#     return Fmag
# end