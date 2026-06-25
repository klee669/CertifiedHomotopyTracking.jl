export PathBox, PathVisualization, path_boxes, export_path_tikz, export_path_obj

struct PathBox
    t::AcbFieldElem
    x::Vector{AcbFieldElem}
    source::Symbol
    metadata::NamedTuple
end

struct PathVisualization
    boxes::Vector{PathBox}
    axes::Tuple
    source::Symbol
    metadata::NamedTuple
end

path_boxes(viz::PathVisualization) = viz.boxes

function _path_box_with_field(box::PathBox, CC::AcbField)
    isempty(box.x) && return box
    parent(box.x[1]) === CC && return box
    return PathBox(CC(box.t), _convert_acb_vector(CC, box.x), box.source, box.metadata)
end
_path_box_with_field(box, ::AcbField) = box

function _normalize_path_axes(axes)
    axes === nothing && return (:t, 1)
    axes isa Tuple && return axes
    axes isa AbstractVector && return Tuple(axes)
    return (axes,)
end

function _path_axis_label(axis)
    axis === :t && return "\$t\$"
    axis isa Integer && return "\$\\mathrm{Re}(x_{$(Int(axis))})\$"
    axis isa Symbol && return String(axis)
    if axis isa Tuple && length(axis) == 2
        a, b = axis
        if a isa Integer
            part = b === :imag || b === :im ? "Im" : b === :abs ? "abs" : "Re"
            part == "abs" && return "\$|x_{$(Int(a))}|\$"
            return "\$\\mathrm{$part}(x_{$(Int(a))})\$"
        elseif b isa Integer
            part = a === :imag || a === :im ? "Im" : a === :abs ? "abs" : "Re"
            part == "abs" && return "\$|x_{$(Int(b))}|\$"
            return "\$\\mathrm{$part}(x_{$(Int(b))})\$"
        end
    end
    return string(axis)
end

function _path_axis_value(box::PathBox, axis)
    axis === :t && return box.t
    if axis isa Integer
        return real(box.x[Int(axis)])
    elseif axis isa Tuple && length(axis) == 2
        a, b = axis
        if a isa Integer
            1 <= Int(a) <= length(box.x) || throw(BoundsError(box.x, Int(a)))
            b === :imag || b === :im ? (return imag(box.x[Int(a)])) : nothing
            b === :abs && return abs(convert_to_double_int(box.x[Int(a)]))
            return real(box.x[Int(a)])
        elseif b isa Integer
            1 <= Int(b) <= length(box.x) || throw(BoundsError(box.x, Int(b)))
            a === :imag || a === :im ? (return imag(box.x[Int(b)])) : nothing
            a === :abs && return abs(convert_to_double_int(box.x[Int(b)]))
            return real(box.x[Int(b)])
        end
    end
    throw(ArgumentError("unsupported path visualization axis $axis. Use :t, an integer coordinate, or (i, :real/:imag)."))
end

_path_box_stage(box::PathBox) = haskey(box.metadata, :stage) ? box.metadata.stage : missing
_is_track_node_box(box::PathBox) =
    box.source === :track_path && _path_box_stage(box) in (:initial_refinement, :refinement, :final_refinement)

function _path_export_boxes(boxes)
    return PathBox[box for box in boxes if !_is_track_node_box(box)]
end

function _real_interval_bounds(value)
    value isa Real && return Float64(value), Float64(value)
    value isa Complex && return Float64(real(value)), Float64(real(value))
    a = value isa AcbFieldElem ? real(value) : value
    mid = Float64(Nemo.midpoint(a))
    rad = Float64(Nemo.radius(a))
    return mid - rad, mid + rad
end

function _path_axis_bounds(box::PathBox, axis)
    value = _path_axis_value(box, axis)
    value isa Real && return Float64(value), Float64(value)
    return _real_interval_bounds(value)
end

function _path_axis_bbox(boxes, axes)
    lows = fill(Inf, length(axes))
    highs = fill(-Inf, length(axes))
    for box in boxes
        for (i, axis) in enumerate(axes)
            lo, hi = _path_axis_bounds(box, axis)
            lows[i] = min(lows[i], lo)
            highs[i] = max(highs[i], hi)
        end
    end
    for i in eachindex(lows)
        if !isfinite(lows[i]) || !isfinite(highs[i])
            lows[i], highs[i] = 0.0, 1.0
        elseif lows[i] == highs[i]
            delta = max(1e-6, abs(lows[i]) * 1e-6)
            lows[i] -= delta
            highs[i] += delta
        end
    end
    return lows, highs
end

_path_scale(v, lo, hi, size) = size * (Float64(v) - lo) / (hi - lo)

function _path_rect(box::PathBox, axes, lows, highs, width, height)
    xlo, xhi = _path_axis_bounds(box, axes[1])
    ylo, yhi = _path_axis_bounds(box, axes[2])
    return (
        _path_scale(xlo, lows[1], highs[1], width),
        _path_scale(ylo, lows[2], highs[2], height),
        _path_scale(xhi, lows[1], highs[1], width),
        _path_scale(yhi, lows[2], highs[2], height),
    )
end

function _path_box_scaled_bounds(box::PathBox, axes, lows, highs, sizes)
    bounds = Tuple{Float64,Float64}[]
    for (i, axis) in enumerate(axes)
        lo, hi = _path_axis_bounds(box, axis)
        push!(
            bounds,
            (
                _path_scale(lo, lows[i], highs[i], sizes[i]),
                _path_scale(hi, lows[i], highs[i], sizes[i]),
            ),
        )
    end
    return bounds
end

function _path_point_scaled(point::PathBox, axes, lows, highs, sizes)
    return Tuple(_path_scale(first(_path_axis_bounds(point, axes[i])), lows[i], highs[i], sizes[i]) for i in eachindex(axes))
end

function _path_trace_points(input, trace)
    trace !== nothing && return trace
    hasproperty(input, :visualization) && return _path_trace_points(getproperty(input, :visualization), nothing)
    input isa PathVisualization || return PathBox[]
    haskey(input.metadata, :trace_points) || return PathBox[]
    return input.metadata.trace_points
end

function _write_tikz_preamble!(io, standalone, border)
    if standalone
        println(io, "\\documentclass[tikz,border=$border]{standalone}")
        println(io, "\\usepackage{tikz}")
        println(io, "\\begin{document}")
    end
    return nothing
end

function _write_tikz_trace_2d!(io, trace_points, axes_tuple, lows, highs; width, height, trace_color)
    length(trace_points) >= 2 || return nothing
    sizes = (Float64(width), Float64(height))
    coords = String[]
    for point in trace_points
        x, y = _path_point_scaled(point, axes_tuple, lows, highs, sizes)
        push!(coords, "($x,$y)")
    end
    println(io, "\\draw[$trace_color, line width=0.55pt] ", join(coords, " -- "), ";")
    return nothing
end

function _write_tikz_2d!(io, boxes, axes_tuple, lows, highs; width, height, color, trace_points=PathBox[], show_trace=false, trace_color="black")
    println(io, "\\begin{tikzpicture}[x=1cm,y=1cm]")
    println(io, "\\draw[->] (0,0) -- ($(width + 0.35),0) node[right] {", _path_axis_label(axes_tuple[1]), "};")
    println(io, "\\draw[->] (0,0) -- (0,$(height + 0.35)) node[above] {", _path_axis_label(axes_tuple[2]), "};")
    println(io, "\\node[below left] at (0,0) {$(round(lows[1]; sigdigits=4)), $(round(lows[2]; sigdigits=4))};")
    println(io, "\\node[below] at ($width,0) {$(round(highs[1]; sigdigits=4))};")
    println(io, "\\node[left] at (0,$height) {$(round(highs[2]; sigdigits=4))};")
    println(io, "\\draw[gray!25] (0,0) grid[xstep=$(width / 8),ystep=$(height / 6)] ($width,$height);")
    for box in boxes
        x1, y1, x2, y2 = _path_rect(box, axes_tuple, lows, highs, width, height)
        println(io, "\\draw[draw=$(color)!70!black, fill=$(color)!18, line width=0.25pt, fill opacity=0.28] ($x1,$y1) rectangle ($x2,$y2);")
    end
    show_trace && _write_tikz_trace_2d!(io, trace_points, axes_tuple, lows, highs; width=width, height=height, trace_color=trace_color)
    println(io, "\\end{tikzpicture}")
    return nothing
end

function _write_tikz_cuboid!(io, bounds, color)
    (x1, x2), (y1, y2), (z1, z2) = bounds
    faces = (
        ((x1, y1, z1), (x2, y1, z1), (x2, y2, z1), (x1, y2, z1)),
        ((x1, y1, z2), (x2, y1, z2), (x2, y2, z2), (x1, y2, z2)),
        ((x1, y1, z1), (x2, y1, z1), (x2, y1, z2), (x1, y1, z2)),
        ((x2, y1, z1), (x2, y2, z1), (x2, y2, z2), (x2, y1, z2)),
        ((x1, y2, z1), (x2, y2, z1), (x2, y2, z2), (x1, y2, z2)),
        ((x1, y1, z1), (x1, y2, z1), (x1, y2, z2), (x1, y1, z2)),
    )
    for face in faces
        coords = join(("($(p[1]),$(p[2]),$(p[3]))" for p in face), " -- ")
        println(io, "\\filldraw[fill=$(color)!16, draw=$(color)!65!black, fill opacity=0.16, line width=0.2pt] $coords -- cycle;")
    end
    return nothing
end

function _write_tikz_trace_3d!(io, trace_points, axes_tuple, lows, highs; width, height, depth, trace_color)
    length(trace_points) >= 2 || return nothing
    sizes = (Float64(width), Float64(depth), Float64(height))
    coords = String[]
    for point in trace_points
        x, y, z = _path_point_scaled(point, axes_tuple, lows, highs, sizes)
        push!(coords, "($x,$y,$z)")
    end
    println(io, "\\draw[$trace_color, line width=0.55pt] ", join(coords, " -- "), ";")
    return nothing
end

function _write_tikz_3d!(io, boxes, axes_tuple, lows, highs; width, height, depth, color, trace_points=PathBox[], show_trace=false, trace_color="black")
    sizes = (Float64(width), Float64(depth), Float64(height))
    println(io, "\\begin{tikzpicture}[x={(1cm,0cm)}, y={(0.45cm,0.28cm)}, z={(0cm,1cm)}]")
    println(io, "\\draw[->] (0,0,0) -- ($(sizes[1] + 0.45),0,0) node[right] {", _path_axis_label(axes_tuple[1]), "};")
    println(io, "\\draw[->] (0,0,0) -- (0,$(sizes[2] + 0.45),0) node[above right] {", _path_axis_label(axes_tuple[2]), "};")
    println(io, "\\draw[->] (0,0,0) -- (0,0,$(sizes[3] + 0.45)) node[above] {", _path_axis_label(axes_tuple[3]), "};")
    println(io, "\\node[below left] at (0,0,0) {$(round(lows[1]; sigdigits=4)), $(round(lows[2]; sigdigits=4)), $(round(lows[3]; sigdigits=4))};")
    println(io, "\\node[below] at ($(sizes[1]),0,0) {$(round(highs[1]; sigdigits=4))};")
    println(io, "\\node[right] at (0,$(sizes[2]),0) {$(round(highs[2]; sigdigits=4))};")
    println(io, "\\node[left] at (0,0,$(sizes[3])) {$(round(highs[3]; sigdigits=4))};")
    println(io, "\\draw[gray!25] (0,0,0) -- ($(sizes[1]),0,0) -- ($(sizes[1]),$(sizes[2]),0) -- (0,$(sizes[2]),0) -- cycle;")
    println(io, "\\draw[gray!20] (0,0,0) -- (0,0,$(sizes[3]));")
    println(io, "\\draw[gray!20] ($(sizes[1]),0,0) -- ($(sizes[1]),0,$(sizes[3]));")
    println(io, "\\draw[gray!20] (0,$(sizes[2]),0) -- (0,$(sizes[2]),$(sizes[3]));")
    println(io, "\\draw[gray!20] ($(sizes[1]),$(sizes[2]),0) -- ($(sizes[1]),$(sizes[2]),$(sizes[3]));")
    for box in boxes
        bounds = _path_box_scaled_bounds(box, axes_tuple, lows, highs, sizes)
        _write_tikz_cuboid!(io, bounds, color)
    end
    show_trace && _write_tikz_trace_3d!(io, trace_points, axes_tuple, lows, highs; width=width, height=height, depth=depth, trace_color=trace_color)
    println(io, "\\end{tikzpicture}")
    return nothing
end

function _as_path_boxes(input)
    input isa PathVisualization && return input.boxes
    input isa AbstractVector{PathBox} && return collect(input)
    hasproperty(input, :path_boxes) && return _as_path_boxes(getproperty(input, :path_boxes))
    if input isa AbstractVector
        return PathBox[box for box in input if box isa PathBox]
    end
    hasproperty(input, :boxes) && return _as_path_boxes(getproperty(input, :boxes))
    throw(ArgumentError("expected PathVisualization, TrackResult, PosterioriPathResult, or a vector of PathBox objects."))
end

path_boxes(input) = _as_path_boxes(input)

function export_path_tikz(
    input,
    filename::AbstractString;
    axes=nothing,
    width=12.0,
    height=7.0,
    depth=6.0,
    color="blue",
    standalone=true,
    border="4pt",
    show_trace=false,
    trace=nothing,
    trace_color="black",
)
    boxes = _path_export_boxes(_as_path_boxes(input))
    isempty(boxes) && throw(ArgumentError("no path boxes available to visualize. Run with store_boxes=:full or visualize=true."))
    trace_points = PathBox[_path_trace_points(input, trace)...]
    all_boxes = show_trace ? PathBox[boxes; trace_points] : boxes
    axes_tuple = _normalize_path_axes(axes)
    length(axes_tuple) in (2, 3) || throw(ArgumentError("TikZ path export expects two or three axes."))
    lows, highs = _path_axis_bbox(all_boxes, axes_tuple)
    open(filename, "w") do io
        _write_tikz_preamble!(io, standalone, border)
        if length(axes_tuple) == 2
            _write_tikz_2d!(
                io,
                boxes,
                axes_tuple,
                lows,
                highs;
                width=width,
                height=height,
                color=color,
                trace_points=trace_points,
                show_trace=show_trace,
                trace_color=trace_color,
            )
        else
            _write_tikz_3d!(
                io,
                boxes,
                axes_tuple,
                lows,
                highs;
                width=width,
                height=height,
                depth=depth,
                color=color,
                trace_points=trace_points,
                show_trace=show_trace,
                trace_color=trace_color,
            )
        end
        standalone && println(io, "\\end{document}")
    end
    return filename
end

function _path_obj_point(box::PathBox, axes, lows, highs)
    coords = Float64[]
    for i in 1:3
        lo, hi = _path_axis_bounds(box, axes[i])
        push!(coords, _path_scale((lo + hi) / 2, lows[i], highs[i], 1.0))
    end
    return coords
end

function _path_obj_corners(box::PathBox, axes, lows, highs)
    bounds = [_path_axis_bounds(box, axis) for axis in axes]
    corners = Vector{Vector{Float64}}()
    for mask in 0:7
        point = Float64[]
        for j in 1:3
            lo, hi = bounds[j]
            v = ((mask >> (j - 1)) & 1) == 1 ? hi : lo
            push!(point, _path_scale(v, lows[j], highs[j], 1.0))
        end
        push!(corners, point)
    end
    return corners
end

function export_path_obj(input, filename::AbstractString; axes=(:t, 1, (1, :imag)))
    boxes = _path_export_boxes(_as_path_boxes(input))
    isempty(boxes) && throw(ArgumentError("no path boxes available to visualize. Run with store_boxes=:full or visualize=true."))
    axes_tuple = _normalize_path_axes(axes)
    length(axes_tuple) == 3 || throw(ArgumentError("OBJ path export expects exactly three axes."))
    lows, highs = _path_axis_bbox(boxes, axes_tuple)
    open(filename, "w") do io
        println(io, "# CertifiedHomotopyTracking.jl path boxes")
        vertex_offset = 0
        faces = ((1, 2, 4, 3), (5, 7, 8, 6), (1, 5, 6, 2), (3, 4, 8, 7), (1, 3, 7, 5), (2, 6, 8, 4))
        for (i, box) in enumerate(boxes)
            println(io, "o path_box_$i")
            for corner in _path_obj_corners(box, axes_tuple, lows, highs)
                println(io, "v ", join(corner, " "))
            end
            for face in faces
                println(io, "f ", join((vertex_offset + j for j in face), " "))
            end
            vertex_offset += 8
        end
    end
    return filename
end

_visualization_option(options, key::Symbol, default=nothing) =
    options isa NamedTuple && haskey(options, key) ? getproperty(options, key) : default

function _visualization_output_path(visualize, filename, options)
    filename !== nothing && return String(filename)
    opt_filename = _visualization_option(options, :filename)
    opt_filename !== nothing && return String(opt_filename)
    visualize isa AbstractString && return String(visualize)
    visualize isa NamedTuple && haskey(visualize, :filename) && return String(visualize.filename)
    return nothing
end

function _visualization_axes(visualize, axes, options)
    axes !== nothing && return axes
    opt_axes = _visualization_option(options, :axes)
    opt_axes !== nothing && return opt_axes
    visualize isa NamedTuple && haskey(visualize, :axes) && return visualize.axes
    return nothing
end

function _visualization_format(path)
    path === nothing && return :none
    lower = lowercase(path)
    endswith(lower, ".obj") && return :obj
    return :tikz
end

function _visualization_show_trace(visualize, show_trace, options)
    show_trace !== nothing && return show_trace
    opt_show_trace = _visualization_option(options, :show_trace)
    opt_show_trace !== nothing && return Bool(opt_show_trace)
    visualize isa NamedTuple && haskey(visualize, :show_trace) && return Bool(visualize.show_trace)
    return false
end

function _maybe_export_path_visualization(
    viz::PathVisualization,
    visualize,
    filename;
    axes=nothing,
    show_trace=nothing,
    visualize_options=(;),
)
    path = _visualization_output_path(visualize, filename, visualize_options)
    path === nothing && return viz
    export_axes = _visualization_axes(visualize, axes, visualize_options)
    if _visualization_format(path) === :obj
        export_path_obj(viz, path; axes = export_axes === nothing ? (:t, 1, (1, :imag)) : export_axes)
    else
        export_path_tikz(
            viz,
            path;
            axes = export_axes,
            width = _visualization_option(visualize_options, :width, 12.0),
            height = _visualization_option(visualize_options, :height, 7.0),
            depth = _visualization_option(visualize_options, :depth, 6.0),
            color = _visualization_option(visualize_options, :color, "blue"),
            standalone = _visualization_option(visualize_options, :standalone, true),
            border = _visualization_option(visualize_options, :border, "4pt"),
            show_trace = _visualization_show_trace(visualize, show_trace, visualize_options),
            trace_color = _visualization_option(visualize_options, :trace_color, "black"),
        )
    end
    return viz
end
