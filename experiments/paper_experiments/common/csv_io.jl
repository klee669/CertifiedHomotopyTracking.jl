function write_csv(path, rows)
    isempty(rows) && return path
    mkpath(dirname(path))
    header = collect(keys(first(rows)))
    open(path, "w") do io
        println(io, join(header, ','))
        for row in rows
            println(io, join((getfield(row, name) for name in header), ','))
        end
    end
    return path
end
