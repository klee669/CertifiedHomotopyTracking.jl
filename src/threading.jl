function _thread_count(threading::Bool, ntasks::Integer)
    ntasks >= 1 || throw(ArgumentError("ntasks must be positive."))
    threading || return 1
    return min(Int(ntasks), Base.Threads.nthreads())
end

_threading_enabled(threading::Bool, ntasks::Integer) =
    _thread_count(threading, ntasks) > 1

function _thread_chunks(items::AbstractVector, ntasks::Integer)
    n = length(items)
    n == 0 && return UnitRange{Int}[]
    k = min(Int(ntasks), n)
    return [(((i - 1) * n) ÷ k + 1):((i * n) ÷ k) for i in 1:k]
end
