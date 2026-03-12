export toCCi,
    max_int_norm, 
    max_norm, 
    midpoint_complex_int, 
    midpoint_complex_box, 
    convert_to_box_int, 
    convert_to_box_matrix, 
    convert_to_double_int, 
    convert_to_double_matrix,
    mag_complex,
    norm_inf,
    to_acb_vec,
    get_mid,
    get_mid_vec

using AbstractAlgebra


function toCCi(x, F::AbstractAlgebra.Field)
    if x isa Number
        return F(x)
    elseif x isa Complex
        return F(real(x), imag(x))
    elseif x isa AbstractVector
        return [toCCi(e, F) for e in x]
    elseif x isa AbstractMatrix
        return [toCCi(e, F) for e in x] |> x -> reshape(x, size(x))
    else
        error("Cannot convert object of type $(typeof(x)) to field element")
    end
end

#### max norm parts
function max_int_norm(interval::AcbFieldElem)
    r = real(interval);
    i = imag(interval);
    rmax = convert(Float64, abs(Nemo.midpoint(r)) + Nemo.radius(r));
    imax = convert(Float64, abs(Nemo.midpoint(i)) + Nemo.radius(i));
    maximum([rmax, imax])
end

function max_norm(intvec::Union{Matrix{AcbFieldElem},AcbMatrix})
    sz = size(intvec);
    nr = sz[1];
    nc = sz[2];
    abs_list = [];
    for i in 1:nr
        for j in 1:nc
            push!(abs_list,max_int_norm(intvec[i,j]));
        end
    end 
    maximum(abs_list)
end



#### midpoint functions for interval and interval vectors
function midpoint_complex_int(int::AcbFieldElem)
    ring = parent(int);
    ring(midpoint(real(int)),midpoint(imag(int)));
end

function midpoint_complex_box(int::Union{Vector{AcbFieldElem},AcbMatrix})
    n = length(int);
    ring = parent(int[1]);
    result = zeros(ring, n);
    for i in 1:n
        result[i]=midpoint_complex_int(int[i]);
    end
    result
end





function convert_to_double_int(interval::AcbFieldElem)
    Base.convert(Float64,real(interval)) + Base.convert(Float64,imag(interval))*im
end

function convert_to_double_matrix(M)
    sz = size(M);
    nr = sz[1];
    nc = sz[2];
    result = zeros(Complex{Float64},nr,nc);
    for j in 1:nr
        for i in 1:nc
            result[j,i] = convert_to_double_int(M[j,i]);
        end
    end
    result
end



function convert_to_box_int(vec,ring)
    n = length(vec);
    result = [];
    for i in 1:n
        push!(result,ring(real(vec[i]),imag(vec[i])));
    end
    result
end


function convert_to_box_matrix(M,ring)
    sz = size(M);
    nr = sz[1];
    nc = sz[2];
    result = zero_matrix(ring,nr,nc);
    for j in 1:nr
        for i in 1:nc
            result[j,i] = ring(real(M[j,i]),imag(M[j,i]));
        end
    end
    Matrix(result)
end


function mag_complex(z::AcbFieldElem)
    val_interval = abs(z)
    return Float64(val_interval)
end
norm_inf(v::AbstractArray{AcbFieldElem}) = maximum(mag_complex.(v))

function to_acb_vec(x::AbstractVector)
    return [CC(xi) for xi in x] 
end

function get_mid(z::AcbFieldElem)
    re = midpoint(real(z)) 
    im = midpoint(imag(z))
    return CC(re, im) 
end
get_mid_vec(v::AbstractVector{AcbFieldElem}) = get_mid.(v)
