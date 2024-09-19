@inline spherical_distance(R::Number, a::Vec, b::Vec) = R*atan(norm(a × b) / (a ⋅ b))
@inline spherical_distance(a::Vec, b::Vec) = spherical_distance(norm(a), a, b)

@inline spherical_midpoint(R::Number, a::Vec, b::Vec) = R * normalize((a + b) / 2)
@inline spherical_midpoint(a::Vec, b::Vec) = spherical_midpoint(norm(a), a, b)

@inline function spherical_angles(R::Number, a::Vec, b::Vec, c::Vec)
    â = a / R
    b̂ = b / R
    ĉ = c / R

    n_ab = â × b̂
    n_ac = â × ĉ

    A = acos(n_ab ⋅ n_ac)

    n_bc = b̂ × ĉ

    B = acos(-(n_bc ⋅ n_ab))

    C = acos(n_ac ⋅ n_bc)

    return Vec(A, B, C)
end

@inline spherical_angles(a::Vec, b::Vec, c::Vec) = spherical_angles(norm(a), a, b, c)

@inline spherical_area(R::Number, a::Vec, b::Vec, c::Vec) = (R * R) * (sum(spherical_angles(R, a, b, c)) - π)

for N in 4:12
    @eval function spherical_area(R::Number, vars::Vararg{Vec,$N})
        $(Expr(:meta,:inline))
        spherical_area(R, ntuple(i -> vars[i], $(Val(N - 1)))...) + spherical_area(R, vars[1], vars[$(N - 1)], vars[$N])
    end
end

@inline spherical_area(v::Vararg{Vec,N}) where {N} = spherical_area(norm(v[1]), v...)

@inline function spherical_area(R::Number, points::Union{<:AbstractVector{T}, NTuple{N, T}}) where {T <: AbstractVec,N}
    @inbounds p1 = points[1]
    @inbounds p2 = points[2]
    a = zero(nonzero_eltype(T))

    @inbounds for i in Iterators.drop(eachindex(points),2)
        p3 = points[i]
        a += spherical_area(R, p1, p2, p3)
        p2 = p3
    end

    return a
end

@inline spherical_area(points::Union{<:AbstractVector{T}, NTuple{N, T}}) where {T <: AbstractVec,N} = spherical_area(norm(@inbounds(points[1])), points)

@inline function spherical_area(R::Number, points::AbstractVector{T}, indices::VecOrTuple) where {T<:AbstractVec}
    @inbounds p1 = points[indices[1]]
    @inbounds p2 = points[indices[2]]
    a = zero(nonzero_eltype(T))

    @inbounds for i in Iterators.drop(indices, 2)
        p3 = points[i]
        a += spherical_area(R, p1, p2, p3)
        p2 = p3
    end

    return a
end

@inline spherical_area(points::AbstractVector{T}, indices::VecOrTuple) where {T<:AbstractVec} = spherical_area(norm(@inbounds(points[indices[1]])), points, indices)

@inline function lonlat_to_position(R::Number, lon::Number, lat::Number)
    sinlon = sin(lon)
    rsinlat = R * sin(lat)
    coslon = cos(lon)
    rcoslat = R * cos(lat)
    return Vec(rcoslat * coslon, rcoslat * sinlon, rsinlat)
end

@inline function in_spherical_triangle(R::T, p::Vec, a::Vec, b::Vec, c::Vec) where {T <: Number}
    ep = eps(T)
    p̂ = p/R
    â = a/R
    b̂ = b/R
    ĉ = c/R
    if signbit(p̂ ⋅ (ĉ × â) + ep) || signbit(p̂ ⋅ (â × b̂) + ep) || signbit(p̂ ⋅ (b̂ × ĉ) + ep) # add eps so we don't get false negative when p is equal to either a, b, or c
        return false
    else 
        return true
    end
end

@inline in_spherical_triangle(p::Vec, a::Vec, b::Vec, c::Vec) = in_spherical_triangle(norm(a), p, a, b, c)

function in_spherical_polygon(R::Number, p::Vec,points::AbstractVector{T},indices::VecOrTuple) where {T<:AbstractVec}
    @inbounds p1 = points[indices[1]]
    @inbounds p2 = points[indices[2]]

    inside = false
    @inbounds for i in Iterators.drop(indices,2)
        p3 = points[i]
        inside = in_spherical_triangle(R,p,p1,p2,p3)
        inside && break
        p2 = p3
    end

    return inside
end

"""
    in_spherical_polygon(point::Vec, vertices::AbstractVector{<:Vec}, indices) -> Bool

Whether `point` is inside a given spherical polygon with vertices given by `getindex.((vertices,),indices)`
"""
@inline in_spherical_polygon(p::Vec, points::AbstractVector{T}, indices::VecOrTuple) where {T <: AbstractVec} = in_spherical_polygon(norm(points[indices[1]]), p, points, indices)
