@inline arc_length(R::Number, a::Vec, b::Vec) = R * angle(a, b)
@inline arc_length(a::Vec, b::Vec) = arc_length(norm(a), a, b)

@inline arc_midpoint(R::Number, a::Vec, b::Vec) = R * normalize((a + b) / 2)
@inline arc_midpoint(a::Vec, b::Vec) = arc_midpoint(norm(a), a, b)

@inline function spherical_triangle_angles(a::Vec, b::Vec, c::Vec)

    n_ab = a × b
    n_ac = a × c

    A = angle(n_ab, n_ac)

    n_bc = b × c

    B = angle(-n_bc, n_ab)

    C = angle(n_ac, n_bc)

    return (A, B, C)
end


@inline spherical_polygon_area(R::Number, a::Vec, b::Vec, c::Vec) = (R * R) * (sum(spherical_triangle_angles(a, b, c)) - π)

for N in 4:12
    @eval function spherical_polygon_area(R::Number, vars::Vararg{Vec, $N})
        $(Expr(:meta, :inline))
        spherical_polygon_area(R, ntuple(i -> vars[i], $(Val(N - 1)))...) + spherical_polygon_area(R, vars[1], vars[$(N - 1)], vars[$N])
    end
end

@inline spherical_polygon_area(v::Vararg{Vec, N}) where {N} = spherical_polygon_area(norm(v[1]), v...)

@inline function spherical_polygon_area(R::Number, points::Union{<:AbstractVector{T}, NTuple{N, T}}) where {T <: AbstractVec, N}
    @inbounds p1 = points[1]
    @inbounds p2 = points[2]
    a = zero(nonzero_eltype(T))

    @inbounds for i in Iterators.drop(eachindex(points), 2)
        p3 = points[i]
        a += spherical_polygon_area(R, p1, p2, p3)
        p2 = p3
    end

    return a
end

@inline spherical_polygon_area(points::Union{<:AbstractVector{T}, NTuple{N, T}}) where {T <: AbstractVec, N} = spherical_polygon_area(norm(@inbounds(points[1])), points)

@inline function spherical_polygon_area(R::Number, points::AbstractVector{T}, indices::VecOrTuple) where {T <: AbstractVec}
    @inbounds p1 = points[indices[1]]
    @inbounds p2 = points[indices[2]]
    a = zero(nonzero_eltype(T))

    @inbounds for i in Iterators.drop(indices, 2)
        p3 = points[i]
        a += spherical_polygon_area(R, p1, p2, p3)
        p2 = p3
    end

    return a
end

@inline spherical_polygon_area(points::AbstractVector{T}, indices::VecOrTuple) where {T <: AbstractVec} = spherical_polygon_area(norm(@inbounds(points[indices[1]])), points, indices)

@inline function lonlat_to_position(R::Number, lon::Number, lat::Number)
    sinlon = sin(lon)
    rsinlat = R * sin(lat)
    coslon = cos(lon)
    rcoslat = R * cos(lat)
    return Vec(rcoslat * coslon, rcoslat * sinlon, rsinlat)
end

@inline function in_spherical_triangle(R::T, p::Vec, a::Vec, b::Vec, c::Vec) where {T <: Number}
    ep = eps(T)
    p̂ = p / R
    â = a / R
    b̂ = b / R
    ĉ = c / R
    if signbit(p̂ ⋅ (ĉ × â) + ep) || signbit(p̂ ⋅ (â × b̂) + ep) || signbit(p̂ ⋅ (b̂ × ĉ) + ep) # add eps so we don't get false negative when p is equal to either a, b, or c
        return false
    else
        return true
    end
end

@inline in_spherical_triangle(p::Vec, a::Vec, b::Vec, c::Vec) = in_spherical_triangle(norm(a), p, a, b, c)

function in_spherical_polygon(R::Number, p::Vec, points::AbstractVector{T}, indices::VecOrTuple) where {T <: AbstractVec}
    @inbounds p1 = points[indices[1]]
    @inbounds p2 = points[indices[2]]

    inside = false
    @inbounds for i in Iterators.drop(indices, 2)
        p3 = points[i]
        inside = in_spherical_triangle(R, p, p1, p2, p3)
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

@inline function spherical_polygon_moment(R::Number, a::Vec, b::Vec, c::Vec)
    ab = normalize(a × b) * angle(a, b)
    bc = normalize(b × c) * angle(b, c)
    ca = normalize(c × a) * angle(c, a)
    return R * (ab + bc + ca) / 2
end

@inline function spherical_polygon_moment(R::Number, points::VecOrTuple{T}) where {T <: AbstractVec}

    @inbounds p_1 = points[1]
    @inbounds p2 = points[2]

    a = normalize(p_1 × p2) * angle(p_1, p2)

    @inbounds for i in Iterators.drop(eachindex(points), 2)
        p1 = p2
        p2 = points[i]
        a += normalize(p1 × p2) * angle(p1, p2)
    end

    a += normalize(p2 × p_1) * angle(p2, p_1)

    return R * a / 2
end

@inline function spherical_polygon_moment(R::Number, points::AbstractVector{T}, indices::VecOrTuple) where {T <: AbstractVec}

    @inbounds p_1 = points[indices[1]]
    @inbounds p2 = points[indices[2]]

    a = normalize(p_1 × p2) * angle(p_1, p2)

    @inbounds for i in Iterators.drop(indices, 2)
        p1 = p2
        p2 = points[i]
        a += normalize(p1 × p2) * angle(p1, p2)
    end

    a += normalize(p2 × p_1) * angle(p2, p_1)

    return R * a / 2
end

@inline spherical_polygon_centroid(R::Number, a::Vec, b::Vec, c::Vec) = R * normalize(spherical_polygon_moment(R, a, b, c))
@inline spherical_polygon_centroid(a::Vec, b::Vec, c::Vec) = spherical_polygon_centroid(norm(a), b, c)

@inline spherical_polygon_centroid(R::Number, points::VecOrTuple{T}) where {T <: AbstractVec} = R * normalize(spherical_polygon_moment(R, points))
@inline spherical_polygon_centroid(points::VecOrTuple{T}) where {T <: AbstractVec} = spherical_polygon_centroid(norm(points[1]), points)

@inline spherical_polygon_centroid(R::Number, points::AbstractVector{T}, indices::VecOrTuple) where {T <: AbstractVec} = R * normalize(spherical_polygon_moment(R, points, indices))
@inline spherical_polygon_centroid(points::AbstractVector{T}, indices::VecOrTuple) where {T <: AbstractVec} = spherical_polygon_centroid(norm(points[indices[1]]), points, indices)

"""
    eastward_vector(position::Vec) -> east_unit_vector::Vec2D

Return the eastward unit vector to the `position` vector.
"""
@inline function eastward_vector(p::Vec3D)
    ϕ = atan(p.y, p.x)
    return cos(ϕ) * 𝐣 - sin(ϕ) * 𝐢
end

"""
    northward_vector(position::Vec) -> northward_unit_vector::Vec3D

Return the northward unit vector to the `position` vector.
"""
@inline function northward_vector(p::Vec3D)
    x = p.x
    y = p.y
    z = p.z
    ϕ = atan(y, x)
    θ = atan(sqrt(x*x + y*y), z)
    sinϕ = sin(ϕ)
    cosϕ = cos(ϕ)
    sinθ = sin(θ)
    cosθ = cos(θ)
    return sinθ*𝐤 - (cosϕ*cosθ)*𝐢 - (sinϕ*cosθ)*𝐣
end

