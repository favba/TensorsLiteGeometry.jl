module TensorsLiteGeometry
using TensorsLite, ImmutableVectors, LinearAlgebra

export circumcenter, closest, possible_positions_periodic, centroid, mass_centroid
export area, is_obtuse, in_triangle, in_polygon
export circle_edge_intersection, polygon_circle_intersection_area

export arc_length, arc_midpoint, spherical_triangle_angles, spherical_polygon_area, in_spherical_triangle
export lonlat_to_position, in_spherical_polygon
export spherical_polygon_moment, spherical_polygon_centroid

export periodic_to_base_point, isapprox_periodic
export eastward_vector, northward_vector

const VecOrTuple{T} = Union{<:(NTuple{N, T} where {N}), <:AbstractVector{T}}

@inline precise_norm(a) = norm(a)

@inline precise_norm(a::Vec2Dxy) = hypot(a.x, a.y)
@inline precise_norm(a::Vec2Dxz) = hypot(a.x, a.z)
@inline precise_norm(a::Vec2Dyz) = hypot(a.y, a.z)

@inline precise_norm(a::Vec3D) = hypot(a.x, a.y, a.z)

#Robust way to compute angles
#atan(norm(a × b), (a ⋅ b))
@inline angle(a, b) = atan(precise_norm(a × b), (a ⋅ b))

"""
    circumcenter(a::Vec,b::Vec,c::Vec) -> Vec

Returns the circumcenter position of the triangle formed by points `a`,`b`,`c`
"""
@inline function circumcenter(a::Vec, b::Vec, c::Vec)

    ab = b - a
    ac = c - a
    A = angle(ab, ac)
    bc = c - b
    B = angle(-ab, bc)
    C = π - A - B

    sin2A = sin(2 * A)
    sin2B = sin(2 * B)
    sin2C = sin(2 * C)

    return (muladd(sin2A, a, muladd(sin2B, b, sin2C * c))) / (sin2A + sin2B + sin2C)
end

@generated function closest(p::Vec, points::Tuple{Vararg{T, N}}) where {T <: Vec, N}
    quote
        $(Expr(:meta, :inline))
        p_i = points[1]
        min_dist = norm(p - points[1])
        min_p = p_i

        Base.Cartesian.@nexprs $(N - 1) i -> begin
            p_i = points[i + 1]
            norm_dist = norm(p - p_i)
            if norm_dist < min_dist
                min_p = p_i
                min_dist = norm_dist
            end
        end

        return min_p
    end
end

@inline closest(points::Tuple, p::Vec) = closest(p, points)

@inline function possible_positions_periodic(p::Vec, xp::Number, yp::Number)
    xpi = xp * 𝐢
    ypj = yp * 𝐣
    p_p_x = p + xpi
    p_m_x = p - xpi
    p_p_y = p + ypj
    p_m_y = p - ypj
    p_p_x_p_y = p_p_x + ypj
    p_m_x_p_y = p_m_x + ypj
    p_p_x_m_y = p_p_x - ypj
    p_m_x_m_y = p_m_x - ypj

    return (p, p_p_x, p_m_x, p_p_y, p_m_y, p_p_x_p_y, p_m_x_p_y, p_p_x_m_y, p_m_x_m_y)
end

#This only work if distance between points is much smaller than periods
#@inline function closest(p::Vec2Dxy, p2::Vec2Dxy, xp::Number, yp::Number)
#    #d = min(xp, yp) / 2
#    #norm(p2 - p) < d && return p2
#    #return closest(p, possible_positions_periodic(p2, xp, yp)[2:9])
#    px = p.x
#    py = p.y
#    p2x = p2.x
#    p2y = p2.y
#    dx = p2x - px
#    dy = p2y - py
#
#    adx = abs(dx)
#    ady = abs(dy)
#
#    pfx = if (adx < 0.5*xp)
#        p2x
#    else
#        p2xp = p2x + xp
#        if (abs(p2xp - px) < adx)
#            p2xp 
#        else 
#            p2x - xp
#        end
#    end
#
#    pfy = if (ady < 0.5*yp)
#        p2y
#    else
#        p2yp = p2y + yp
#        if (abs(p2yp - py) < ady)
#            p2yp 
#        else 
#            p2y - yp
#        end
#    end
#
#    return Vec(x=pfx, y=pfy)
#end

@inline function closest(base_val::Real, val::T, period=T(360)) where {T<:Real}
    vals1 = val - period
    vals2 = val
    vals3 = val + period

    min_diff = abs(vals1 - base_val)
    result = vals1

    diff = abs(vals2 - base_val)
    if diff < min_diff
        result = vals2
        min_diff = diff
    end

    if abs(vals3 - base_val) < min_diff
        result = vals3
    end
    return result
end

@inline closest(p::Vec2Dxy, p2::Vec2Dxy, xp::Number, yp::Number) = Vec(x = closest(p.x, p2.x, xp), y = closest(p.y, p2.y, yp))

@inline area(a::Vec, b::Vec, c::Vec) = 0.5 * norm((b - a) × (c - b))

for N in 4:12
    @eval function area(vars::Vararg{Vec, $N})
        $(Expr(:meta, :inline))
        area(ntuple(i -> vars[i], $(Val(N - 1)))...) + area(vars[1], vars[$(N - 1)], vars[$N])
    end
end

@inline function in_triangle(p::Vec, at::Number, a::Vec, b::Vec, c::Vec)

    aa = area(p, b, c)
    ab = area(p, c, a)
    ac = area(p, a, b)

    inside = (aa + ab + ac) ≈ at

    return inside, aa, ab, ac, at
end

""" Checks if point `p` is inside the triangle given by the vertices `a`,`b`,`c`
"""
@inline in_triangle(p::Vec, a::Vec, b::Vec, c::Vec) = in_triangle(p, area(a, b, c), a, b, c)

function area(points::Union{<:AbstractVector{T}, NTuple{N, T}}) where {T <: AbstractVec, N}
    @inbounds p1 = points[1]
    @inbounds p2 = points[2]
    a = zero(nonzero_eltype(T))

    @inbounds for i in Iterators.drop(eachindex(points), 2)
        p3 = points[i]
        a += area(p1, p2, p3)
        p2 = p3
    end

    return a
end

function area(points::AbstractVector{T}, indices::VecOrTuple) where {T <: AbstractVec}
    @inbounds p1 = points[indices[1]]
    @inbounds p2 = points[indices[2]]
    a = zero(nonzero_eltype(T))

    @inbounds for i in Iterators.drop(indices, 2)
        p3 = points[i]
        a += area(p1, p2, p3)
        p2 = p3
    end

    return a
end

function area(points::AbstractVector{T}, indices::VecOrTuple, x_period::Number, y_period::Number) where {T <: AbstractVec}
    @inbounds p1 = points[indices[1]]
    @inbounds p2 = closest(p1, points[indices[2]], x_period, y_period)
    a = zero(nonzero_eltype(T))

    @inbounds for i in Iterators.drop(indices, 2)
        p3 = closest(p1, points[i], x_period, y_period)
        a += area(p1, p2, p3)
        p2 = p3
    end

    return a
end

"""
    centroid(a::Vec,b::Vec,c::Vec) -> Vec

Returns the centroid (mass center) position of the triangle formed by points `a`,`b`,`c`
"""
@inline centroid(a::Vec, b::Vec, c::Vec) = (a + b + c) / 3

"""
    mass_centroid(ρ::Function, a::Vec,b::Vec,c::Vec) -> Vec

Returns the centroid (mass center) position of the triangle formed by points `a`,`b`,`c` and density funciton `ρ(𝐱)`
The result is an approximation.
"""
@inline function mass_centroid(ρ::F, a::Vec, b::Vec, c::Vec) where {F <: Function}
    ρa = ρ(a)
    ρb = ρ(b)
    ρc = ρ(c)
    ρm = (ρa + ρb + ρc) / 3
    return ((1 + ((2ρa) / ρm)) * a + (1 + ((2ρb) / ρm)) * b + (1 + ((2ρc) / ρm)) * c) / 9
end

"""
    centroid(points) -> Vec

Returns the centroid (mass center) position of the polygon formed by `points`.
"""
@inline function centroid(points::AbstractVector{T}) where {T <: AbstractVec}

    @inbounds p1 = points[1]
    @inbounds p2 = points[2]
    @inbounds p3 = points[3]

    ta = area(p1, p2, p3)
    c = centroid(p1, p2, p3)
    p2 = p3

    @inbounds for i in Iterators.drop(eachindex(points), 3)
        p3 = points[i]
        a = area(p1, p2, p3)
        ta += a
        w = a / ta
        c = (1 - w) * c + w * centroid(p1, p2, p3)
        p2 = p3
    end

    return c
end

"""
    mass_centroid(ρ, points) -> Vec

Returns the centroid (mass center) position of the polygon formed by `points` with density function `ρ(𝐱)`.
The result is an approximation due to the assumption that `ρ` is linear inside each triangle that forms the polygon.
"""
@inline function mass_centroid(ρ::F, points::AbstractVector{T}) where {F <: Function, T <: AbstractVec}

    @inbounds p1 = points[1]
    @inbounds p2 = points[2]
    @inbounds p3 = points[3]

    rp1 = ρ(p1)
    rp2 = ρ(p2)
    rp3 = ρ(p2)
    rm = (rp1 + rp2 + rp3) / 3

    tmass = rm * area(p1, p2, p3)
    c = ((1 + ((2rp1) / rm)) * p1 + (1 + ((2rp2) / rm)) * p2 + (1 + ((2rp3) / rm)) * p3) / 9

    p2 = p3
    rp2 = rp3

    @inbounds for i in Iterators.drop(eachindex(points), 3)
        p3 = points[i]
        rp3 = ρ(p3)
        rm = (rp1 + rp2 + rp3) / 3
        mass = rm * area(p1, p2, p3)
        tmass += mass
        w = mass / tmass
        c = (1 - w) * c + w * (((1 + ((2rp1) / rm)) * p1 + (1 + ((2rp2) / rm)) * p2 + (1 + ((2rp3) / rm)) * p3) / 9)
        p2 = p3
        rp2 = rp3
    end

    return c
end


"""
    centroid(points, indices) -> Vec

Returns the centroid (mass center) position of the polygon formed by points `getindex.((points,), indices)`
"""
@inline function centroid(points::AbstractVector{T}, indices::VecOrTuple) where {T <: AbstractVec}

    @inbounds p1 = points[indices[1]]
    @inbounds p2 = points[indices[2]]
    @inbounds p3 = points[indices[3]]

    ta = area(p1, p2, p3)
    c = centroid(p1, p2, p3)
    p2 = p3

    @inbounds for i in Iterators.drop(indices, 3)
        p3 = points[i]
        a = area(p1, p2, p3)
        ta += a
        w = a / ta
        c = (1 - w) * c + w * centroid(p1, p2, p3)
        p2 = p3
    end

    return c
end

"""
    mass_centroid(ρ, points, indices) -> Vec

Returns the centroid (mass center) position of the polygon formed by points `getindex.((points,), indices)` with density funciton `ρ(𝐱)`.
The result is an approximation due to the assumption that `ρ` is linear inside each triangle that forms the polygon.
"""
@inline function mass_centroid(ρ::F, points::AbstractVector{T}, indices::VecOrTuple) where {F <: Function, T <: AbstractVec}

    @inbounds p1 = points[indices[1]]
    @inbounds p2 = points[indices[2]]
    @inbounds p3 = points[indices[3]]

    rp1 = ρ(p1)
    rp2 = ρ(p2)
    rp3 = ρ(p2)
    tr = rp1 + rp2 + rp3

    tmass = (tr / 3) * area(p1, p2, p3)
    c = (rp1 / tr) * p1 + (rp2 / tr) * p2 + (rp3 / tr) * p3

    p2 = p3
    rp2 = rp3

    @inbounds for i in Iterators.drop(indices, 3)
        p3 = points[i]
        rp3 = ρ(p3)
        tr = rp1 + rp2 + rp3
        mass = (tr / 3) * area(p1, p2, p3)
        tmass += mass
        w = mass / tmass
        c = (1 - w) * c + w * ((rp1 / tr) * p1 + (rp2 / tr) * p2 + (rp3 / tr) * p3)
        p2 = p3
        rp2 = rp3
    end

    return c
end

@inline function centroid(points::AbstractVector{T}, indices::VecOrTuple, x_period::Number, y_period::Number) where {T <: AbstractVec}

    @inbounds p1 = points[indices[1]]
    @inbounds p2 = closest(p1, points[indices[2]], x_period, y_period)
    @inbounds p3 = closest(p1, points[indices[3]], x_period, y_period)

    ta = area(p1, p2, p3)
    c = centroid(p1, p2, p3)
    p2 = p3

    @inbounds for i in Iterators.drop(indices, 3)
        p3 = closest(p1, points[i], x_period, y_period)
        a = area(p1, p2, p3)
        ta += a
        w = a / ta
        c = (1 - w) * c + w * centroid(p1, p2, p3)
        p2 = p3
    end

    return c
end


"""
    in_polygon(point::Vec, vertices::AbstractVector{<:Vec}, indices) -> Bool

Whether `point` is inside a polygon with vertices given by `getindex.((vertices,),indices)`
"""
function in_polygon(p::Vec, points::AbstractVector{T}, indices::VecOrTuple) where {T <: AbstractVec}
    @inbounds p1 = points[indices[1]]
    @inbounds p2 = points[indices[2]]

    inside = false
    @inbounds for i in Iterators.drop(indices, 2)
        p3 = points[i]
        inside = in_triangle(p, p1, p2, p3)[1]
        inside && break
        p2 = p3
    end

    return inside
end

"""
    in_polygon_periodic(point::Vec, vertices::AbstractVector{<:Vec}, indices,x_period::Number,y_period::Number) -> Bool

Whether `point` is inside a polygon with vertices given by `getindex.((vertices,),indices) adjusted to the fact that some vertices positions might overflow the periodic domain.`
"""
function in_polygon(p::Vec, points::AbstractVector{T}, indices::VecOrTuple, x_period::Number, y_period::Number) where {T <: AbstractVec}
    @inbounds p1 = points[indices[1]]
    @inbounds p2 = closest(p1, points[indices[2]], x_period, y_period)

    inside = false
    @inbounds for i in Iterators.drop(indices, 2)
        p3 = closest(p1, points[i], x_period, y_period)
        inside = in_triangle(p, p1, p2, p3)[1]
        inside && break
        p2 = p3
    end

    return inside
end

"""
    is_obtuse(a::Vec,b::Vec,c::Vec) -> Bool

Whether the triangle formed by points `a`,`b`,`c` is obtuse
"""
function is_obtuse(a::Vec, b::Vec, c::Vec)
    ab = b - a
    ac = c - a

    #cos(0) <= 0
    A = angle(ac, ab)
    A >= π / 2 && return true

    bc = c - b
    C = angle(bc, ac)

    r = if ((C >= π / 2) || (C + A <= π / 2))
        true
    else
        false
    end
    return r
end

"""
    circle_edge_intersection(vertex_1::Vec2D, vertex_2::Vec2D, center::Vec2D, r2::Number) -> Vec2D

Given edge vertices `vertex_1` and `vertex_2` returns the intersection point of the edge with a circle with radius `sqrt(r2)` centered at `center` .
"""
function circle_edge_intersection(p1, p2, center, r2::Number)
    term1 = p1 - p2
    term2 = p2 - center
    a = term1 ⋅ term1
    b = 2 * term1 ⋅ term2
    c = term2 ⋅ term2 - r2
    sqrt_bterm = sqrt(b * b - 4 * a * c)

    α1, α2 = ((-b + sqrt_bterm) / 2a, (-b - sqrt_bterm) / 2a)

    if ((0.0 <= α1) && (α1 <= 1.0))
        α = α1
    elseif ((0.0 <= α2) && (α2 <= 1.0))
        α = α2
    else
        error("Circle doesn't cross edge")
    end

    return α * p1 + (1 - α) * p2
end

function is_rev_sorted(in_circle::AbstractVector)
    r = true
    for i in Iterators.drop(eachindex(in_circle), 1)
        r = in_circle[i - 1] >= in_circle[i]
        r || break
    end
    return r
end

function is_in_circle(c::Vec, r2::Number, point::Vec)
    cp = point - c
    return cp ⋅ cp <= r2
end

function polygon_circle_intersection_area(center, r2, vertices_positions::ImmutableVector{N}, in_circle::ImmutableVector{N, Bool}) where {N}
    l_in_c = in_circle
    lv = vertices_positions

    # Rearange terms to get all (and only) terms inside the circle placed at the beginning of the ImmutableVector
    while !is_rev_sorted(l_in_c)
        l_in_c = circshift(l_in_c, 1)
        lv = circshift(lv, 1)
    end

    nv_in_circle = findfirst(isequal(false), l_in_c) - 1

    first_intersection = @inbounds circle_edge_intersection(lv[nv_in_circle], lv[nv_in_circle + 1], center, r2)
    second_intersection = @inbounds circle_edge_intersection(lv[end], lv[1], center, r2)

    new_polygon = @inbounds push(ImmutableVector{N + 1}(lv[1:nv_in_circle]), first_intersection, second_intersection)

    area1 = area(new_polygon)

    cp1 = first_intersection - center
    cp2 = second_intersection - center

    #Arc angle
    θ = angle(cp1, cp2)

    area2 = θ * r2 / 2 - area(center, first_intersection, second_intersection)

    return area1 + area2
end

function polygon_circle_intersection_area(center, r2, vertices_positions::ImmutableVector{N}) where {N}
    in_circle = map(v -> (is_in_circle(center, r2, v)), vertices_positions)
    return polygon_circle_intersection_area(center, r2, vertices_positions, in_circle)
end

"""
    periodic_to_base_point(p::Vec2Dxy, lx, ly) -> base_point::Vec2Dxy

Given a point `p` returns the equivalent point that lies in `[0, lx) × [0, ly)` (assuming the domain is bi-periodic with periods `lx` and `ly`).
"""
function periodic_to_base_point(p::Vec2Dxy, lx::Number, ly::Number)
    r = p
    while r.x >= lx
        r -= lx*𝐢
    end
    while r.x < 0.0
        r += lx*𝐢
    end
    while r.y >= ly
        r -= ly*𝐣
    end
    while r.y < 0.0
        r += ly*𝐣
    end
    return r
end

"""
    isapprox_periodic(a::Vec2Dxy, b::Vec2Dxy, lx, ly) -> Bool

Check if `a` and `b` represent the same base point (to [`isapprox`](@ref) tolerance) in a periodic box with `lx` and `ly` periods.
"""
function isapprox_periodic(p1::Vec2Dxy{T}, p2::Vec2Dxy, lx::Number, ly::Number; atol = 0, rtol = sqrt(eps(T))) where {T} 
    bp1 = periodic_to_base_point(p1, lx, ly)
    bp2 = periodic_to_base_point(p2, lx, ly)
    # using `closest` because we can't compare bp1 and bp2 directly when they fall "exactly" at the periodic boxes edges
    return isapprox(bp1, closest(bp1, bp2, lx, ly), atol = atol, rtol = rtol)
end

include("spherical_geometry.jl")

end
