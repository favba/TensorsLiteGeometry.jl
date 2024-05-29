module TensorsLiteGeometry
using TensorsLite

export circumcenter, closest, possible_positions_periodic
export area, is_obtuse, in_triangle, in_polygon
export circle_edge_intersection

const VecOrTuple = Union{<:Tuple,<:AbstractVector}

"""
    circumcenter(a::Vec,b::Vec,c::Vec) -> Vec

Returns the circumcenter position of the triangle formed by points `a`,`b`,`c`
"""
@inline function circumcenter(a::Vec,b::Vec,c::Vec)

    ab = b-a
    ac = c-a
    norm_ab = norm(ab)
    norm_ac = norm(ac)
    A = acos((ab⋅ac)/(norm_ab*norm_ac))
    bc = c-b
    norm_bc = norm(bc)
    B = acos(-(ab⋅bc)/(norm_ab*norm_bc))
    C = π - A - B

    sin2A = sin(2*A)
    sin2B = sin(2*B)
    sin2C = sin(2*C)

    return (muladd(sin2A, a, muladd(sin2B, b, sin2C*c)))/(sin2A+sin2B+sin2C)
end

@generated function closest(p::Vec,points::Tuple{Vararg{T,N}}) where {T<:Vec,N}
    quote
        $(Expr(:meta,:inline))
        p_i = points[1]
        min_dist = norm(p-points[1])
        min_p = p_i

        Base.Cartesian.@nexprs $(N-1) i -> begin
            p_i = points[i+1]
            norm_dist = norm(p-p_i)
            if norm_dist < min_dist
                min_p = p_i
                min_dist = norm_dist
            end
        end
    
        return min_p
    end
end

@inline closest(points::Tuple,p::Vec) = closest(p,points)

@inline function possible_positions_periodic(p::Vec,xp::Number,yp::Number)
    xpi = xp*𝐢
    ypj = yp*𝐣
    p_p_x = p + xpi
    p_m_x = p - xpi
    p_p_y = p + ypj
    p_m_y = p - ypj
    p_p_x_p_y = p_p_x + ypj
    p_m_x_p_y = p_m_x + ypj
    p_p_x_m_y = p_p_x - ypj
    p_m_x_m_y = p_m_x - ypj

    return (p,p_p_x,p_m_x,p_p_y,p_m_y,p_p_x_p_y,p_m_x_p_y,p_p_x_m_y,p_m_x_m_y)
end

#@inline closest(p::Vec,p2::Vec,xp::Number,yp::Number) = closest(p,possible_positions_periodic(p2,xp,yp))

#This only work if distance between points is much smaller than periods
@inline function closest(p::Vec,p2::Vec,xp::Number,yp::Number)
    d = min(xp,yp)/2
    norm(p2-p) < d && return p2
    return closest(p,possible_positions_periodic(p2,xp,yp)[2:9])
end

@inline area(a::Vec,b::Vec,c::Vec) = 0.5*norm((b-a)×(c-b))

for N in 4:12
    @eval function area(vars::Vararg{Vec,$N})
        $(Expr(:meta,:inline))
        area(ntuple(i->vars[i],$(Val(N-1)))...) + area(vars[1],vars[$(N-1)],vars[$N])
    end
end

@inline function in_triangle(p::Vec,at::Number,a::Vec,b::Vec,c::Vec)

    aa = area(p,b,c)
    ab = area(p,c,a)
    ac = area(p,a,b)

    inside = (aa+ab+ac) ≈ at

    return inside,aa,ab,ac,at
end

""" Checks if point `p` is inside the triangle given by the vertices `a`,`b`,`c`
"""
@inline in_triangle(p::Vec,a::Vec,b::Vec,c::Vec) = in_triangle(p,area(a,b,c),a,b,c)

function area(points::Union{<:AbstractVector{T},NTuple{N,T}}) where {T<:AbstractVec,N}
    @inbounds p1 = points[1]
    @inbounds p2 = points[2]
    a = zero(nonzero_eltype(T))

    @inbounds for i in Iterators.drop(eachindex(points),2)
        p3 = points[i]
        a += area(p1,p2,p3)
        p2 = p3
    end

    return a
end

function area(points::AbstractVector{T},indices::VecOrTuple) where {T<:AbstractVec}
    @inbounds p1 = points[indices[1]]
    @inbounds p2 = points[indices[2]]
    a = zero(nonzero_eltype(T))

    @inbounds for i in Iterators.drop(indices,2)
        p3 = points[i]
        a += area(p1,p2,p3)
        p2 = p3
    end

    return a
end

function area(points::AbstractVector{T},indices::VecOrTuple,x_period::Number,y_period::Number) where {T<:AbstractVec}
    @inbounds p1 = points[indices[1]]
    @inbounds p2 = closest(p1,points[indices[2]],x_period,y_period)
    a = zero(nonzero_eltype(T))

    @inbounds for i in Iterators.drop(indices,2)
        p3 = closest(p1,points[i],x_period,y_period)
        a += area(p1,p2,p3)
        p2 = p3
    end

    return a
end

"""
    in_polygon(point::Vec, vertices::AbstractVector{<:Vec}, indices) -> Bool

Whether `point` is inside a polygon with vertices given by `getindex.((vertices,),indices)`
"""
function in_polygon(p::Vec,points::AbstractVector{T},indices::VecOrTuple) where {T<:AbstractVec}
    @inbounds p1 = points[indices[1]]
    @inbounds p2 = points[indices[2]]

    inside = false
    @inbounds for i in Iterators.drop(indices,2)
        p3 = points[i]
        inside = in_triangle(p,p1,p2,p3)[1]
        inside && break
        p2 = p3
    end

    return inside
end

"""
    in_polygon_periodic(point::Vec, vertices::AbstractVector{<:Vec}, indices,x_period::Number,y_period::Number) -> Bool

Whether `point` is inside a polygon with vertices given by `getindex.((vertices,),indices) adjusted to the fact that some vertices positions might overflow the periodic domain.`
"""
function in_polygon(p::Vec,points::AbstractVector{T},indices::VecOrTuple,x_period::Number,y_period::Number) where {T<:AbstractVec}
    @inbounds p1 = points[indices[1]]
    @inbounds p2 = closest(p1,points[indices[2]],x_period,y_period)

    inside = false
    @inbounds for i in Iterators.drop(indices,2)
        p3 = closest(p1,points[i],x_period,y_period)
        inside = in_triangle(p,p1,p2,p3)[1]
        inside && break
        p2 = p3
    end

    return inside
end

"""
    is_obtuse(a::Vec,b::Vec,c::Vec) -> Bool

Whether the triangle formed by points `a`,`b`,`c` is obtuse
"""
function is_obtuse(a::Vec,b::Vec,c::Vec)
    ab = b-a
    nab = norm(ab)
    ac = c-a
    nac = norm(ac)

    #cos(0) <= 0
    A = acos((ac⋅ab)/(nab*nac))
    A >= π/2 && return true

    bc = c-b
    nbc = norm(bc)
    C = acos((bc⋅ac)/(nbc*nac))

    r = if ((C >= π/2) || (C+A <= π/2))
        true
    else
        false
    end
    return r
end

"""
    circle_edge_intersection(vertex_1::Vec2D, vertex_2::Vec2D, center::Vec2D, r::Number) -> Vec2D

Given edge vertices `vertex_1` and `vertex_2` returns the intersection point of the edge with a circle with radius `r` centered at `center` .
"""
function circle_edge_intersection(p1,p2,center,radius::Number)
    term1 = p1-p2
    term2 = p2-center
    a = term1⋅term1
    b = 2*term1⋅term2
    c = term2⋅term2 - radius*radius
    sqrt_bterm = sqrt(b*b -4*a*c)

    α1,α2 = ((-b + sqrt_bterm)/2a, (-b - sqrt_bterm)/2a )

    if ((0.0 <= α1) && (α1 <= 1.0)) 
        α = α1 
    elseif ((0.0 <= α2) && (α2 <= 1.0)) 
        α = α2
    else
        error("Circle doesn't cross edge")
    end

    return α*p1 + (1-α)*p2
end

end