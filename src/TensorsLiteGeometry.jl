module TensorsLiteGeometry
using TensorsLite

export circumcenter, closest, possible_positions_periodic

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

@inline possible_positions_periodic(p::Vec,periods::Tuple{Vararg{T,N}}) where {T<:Vec,N} = @inline map(+,ntuple(x->p,Val{N}()), periods)

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

@inline closest(p::Vec,p2::Vec,xp::Number,yp::Number) = closest(p,possible_positions_periodic(p2,xp,yp))

end
