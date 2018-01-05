module KeplerCoords

    export coordOrb, coordCart, OrbtoCart, CarttoOrb, +

    mutable struct coordCart
        x::Float64
        y::Float64
        z::Float64

        vx::Float64
        vy::Float64
        vz::Float64
    end

    mutable struct coordOrb
        a::Float64 # Semimajor Axis
        e::Float64  # eccentricity
        i::Float64  # inclination

        ω::Float64 # \omega Argument of periapsis
        Ω::Float64 # \Omega Longitude of ascending node

        ν::Float64 # \nu True Anomaly
    end

    import Base.+, Base.*

    function Base.:+(C1::coordCart,C2::coordCart)
        xn= C1.x+C2.x
        yn= C1.y+C2.y
        zn= C1.z+C2.z

        vxn=C1.vx+C2.vx
        vyn=C1.vy+C2.vy
        vzn=C1.vz+C2.vz

        return coordCart(xn,yn,zn,vxn,vyn,vzn)
    end

    function Base.:*(c::Float64,C::coordCart)

        return coordCart(c*C.x,c*C.y,c*C.z,c*C.vx,c*C.vy,c*C.vz)
    end

    function Rotmat(θ::Float64)
        R=Array{Float64,2}(2,2)
        R[1,1]=cos(θ)
        R[2,2]=cos(θ)
        R[1,2]=-sin(θ)
        R[2,1]=sin(θ)
        return R
    end

    function OrbtoCart(O::coordOrb,μ::Float64)
        r=(O.a*(1-O.e^2))/(1+O.e*cos(O.ν))
        h=sqrt(μ*O.a*(1-O.e^2))
        v=h/r

        R1=eye(Float64,3)
        R1[[1,2],[1,2]]=Rotmat(O.ω+O.ν)

        R2=eye(Float64,3)
        R2[[1,3],[1,3]]=Rotmat(O.i)

        R3=eye(Float64,3)
        R3[[1,2],[1,2]]=Rotmat(O.Ω)

        pos=r*R3*R2*R1*[1,0,0]
        vel=v*R3*R2*R1*[0,1,0]

       return coordCart(pos[1],pos[2],pos[3],vel[1],vel[2],vel[3])
    end

    function CarttoOrb(C::coordCart,μ::Float64)
        r=sqrt(C.x^2+C.y^2+C.z^2)
        v=sqrt(C.vx^2+C.vy^2+C.vz^2)
        hvec=cross([C.x,C.y,C.z],[C.vx,C.vy,C.vz])
        h=norm(hvec)

        E=v^2/2-μ/r

        a=-μ/(2*E)
        e=sqrt(1-(h^2)/(a*μ))

        i=acos(hvec[3]/h)

        Ω=atan2(hvec[1],-hvec[2])

        ων=atan2(C.z/sin(i),C.x*cos(Ω)+C.y*sin(Ω))
        ν=acos((a*(1-e^2)-r)/(e*r))


       return coordOrb(a,e,i,Ω,ω,ν)
    end
end
