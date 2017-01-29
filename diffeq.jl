module diffeq

  export Euler, RK2, RK4, Solver

  function Euler(f::Array{Function,1},t0::Float64,x::Array{Float64,1},h::Float64)
      d=length(f)
      xp=copy(x)
      for ii in 1:d
          xp[ii]+=h*f[ii](t0,x)
      end

      return t0+h,xp
  end

  function RK2(f::Array{Function,1},t0::Float64,x::Array{Float64,1},h::Float64)
      d=length(f)

      xp=copy(x)
      xk1=copy(x)

      for ii in 1:d
          xk1[ii]+=f[ii](t0,x)*h/2
      end
      for ii in 1:d
          xp[ii]+=f[ii](t0+h/2,xk1)*h
      end

      return t0+h,xp
  end

  function RK4(f::Array{Function,1},t0::Float64,x::Array{Float64,1},h::Float64)
      d=length(f)

      k1=zeros(x)
      k2=zeros(x)
      k3=zeros(x)
      k4=zeros(x)

      for ii in 1:d
          k1[ii]=h*f[ii](t0,x)
      end
      for ii in 1:d
          k2[ii]=h*f[ii](t0+h/2,x+k1/2)
      end
      for ii in 1:d
          k3[ii]=h*f[ii](t0+h/2,x+k2/2)
      end
      for ii in 1:d
          k4[ii]=h*f[ii](t0+h,x+k3)
      end

      return t0+h,x+(k1+2*k2+2*k3+k4)/6
  end

  function Solver(f::Array{Function,1},Method::Function,t0::Float64,
          x0::Array{Float64,1},h::Float64,N::Int64)
      d=length(f)
      ts=zeros(Float64,N+1)
      xs=zeros(Float64,d,N+1)

      ts[1]=t0
      xs[:,1]=x0

      for i in 2:(N+1)
          ts[i],xs[:,i]=Method(f,ts[i-1],xs[:,i-1],h)
      end

      return ts,xs
  end

end
