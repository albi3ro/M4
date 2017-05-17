function GenerateHermite(n)
    Hermite=Function[]

    push!(Hermite,x->1);
    push!(Hermite,x->2*x);

    for ni in 3:n
        push!(Hermite,x->2.*x.*Hermite[ni-1](x).-2.*n.*Hermite[ni-2](x))
    end
    return Hermite
end

m=1
ω=1
ħ=1

Ψ(n,x)=1/sqrt(factorial(n)*2^n)*(m*ω/(ħ*π))^(1/4)*exp(-m*ω*x^2/(2*ħ))*Hermite[n](sqrt(m*ω/ħ)*x)
