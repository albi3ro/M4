function findstate(state::Int,set::Array{Int})
    iimin=1
    iimax=length(set)

    if set[iimin] == state
        return iimin
    end
    if set[iimax] == state
        return iimax
    end


    found=false
    count=0

    while found==false && count < length(set)
        count+=1
        index=floor(Int,iimin+(iimax-iimin)/2)
        if state < set[index]
            iimax=index-1
        elseif state > set[index]
            iimin=index+1
        else
            found=true
            return index
        end
    end

    if found == false
        println("findstate never found a match")
        println("Are you sure the state is in that Array?")
    end

    return 0
end

function bulk(n::Int,mz::Int)

    nstates=2^n

    psi=collect(0:(nstates-1))

    powers2=collect(0:(n-1))
    powers2=2.^powers2

    mask=[0;powers2]+[powers2;0]
    mask=mask[2:end-1]

    m=zeros(psi)
    for ii in 1:nstates
        m[ii]=sum((psi[ii]&powers2)./(powers2))
    end

    ind=sortperm(m)
    m=(m[ind]-n/2)/2
    psi=psi[ind]

    ma=collect(0:.5:n/2.)-n/4

    psia=Array{Array{Int64}}(n+1)
    first=1
    last=1
    for ii in 1:(n+1)
        psia[ii]=psi[first:last]

        #now we compute them for the next round
        first=last+1
        last=last+binomial(n,ii)
    end


    dim=binomial(n,mz-1)
    M=spzeros(Float64,dim,dim);
    Psi=psia[mz]

    mp=sum(Psi[1]&powers2./powers2)
    for ii in 1:length(Psi)
        p=Psi[ii]
        for jj in 1:(n-1)
            flipped=p$mask[jj]
            if sum((flipped&powers2)./powers2) == mp
                index=findstate(flipped,Psi)
                M[ii,index]=.5
                M[index,ii]=.5
            end

        end
    end

    ei=eigvals(full(M))
    println(size(M))

    writedlm("eigvals.dat",ei)
    return 0
end

@allocated @time bulk(17,8)
