module Lattices
export Lattice, MakeLattice, PlotNeighbors
    using PyPlot

    type Lattice
        name::ASCIIString
        l::Int16
        dim::Int8
        a::Array
        unit::Array
        N::Int16
        X::Array
        nnei::Int8
        neigh::Array
    end
    
    function armod(x,y)
        return mod(x-1,y)+1
    end
    
    
    function MakeLattice(name,l)
        if(name=="Square")
            return MakeSquare(l)
        elseif(name=="Chain")
            return MakeChain(l)
        elseif(name=="Honeycomb")
            return MakeHoneycomb(l)
            
            
            #elseif(name=="Checkerboard")
            #elseif(name=="Pyrochlore")
            #elseif(name=="Kagome")
        else
            println("I didnt give a good name, making a 2x2 square")
            return MakeSquare(2)
        end
    
    end
    
    function MakeSquare(l)
        N=l^2;
        d=2;
        X=Array{Int16}(N,2);
        nnei=4;
        neigh=Array{Int16}(N,4);
    
        X[:,1]=armod(collect(1:N),l);
        X[:,2]=ceil(collect(1:(N))/l);
        a=[[1 0]
            [0 1]];
        unit=[0,0];
        
        neigh[:,1]=armod(X[:,1]+1,l)+l*(X[:,2]-1);
        neigh[:,2]=armod(X[:,1]+l-1,l)+l*(X[:,2]-1);
        neigh[:,3]=X[:,1]+mod(X[:,2],l)*l;
        neigh[:,4]=X[:,1]+mod(X[:,2]+l-2,l)*l;
                
        return Lattice("Square",l,d,a,unit,N,X,nnei,neigh)
    end 
    
    function MakeChain(l)
        N=l
        d=1
        nnei=2
        neigh=Array{Int16}(N,2)
        a=[1];
        unit=[0];
                
        X=collect(1:l);
        
        neigh[:,1]=armod(collect(2:(l+1)),l);
        neigh[:,2]=armod(collect(l:(2*l-1)),l);
            
        return Lattice("Chain",l,d,a,unit,N,X,nnei,neigh)
    end
    
    function MakeHoneycomb(l)
        d=2;
        nnei=3;
        N=2*l^2;
        a=[[2*cos(pi/6) 0]
            [cos(pi/6) 1+sin(pi/6)]];
        unit=[[0 0]
                [cos(π/6) sin(π/6)]];

        return GeneralMultiUnit2D("Honeycomb",l,d,a,unit,nnei)
    end
    
    function GeneralMultiUnit2D(name,l,d,a,unit,nnei)
        ncell=size(a)[1];
        N=l*l*ncell;
        
        a1=repeat(a[1,:],outer=[ncell,1]);
        a2=repeat(a[2,:],outer=[ncell*l,1]);
        
        X=Array{Float64}(N,2);
        # Here we are actually calculating the positions for every site
        for i in 1:l    #for the first row
            X[ncell*i-ncell+1:ncell*i,:]=unit+(i-1)*a1;
        end
        
        for j in 2:l    #copying the first row into the first layer
            X[(ncell*l*(j-1)+1):(ncell*l*(j-1)+ncell*l),:]=X[1:ncell*l,:]+(j-1)*a2;
        end
        
        neigh=Array{Int16}(N,nnei)
        
        return Lattice(name,l,d,a,unit,N,X,nnei,neigh)
    end
    
    
    
    function PlotNeighbors(lt::Lattice)
        fig=gcf()
        if(lt.dim==2)
            for i in 1:lt.N
                for j in 1:lt.nnei
                    xx=[lt.X[i,1], lt.X[lt.neigh[i,j],1] ]
                    yx=[lt.X[i,2], lt.X[lt.neigh[i,j],2] ]
                    println(xx,' ',yx)
                    plot(xx,yx)
                end
            end
        elseif(lt.dim==3)
            for i in 1:lt.N
                for j in 1:lt.nnei
                    xx=[lt.X[i,1], lt.X[lt.neigh[i,j],1] ]
                    yx=[lt.X[i,2], lt.X[lt.neigh[i,j],2] ]
                    zx=[lt.X[i,3], lt.X[lt.neigh[i,j],3] ]
                    println(xx,' ',yx)
                    plot3D(xx,yx)
                end
            end
        else
            println("Dimension not 2 or 3")
            println("Not plotting")
        end
    end

end






