export convn, fastconv

##############################################
# Generic convn function using direct method for computing convolutions:
# Accelerated Convolutions for Efficient Multi-Scale Time to Contact Computation in Julia
# Alexander Amini, Alan Edelman, Berthold Horn
##############################################

@generated function convn(E::Array{T,N}, k::Array{T,N}) where {T,N}
    quote
        sizeThreshold = 21;
        if length(k) <= sizeThreshold || $N > 2
            #println("using direct")
            retsize = [size(E)...] + [size(k)...] .- 1
            retsize = tuple(retsize...)
            ret = zeros(T, retsize)

            convn!(ret,E,k)
            return ret
        else
            #println("using fft1")
            return conv(E,k)
        end
    end
end

# direct version (do not check if threshold is satisfied)
@generated function fastconv(E::Array{T,N}, k::Array{T,N}) where {T,N}
    quote

        retsize = [size(E)...] + [size(k)...] .- 1
        retsize = tuple(retsize...)
        ret = zeros(T, retsize)

        convn!(ret,E,k)
        return ret

    end
end


# in place helper operation to speedup memory allocations
@generated function convn!(out::Array{T}, E::Array{T,N}, k::Array{T,N}) where {T,N}
    quote
        @inbounds begin
            @nloops $N x E begin
                @nloops $N i k begin
                    (@nref $N out d->(x_d + i_d - 1)) += (@nref $N E x) * (@nref $N k i)
                end
            end
        end
        return out
    end
end
