abstract type MyAbstractComplex end

struct MyComplex <: MyAbstractComplex
    x
    y
end 

struct MyComplex_better <: MyAbstractComplex
    x::Float64
    y::Float64
end

struct MyComplex_better_better{T} <: MyAbstractComplex
    x::T
    y::T
end

function Base.show(io::IO, c::MyComplex)
    print(io, "$(c.x) + $(c.y)i")
end

function Base.:(+)(l::MyComplex, r::MyComplex)
    x = l.x+r.x
    y = l.y + r.y
    return MyComplex(x,y)
end

function Base.:(+)(l::MyAbstractComplex, r::Int)
    x = l.x+r.x
    y = l.y
    return MyComplex_better(x,y)
end

function silly_add(l::MyAbstractComplex, r::Int)
    x = l.x + r
    y = l.y
    return typeof(l)(x,y)
end

function Base.:(+)(l::Int, r::Int)
    x = l.x+r
    y = l.y + r
    return MyComplex(x,y)
end


input_array = rand(100, 100)

"""for a 2D array output a 1D array of averages for that dimension"""
function small_array_calculation(inputmatrix)
    avg = Float64[]
    x,y = size(inputmatrix)
    for i in 1:x
        vector = @view inputmatrix[i,:]
        a = sum(vector)/length(vector)
        push!(avg, a)

    end
end 
