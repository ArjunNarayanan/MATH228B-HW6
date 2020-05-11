using PyPlot

function characteristic_eq(xi,x,t)
    return t*sin(pi*xi) + xi - x
end

function derivative_characteristic_eq(xi,x,t)
    return pi*t*cos(pi*xi) + 1.0
end

function newton_iterate(F,Fprime,x0;atol=1e-12,maxiter=100)
    xnext = x0
    val = F(xnext)
    dval = Fprime(xnext)
    count = 1
    while abs(val) > atol && count < maxiter
        xnext = xnext - val/dval
        val = F(xnext)
        dval = Fprime(xnext)
        count += 1
    end
    if count == maxiter
        error("Newton iteration failed to converge")
    end
    return xnext
end

function initial_condition(x)
    return sin(pi*x)
end

function characteristic_base_points(xrange,t)
    @assert xrange[1] ≈ -1.0

    xirange = similar(xrange)
    xirange[1] = newton_iterate(xi->characteristic_eq(xi,xrange[1],t),
        xi->derivative_characteristic_eq(xi,xrange[1],t),xrange[1])

    for i = 2:length(xrange)
        x = xrange[i]
        xirange[i] = newton_iterate(xi->characteristic_eq(xi,x,t),
            xi->derivative_characteristic_eq(xi,x,t),xirange[i-1])
    end
    return xirange
end

xrange = -1:1e-1:0
t = 0.5
xirange = characteristic_base_points(xrange,t)
