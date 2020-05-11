using PyPlot


function initial_condition(x)
    return sin(pi*x)
end



fig,ax = PyPlot.subplots()
x = range(-1,stop=1,length=30)
N = length(x)
y = zeros(N)
u = initial_condition.(x)
v = ones(N)
ax.quiver(x,y,u,v,scale=3)
fig
