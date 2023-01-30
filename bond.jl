using PyPlot

r = LinRange(0.97,2,100);
p = -1 ./r.^6+1 ./r.^12;plot(r,p)
f = -diff(p)
push!(f,f[end])

plot(r,p)
ax1 = gca()
ax2 = ax1.twinx()
ax2.plot(r,f)
show()
