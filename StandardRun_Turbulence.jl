################################################################
# IBM for evolution of cells producing toxins as a public good #
#               Author: Elias Ehrlich                          #
#                 Date: 18/11/2019                             #
################################################################
# Packages
using Plots
using Colors
using Random
using Distances
using Statistics
using StatsBase
using LaTeXStrings
using Measures
using DelimitedFiles
using ColorSchemes

include("./FastConv/FastConv")
include("./IBM_Functions_Turbulence")

# Random.seed!(1234)

# Parameters
## Resolution: Time step [days], grid cell size [cm]
dt=1/12
dx=0.5

## Space width [cm], number of simulation days for phase before and after invasion
L=50
tmax1=1000
tmax2=11000

## Number of grid cells on each axis, number of non-toxic and toxic cells, number of time steps, time step at start and end of each phase [dt]
Ngrid=ceil(Int, L/dx)
Nnon=10000
Ntox=0
Nt1=ceil(Int, tmax1/dt)
Nt2=ceil(Int, tmax2/dt)

t_start1=1
t_end1=copy(Nt1)
t_start2=Nt1+1
t_end2=Nt1+Nt2

## Initial, spatially averaged nutrient concentration [fmol N cm^-2]
K=9.6*dx^2

## Decay rate of toxins [d^-1], toxin leakage rate [amol T cell^-1 d^-1]
λ=0.2
Q=3.4

## Max. growth rate, non-grazing mortality rate, max. grazing (loss) rate [d^-1]
rmax_tox=0.68
rmax_non=0.7
dtox=0.05
dnon=0.05
Gmax_tox=0.3
Gmax_non=0.3

## Non-grazing mortality probability
mtox=-exp(-dtox*dt)+1
mnon=-exp(-dnon*dt)+1

## Half-saturation constants for nutrients [fmol N cm^-2] and toxin effect [amol T cm^-2]
HC=6.1*dx^2
HT=30.0*dx^2

## Cell quota [fmol N cell^-1]
qN=0.97

## Diffusivity of nutrients, toxins and cells [cm^2 d^-1]
Dn=0.864
Dc=0.05

## Turbulence amplitude [unit: number of grid cells]
U=10.0/dx # maximum movement by turbulence per time step dt = U/2 --> U in [cm]: U=20 cm

## Average diffusion distance of cells [units: dx=number of grid cells] per time step dt
Δ=sqrt(2*Dc*dt)/dx

## Diffusion kernel (1 or 2 dimension). Warning! Dn/dx/dx*dt must be <1/3
α=Dn/dx/dx*dt

if α>1/3
    error("Time step to large")
end
kernel=α*[1,-2,1] + [0,1,0]

kernel_2D=zeros(3,3)
kernel_2D[2,:]=kernel
for i in 1:3
    k2D = FastConv.convn(kernel_2D[:,i],kernel) # convn is from Package FastConv
    kernel_2D[:,i]=k2D[2:end-1]
end

# Initialization of nutrient conc., toxin conc. and cell positions
InitC=K-qN*(Nnon+Ntox)/(Ngrid*Ngrid)
C=fill(InitC,Ngrid,Ngrid)#2*K*rand(Ngrid,Ngrid)
T=zeros(Ngrid,Ngrid)

xnon=ceil.(Int,rand(Nnon)*Ngrid)
ynon=ceil.(Int,rand(Nnon)*Ngrid)
xtox=ceil.(Int,rand(Ntox)*Ngrid)
ytox=ceil.(Int,rand(Ntox)*Ngrid)

PopDens=zeros(Int64,t_end2+1,2)     # Population density of non-toxic and toxic cells
PopDens[1,:]=[Nnon,Ntox]            # Initial population densities

Patchiness=zeros(Float64,t_end2+1,2)
NON=IBM_Functions_Turbulence.CellCount(xnon,ynon,Ngrid)
TOX=IBM_Functions_Turbulence.CellCount(xtox,ytox,Ngrid)
Patchiness[1,1]=IBM_Functions_Turbulence.LloydPatchiness(NON) # Initial patchiness nontoxics
Patchiness[1,2]=IBM_Functions_Turbulence.LloydPatchiness(TOX) # Initial patchiness nontoxics

Correlation=zeros(Float64,t_end2+1,6)
Correlation[1,:]=IBM_Functions_Turbulence.CorrAllCombinations(NON,TOX,T,C)

T_exposed=zeros(Float64,t_end2+1,2)
mean_T=mean(T)
T_exposed[1,1]=cov(vec(NON),vec(T))/mean(NON)+mean_T # Mean toxin conc. experienced by a non-toxic cell
T_exposed[1,2]=cov(vec(TOX),vec(T))/mean(TOX)+mean_T # Mean toxin conc. experienced by a toxic cell

G_reduced=zeros(Float64,t_end2+1,2)
PosTox=IBM_Functions_Turbulence.LinIndex(xtox,ytox,Ngrid)
PosNon=IBM_Functions_Turbulence.LinIndex(xnon,ynon,Ngrid)
Gtox=@. Gmax_tox*(1.0-T[PosTox]/(HT+T[PosTox]))
Gnon=@. Gmax_non*(1.0-T[PosNon]/(HT+T[PosNon]))
G_reduced[1,1]=mean(1 .- Gnon ./Gmax_non) # Mean reduction of grazing rate experienced by non-toxic cells
G_reduced[1,2]=mean(1 .- Gtox ./Gmax_tox) # Mean reduction of grazing rate experienced by toxic cells

# Run simulation
## before invasion, only resident non-toxic cells
@time C, T, xnon, ynon, xtox, ytox, PopDens, Patchiness, Correlation, T_exposed, G_reduced = IBM_Functions_Turbulence.Sim_Tox_Graz_PDE(C,T,xnon,ynon,xtox,ytox,PopDens,Patchiness,Correlation,T_exposed,G_reduced,Ngrid,t_start1,t_end1,λ,Q,
                                                            rmax_tox,rmax_non,mtox,mnon,Gmax_tox,Gmax_non,HC,HT,qN,dt,kernel_2D,Δ,U)

## invasion event
if (tmax2!=0)
    ReplacePos=sample(1:length(ynon),100,replace=false) #1000
    ReplacePos=sort(ReplacePos)
    xtox=copy(xnon[ReplacePos])
    ytox=copy(ynon[ReplacePos])
    deleteat!(xnon,ReplacePos)
    deleteat!(ynon,ReplacePos)
    PopDens[t_start2,:]=[length(xnon),length(xtox)]                    # Invasion event - Overwrite population densities at t_end1+1=t_start2 (1000 days)

    NON=IBM_Functions_Turbulence.CellCount(xnon,ynon,Ngrid)
    TOX=IBM_Functions_Turbulence.CellCount(xtox,ytox,Ngrid)
    Patchiness[t_start2,1]=IBM_Functions_Turbulence.LloydPatchiness(NON) # Invasion event - overwrite patchiness nontoxics at t_start2
    Patchiness[t_start2,2]=IBM_Functions_Turbulence.LloydPatchiness(TOX) # Invasion event - overwrite patchiness toxics at t_start2

    Correlation[t_start2,:]=IBM_Functions_Turbulence.CorrAllCombinations(NON,TOX,T,C) # Invasion event - overwrite Pearson's correlation coefficients at t_start2

    mean_T=mean(T)
    T_exposed[t_start2,1]=cov(vec(NON),vec(T))/mean(NON)+mean_T # Mean toxin conc. experienced by a non-toxic cell
    T_exposed[t_start2,2]=cov(vec(TOX),vec(T))/mean(TOX)+mean_T # Mean toxin conc. experienced by a toxic cell

    PosTox=IBM_Functions_Turbulence.LinIndex(xtox,ytox,Ngrid)
    PosNon=IBM_Functions_Turbulence.LinIndex(xnon,ynon,Ngrid)
    Gtox=@. Gmax_tox*(1.0-T[PosTox]/(HT+T[PosTox]))
    Gnon=@. Gmax_non*(1.0-T[PosNon]/(HT+T[PosNon]))
    G_reduced[t_start2,1]=mean(1 .- Gnon ./Gmax_non) # Mean reduction of grazing rate experienced by non-toxic cells
    G_reduced[t_start2,2]=mean(1 .- Gtox ./Gmax_tox) # Mean reduction of grazing rate experienced by toxic cells

    @time C, T, xnon, ynon, xtox, ytox, PopDens, Patchiness, Correlation, T_exposed, G_reduced = IBM_Functions_Turbulence.Sim_Tox_Graz_PDE(C,T,xnon,ynon,xtox,ytox,PopDens,Patchiness,Correlation,T_exposed,G_reduced,Ngrid,t_start2,t_end2,λ,Q,
                                                                 rmax_tox,rmax_non,mtox,mnon,Gmax_tox,Gmax_non,HC,HT,qN,dt,kernel_2D,Δ,U)
end

# Plot
p1=scatter(xnon.*dx,ynon.*dx,xlabel=L"x\, [cm]",ylabel=L"y\, [cm]",
    mc=:steelblue,msc=:black,ms=1.2,msw=0.0,
    reuse=false,legend = false,
    xlims=(0,L),ylims=(0,L),title=latexstring("\$t=$(tmax1+tmax2) d\$"),
    guidefontsize=8,titlefontsize=8,left_margin=2mm)
scatter!(xtox.*dx,ytox.*dx,mc=:red,msc=:black,ms=0.8,msw=0.0)
p2=heatmap(dx:dx:L,dx:dx:L,transpose(C)./(dx^2),clims=(minimum(C)/(dx^2),maximum(C)/(dx^2)),title=L"N-conc.\, [fmol\,N\,cm^{-2}]",titlefontsize=8,xlims=(0,L),ylims=(0,L),
    xlab=L"x\,[cm]",ylab=L"y\,[cm]",guidefontsize=8,left_margin=2mm,c=:inferno)
p3=heatmap(dx:dx:L,dx:dx:L,transpose(T)./(dx^2),clims=(0,maximum(T)/(dx^2)),title=L"T-conc.\, [amol\,T\,cm^{-2}]",titlefontsize=8,xlims=(0,L),ylims=(0,L),
    xlab=L"x\,[cm]",ylab=L"y\,[cm]",guidefontsize=8,left_margin=2mm,c=:inferno)
p4=plot(0:dt:(tmax1+tmax2),PopDens[1:(t_end2+1),:],reuse=false,legend=false,xlabel=L"t\, [d]",ylabel=L"Abundance",
    xlims=(0,tmax1+tmax2),ylim=(10^2,10^5),yaxis=:log,
    linecolor=[:steelblue :red],linewidth=0.5,guidefontsize=8, label = ["non-toxic" "toxic"],
    left_margin=2mm,xticks = (0:5000:(tmax1+tmax2),["0","5000","10000","","20000"]))

l=@layout [a{0.44w} b
           c{0.44w} d]
pall=plot(p4,p3,p1,p2,layout=l)
plot!(size=(500,400),dpi=600)
fontsize=10
annotate!( -11, 131, text("A", :left, fontsize), subplot=3)
annotate!( -11, 55, text("B", :left, fontsize), subplot=3)
annotate!( -11, 55, text("C", :left, fontsize), subplot=2)
annotate!( -11, 55, text("D", :left, fontsize), subplot=4)
annotate!( -60, 55, text(L"non-toxic", :left, 8, color=:black), subplot=2)
annotate!( -30, 55, text(L"toxic", :left, 8, color=:black), subplot=2)
annotate!(-65, 54, text(L"\textbf{-}", :left, 20, color=:steelblue), subplot=2)
annotate!(-35, 54, text(L"\textbf{-}", :left, 20, color=:red), subplot=2)

display(pall)
