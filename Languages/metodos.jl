using Pkg
Pkg.add("SpecialFunctions")
Pkg.add("DataFrames")

using SpecialFunctions
using DataFrames

f(x,y) = -(100+(1/(x^2)))*y
y_0 = -0.24593576
z_0 = -0.55769344

#----------------Metodo RK4----------------------
function coeficientes(h, x, y, z)
  k1 = h*z
  l1 = h*f(x,y)
  k2 = h*(z+l1/2)
  l2 = h*f(x+h/2, y+k1/2)
  k3 = h*(z+l2/2)
  l3 = h*f(x+h/2, y+k2/2)
  k4 = h*(z+l3)
  l4 = h*f(x+h, y+k3)

  c = [l1+2*l2+2*l3+l4, k1+2*k2+2*k3+k4]
  return c
end
#------------------------------------------------
#----------------Metodo RK6----------------------
function termos(h, x, y, z)
  k0 = (h^2)*f(x, y)
  k1 = (h^2)*f(x+(h/4), y+(h*z/4)+(k0/32))
  k2 = (h^2)*f(x+(h/2), y+(h*z/2)-(k0/24)+(k1/6))
  k3 = (h^2)*f(x+(3*h/4), y+(3*h*z/4)+(3*k0/32)+(k1/8)+(k2/16))
  k4 = (h^2)*f(x+h, y+h*z+(3*k1/7)-(k2/14)+(k3/7))

  t = [7*k0+32*k1+12*k2+32*k3+7*k4, 7*k0+24*k1+6*k2+8*k3]

  return t
end
#-----------------------------------------------

p = [0.025, 0.25, 0.5]

#Tabelas para armazenar os resultados
df_025 = DataFrame()
df_25 = DataFrame()
df_5 = DataFrame()
#------------------------------------

for h in p
  x = range(1,10*pi,step=h)
  y_rk4 = zeros(length(x))
  z_rk4 = zeros(length(x))
  y_rk6 = zeros(length(x))
  z_rk6 = zeros(length(x))
  j = zeros(length(x))
  
  y_rk4[1] = y_0
  z_rk4[1] = z_0
  y_rk6[1] = y_0
  z_rk6[1] = z_0
  
  for i in 1:length(x)
    j[i] = sqrt(x[i])*besselj0(10*x[i])
  end
  
  for n in 2:length(x)
    c = coeficientes(h, x[n-1], y_rk4[n-1], z_rk4[n-1])
    z_rk4[n] = z_rk4[n-1] + (1/6)*c[1]
    y_rk4[n] = y_rk4[n-1] + (1/6)*c[2]
    
    t = termos(h, x[n-1], y_rk6[n-1], z_rk6[n-1])
    z_rk6[n] = z_rk6[n-1] + (1/(90*h))*t[1]
    y_rk6[n] = y_rk6[n-1] + h*z_rk6[n-1] + (1/90)*t[2]
  end


  #Apenas armazenando os dados para serem printados em forma de tabela
  if h==0.025
    df_025[!,"x"] = x
    df_025[!,"rk4"] = y_rk4
    df_025[!,"rk6"] = y_rk6
    df_025[!,"bessel"] = j
  end
  if h==0.25
    df_25[!,"x"] = x
    df_25[!,"rk4"] = y_rk4
    df_25[!,"rk6"] = y_rk6
    df_25[!,"bessel"] = j
  end
  if h==0.5
    df_5[!,"x"] = x
    df_5[!,"rk4"] = y_rk4
    df_5[!,"rk6"] = y_rk6
    df_5[!,"bessel"] = j
  end
  #------------------------------------------------------------------
end

#Printando os resultados
println("###############################Resultado para h=0.025###########################################")
println(df_025)
println("###############################Resultado para h=0.25############################################")
println(df_25)
println("###############################Resultado para h=0.5#############################################")
println(df_5)
#----------------------------