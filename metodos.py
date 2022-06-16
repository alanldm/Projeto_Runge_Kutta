import streamlit as st
st.set_page_config(layout="wide")
import numpy as np
import matplotlib.pyplot as plt
import scipy.special
import plotly.graph_objects as go
import plotly.express as px


#------------------------------------Texto em LaTeX-----------------------------------------------------------
st.markdown("<h1 style='text-align: center;'>Métodos computacionais</h1>", unsafe_allow_html=True)

st.latex(r'''
\frac{d^2y(x)}{dx^2}=-\left(100+\frac{1}{x^2}\right)y(x)=f(x,y)
''')

st.latex(r'''
\text{Dado que \textbf{y(1)=-0.24593576} e \textbf{y'(1)=-0.55769344}}
''')

st.latex(r'''
\text{\textbf{Intervalo: [1,10π]}}
''')

st.latex(r'''
z=y' \rightarrow z(1)=y'(1)=-0.55769344

\\

z'=y''=-\left(100+\frac{1}{x^2}\right)y(x)
''')

st.markdown("<h1 style='text-align: left;'>Equações RK4:</h1>", unsafe_allow_html=True)

col1, col2 = st.columns(2)

with col1:
  st.latex(r'''
  y_{n+1}=y_n+\frac{1}{6}(k_1+2k_2+2k_3+k_4)
  ''')
  st.latex(r'''
  k_1=hz_n
  \newline
  k_2=h\left(z_n+\frac{l_1}{2}\right)
  \newline
  k_3=h\left(z_n+\frac{l_2}{2}\right)
  \newline
  k_4=h(z_n+l_3)
  ''')

with col2:
  st.latex(r'''
  z_{n+1}=z_n+\frac{1}{6}(l_1+2l_2+2l_3+l_4)
  ''')
  
  st.latex(r'''
  l_1=hf(x_n,y_n)
  \newline
  l_2=hf\left(x_n+\frac{h}{2},y_n+\frac{k_1}{2}\right)
  \newline
  l_3=hf\left(x_n+\frac{h}{2},y_n+\frac{k_2}{2}\right)
  \newline
  l_4=hf(x_n+h,y_n+k_3)
  ''')

st.markdown("<h1 style='text-align: left;'>Equações RK6:</h1>", unsafe_allow_html=True)

st.latex(r'''
y_{n+1}=y_n+hz_n+\frac{1}{90}(7k_0+24k_1+6k_2+8k_3)
''')

st.latex(r'''
z_{n+1}=z_n+\frac{1}{90h}(7k_0+32k_1+12k_2+32k_3+7k_4))
''')

st.latex(r'''
k_0=h^2f(x_n, y_n)
\\
k_1=h^2f\left(x_n+\frac{h}{4}, y_n+\frac{hz_n}{4}+\frac{k_0}{32}\right)
\\
k_2=h^2f\left(x_n+\frac{h}{2}, y_n+\frac{hz_n}{2}-\frac{k_0}{24}+\frac{k_1}{6}\right)
\\
k_3=h^2f\left(x_n+\frac{3h}{4}, y_n+\frac{3hz_n}{4}+\frac{3k_0}{32}+\frac{k_1}{8}+\frac{k_2}{16}\right)
\\
k_4=h^2f\left(x_n+h, y_n+hz_n+\frac{3k_1}{7}-\frac{k_2}{14}+\frac{k_3}{7}\right)
''')

#-------------------------------------------------------------------------------------------------------------------
y_0 = -0.24593576
z_0 = -0.55769344

def f(x, y):
  return -(100+(1/x**2))*y

def coeficientes(h, x, y, z):
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

#Função dos coeficientes para o RK6
def termos(h, x, y, z):
  k0 = (h**2)*f(x, y)
  k1 = (h**2)*f(x+(h/4), y+(h*z/4)+(k0/32))
  k2 = (h**2)*f(x+(h/2), y+(h*z/2)-(k0/24)+(k1/6))
  k3 = (h**2)*f(x+(3*h/4), y+(3*h*z/4)+(3*k0/32)+(k1/8)+(k2/16))
  k4 = (h**2)*f(x+h, y+h*z+(3*k1/7)-(k2/14)+(k3/7))

  t = [7*k0+32*k1+12*k2+32*k3+7*k4, 7*k0+24*k1+6*k2+8*k3]

  return t

hs = np.arange(0.01, 1.01, 0.01)
x = []
yrk4 = []
yrk6 = []
y_exato = []

#Gerando o y exato com a função de bessel
x_01 = np.arange(1, 10*np.pi, 0.01)
for i in x_01:
  y_exato.append(np.sqrt(i)*scipy.special.j0(10*i))
#-----------------------------------------

#Gerando todas as variações de x e os valores de y
for h in hs:
  x_ = np.arange(1, 10*np.pi, h)
  x.append(x_)

  y_rk4 = np.zeros([len(x_),1], dtype=float)
  z_rk4 = np.zeros([len(x_),1], dtype=float)
  y_rk4[0] = y_0
  z_rk4[0] = z_0

  y_rk6 = np.zeros([len(x_),1], dtype=float)
  z_rk6 = np.zeros([len(x_),1], dtype=float)
  y_rk6[0] = y_0
  z_rk6[0] = z_0

  for n in range(1, len(x_)):
    c = coeficientes(h, x_[n-1], y_rk4[n-1], z_rk4[n-1])
    z_rk4[n] = z_rk4[n-1] + (1/6)*c[0]
    y_rk4[n] = y_rk4[n-1] + (1/6)*c[1]

    t = termos(h, x_[n-1], y_rk6[n-1], z_rk6[n-1])
    z_rk6[n] = z_rk6[n-1] + (1/(90*h))*t[0]
    y_rk6[n] = y_rk6[n-1] + h*z_rk6[n-1] + (1/90)*t[1]

  yrk4.append(y_rk4)
  yrk6.append(y_rk6)
#------------------------------------------------

def gera_grafico(level):
    y_rk4 = []
    y_rk6 = []

    for i in range(0,len(yrk4[level])):
        y_rk4.append(float(yrk4[level][i]))
        y_rk6.append(float(yrk6[level][i]))
        
    fig = go.Figure()
    fig.add_trace(go.Scatter(name='Exato', x=x[0],y=y_exato, marker_color='gray'))
    fig.add_trace(go.Scatter(name='RK4', x=x[level],y=y_rk4, marker_color='blue'))
    fig.update_layout(title='Solução pelo RK4',  xaxis_title='x', yaxis_title='y', yaxis_range=[-0.3,0.3])

    fig2 = go.Figure()
    fig2.add_trace(go.Scatter(name='Exato', x=x[0],y=y_exato, marker_color='gray'))
    fig2.add_trace(go.Scatter(name='RK6', x=x[level],y=y_rk6, marker_color='blue'))
    fig2.update_layout(title='Solução pelo RK6', xaxis_title='x', yaxis_title='y', yaxis_range=[-0.3,0.3])

    col1, col2 = st.columns(2)

    with col1:
      st.plotly_chart(fig, use_container_width=True, width=1200)

    with col2:
      st.plotly_chart(fig2, use_container_width=True, width=900)

level = st.slider("Selecione o valor de h", 0, 99)
texto_level = [str(round(i, 2)) for i in np.arange(0.01, 1.01, 0.01)]
st.text('h: {}'.format(texto_level[level]))

gera_grafico(level)


col1, col2, col3 = st.columns([1,2,1])

with col2:
  video_file = open('video_rk46.mp4', 'rb')
  video_bytes = video_file.read()

  st.video(video_bytes)