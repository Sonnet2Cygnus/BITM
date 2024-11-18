# BITM

This code demonstrates a bubble-induced turbulence model. The main functions and related equationss of the code are as follows:

## Equations

Blending function F1
```math
\mathrm{arg}1=\min \left( \min \left( \max \left( \frac{1}{\beta ^*}\frac{\sqrt{k}}{\omega y},500\frac{\nu}{y^2\omega} \right) ,\frac{4\alpha _{\omega 2}k}{CD_{k\omega}y^2} \right) ,10 \right) \,\, F1=\tanh\mathrm{(arg}1^4)
```

F2
```math
\mathrm{arg}2=\min \left( \max \left( \frac{2}{\beta ^*}\frac{\sqrt{k}}{\omega y},500\frac{\nu}{y^2\omega} \right) ,100 \right) \,\, F2=\tanh\mathrm{(arg}2^2)
```

F3
```math
\mathrm{arg}3=\min \left( 150\frac{\nu}{\omega y^2},10 \right) \qquad \,\,F3=1-\tanh\mathrm{(arg}3^4)
```

Turbulent viscosity
```math
\nu _t=\frac{a_1\cdot k}{\max\mathrm{(}a_1\cdot \omega ,b_1\cdot F_2\cdot \sqrt{S^2})}
```

Production
```math
P_k=\min\mathrm{(}G,c_1\cdot \beta ^*\cdot k\cdot \omega )
```
  
k production or dissipation matrix
```math
\left[ \mathrm{L} \right] ^3\times \frac{\left[ \mathrm{M} \right]}{\left[ \mathrm{L} \right] ^3}\times \frac{\left[ \mathrm{L} \right] ^2}{\left[ \mathrm{T} \right] ^2}/\left[ \mathrm{T} \right] =\frac{\left[ \mathrm{M} \right] \cdot \left[ \mathrm{L} \right] ^2}{\left[ \mathrm{T} \right] ^3}
```

omega matrix
```math
\left[ \mathrm{L} \right] ^3\times \frac{\left[ \mathrm{M} \right]}{\left[ \mathrm{L} \right] ^3}\times \frac{1}{\left[ \mathrm{T} \right]}/\left[ \mathrm{T} \right] =\frac{\left[ \mathrm{M} \right]}{\left[ \mathrm{T} \right] ^2}
```

k equation
```math
\frac{\partial (\alpha \rho k)}{\partial t}+\nabla \cdot (\alpha \rho k\mathbf{U})=\nabla \cdot (\alpha \rho D_{k}^{\mathrm{eff}}\nabla k)+\alpha \rho P_k(G)-\frac{2}{3}\alpha \rho k(\nabla \cdot \mathbf{U})-\alpha C_{\mu}\rho ^L\omega k+C_I\frac{3}{4d_p}C_D\rho ^L\alpha ^G\left| u^G-u^L \right|\left( u^G-u^L \right) \left( u^G-u^L \right) \qquad 
```
```math
S_k=C_I\frac{3}{4d_p}C_D\rho ^L\alpha ^G\left| u^G-u^L \right|\left( u^G-u^L \right) \left( u^G-u^L \right) 
```

omega equation
```math
\frac{\partial (\alpha \rho \omega )}{\partial t}+\nabla \cdot (\alpha \rho \omega \mathbf{U})=\nabla \cdot (\alpha \rho D_{\omega}^{\mathrm{eff}}\nabla \omega )+\alpha \rho \gamma \min \left( GbyNu,\frac{c_1}{a_1}\beta ^*\omega \cdot \max\mathrm{(}a_1\omega ,b_1F_{23}\sqrt{S^2}) \right) -\frac{2}{3}\alpha \rho \gamma \omega (\nabla \cdot \mathbf{U})-\alpha ^LC_{\omega D}\rho ^L\omega ^2+\left( \frac{1}{C_{\mu}k}S_{\epsilon}-\frac{\omega}{k}S_k \right) \qquad where\,\,S_{\epsilon}=C_{\epsilon}\frac{S_k}{\tau}
```
