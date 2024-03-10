$$
\begin{align*}
&\begin{cases}
\partial_t c - \nabla(D \cdot \nabla c) - \alpha \cdot c(1-c) = 0 \quad \text{in} \quad \Omega \\
D \cdot \nabla c \cdot n = 0 \quad \text{on} \quad \partial\Omega \\
c(t=0) = c_0 \quad \text{in} \quad \Omega \\
\end{cases} \\
&\text{Weak formulation:} \\
&\text{Find } c \in V \text{ such that } R(c)(v) = 0 \quad \forall v \in V' \\
&R(c)(v) := \int_\Omega (\partial_t c + \nabla(D \cdot \nabla c) - \alpha \cdot c(1-c))v \, dV \\
&\text{Using Green's formula, the second term becomes:} \\
&R(c)(v) := \int_\Omega (\partial_t c + D \cdot \nabla c \cdot \nabla v - \alpha \cdot c(1-c))v \, dV \\
&\text{The residue is now linear, so we calculate the Fréchet derivative to obtain a bilinear form:} \\
&a(c)(\delta,v) := \lim_{h \to 0} \frac{R(c+\delta \cdot h)(v) - R(c)(v)}{h} \\
&= \lim_{h \to 0} \frac{1}{h} \left[ \int_\Omega \partial_t c \cdot v \, dV + \int_\Omega h \cdot \partial_t\delta \cdot v \, dV + \int_\Omega D \cdot \nabla c \cdot \nabla v \, dV + \int_\Omega D \cdot \nabla \delta \cdot \nabla v \, dV - \int_\Omega \alpha \cdot (c + h \cdot \delta)(1-(c + h \cdot \delta)) \cdot v \, dV - R(c)(v) \right] \\
&= \lim_{h \to 0} \frac{1}{h} \left[ \int_\Omega \partial_t c \cdot v \, dV + \int_\Omega h \cdot \partial_t\delta \cdot v \, dV + \int_\Omega D \cdot \nabla c \cdot \nabla v \, dV + \int_\Omega h \cdot D \cdot \nabla \delta \cdot \nabla v \, dV - \int_\Omega -\alpha \cdot v(c-c^2-2 \cdot \delta \cdot h \cdot c +\delta \cdot h-\delta^2 \cdot h^2) \, dV - R(c)(v) \right] \\
&= \lim_{h \to 0} \frac{1}{h} \left[ \int_\Omega \partial_t c \cdot v \, dV + \int_\Omega h \cdot \partial_t\delta \cdot v \, dV + \int_\Omega D \cdot \nabla c \cdot \nabla v \, dV + \int_\Omega h \cdot D \cdot \nabla \delta \cdot \nabla v \, dV - \int_\Omega \alpha \cdot c \cdot v(1-c) \, dV - \int_\Omega -\alpha \cdot h \cdot \delta \cdot v(1-2c-h \cdot \delta) \, dV - R(c)(v) \right] \\
&\text{Since } \int_\Omega \partial_t c \cdot v \, dV + \int_\Omega D \cdot \nabla c \cdot \nabla v \, dV - \int_\Omega \alpha \cdot c \cdot v(1-c) \, dV = R(c)(v) \text{, this cancels with } -R(c)(v) \text{, simplifying the expression :} \\
&= \lim_{h \to 0} \frac{1}{h} \left[ \int_\Omega h \cdot \partial_t\delta \cdot v \, dV + \int_\Omega h \cdot D \cdot \nabla \delta \cdot \nabla v \, dV - \int_\Omega -\alpha \cdot h \cdot \delta \cdot v(1-2c-h \cdot \delta) \, dV \right] \\

&\text{Since } h \text{ in the denominator outside the integrals simplifies with those inside, and the term } -h \cdot \delta \text{ in the last integral goes to } 0 \text{ for } h \text{ going to } 0, \text{ the final result becomes:} \\
&a(c)(\delta,v)= \int_\Omega (\partial_t\delta \cdot v + D \cdot \nabla\delta \cdot \nabla v - \alpha \cdot v (1-2c) \cdot \delta) \, dV \\
&\text{Spatial and Temporal Discretization:} \\
&\text{We discretize the spatial domain using basis functions } N_i \text{ for } v \text{ and } N_j \text{ for }  \delta. \text{ We discretize the time domain using the backward Euler method.} \\
&\text{The discretized form of the equation becomes:} \\
&\sum_{i=1}^{n_{A}} \int_\Omega \left( N_{i} \frac{1}{\Delta t} N_{j} + \nabla N_{i} : D \cdot \nabla N_{j} - N_{i} \alpha (1 - 2 c) N_{j} \right) dV \\

\end{align*}
$$

