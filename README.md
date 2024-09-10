# FruitFlies
Gradient echo memory simulations based on the  Maxwell-Bloch equations

The ability to simulate a gradient echo memory is essential for studying
the protocol dynamics and optimising the performance criteria, such as efficiency.
The GEM protocol can be described by the Maxwell Bloch equations
(MBEs),


$\frac{\partial}{\partial t}\alpha(t,z,\Delta) = - [\frac{\gamma}{2} + i(\Delta + \Delta_{ext})]\alpha(t,z,\Delta) - iwgE(t,z)$,

$\frac{\partial}{\partial z}E(t,z) = i\int^{\infty}_{- \infty} N(\Delta)\alpha(t,z,\Delta)d\Delta$,

where $\alpha$ is the atomic polarisation, $\gamma$ is the decay rate from the excited state, $\Delta_{ext}$ is the detuning gradient, $w$ is the atomic excitation, $g$ is the atomic transition coupling strength, $E$ is the slowly varying envelope of the electric field, $N$ is the atomic density, $t$ is time, $z$ is position and $\Delta$ is detuning from resonance. 