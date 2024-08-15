import numpy as np
from scipy.integrate import solve_ivp, simps

def lorentz(x,delta,lw,grad,stark,length):
    line_shape = 1/((np.pi * (lw + grad*length/2*stark)*(1 + ((delta + grad*stark*x)/(lw + grad*length/2*stark))**2)))
    return line_shape

def lorentz_abs_prof(x,delta,lw,OD,grad,stark,length):
    line_shape = 1/(2*np.pi) * (lw + grad*length/2*stark)/((delta + grad*stark*x)**2 + (1/2*(lw + grad*length/2*stark))**2)
    return OD*line_shape/max(1/(2*np.pi) * lw/(delta**2 + (1/2*lw)**2))

def rect_abs_prof(x,delta,lw,OD,grad,stark,length):
    line_shape = np.where(abs(delta + grad*stark*x) <= (lw + grad*length/2*stark), OD*lw/(lw + grad*length/2*stark), 0)
    return line_shape

def gauss_abs_prof(x,delta,lw,OD,grad,stark,length):
    line_shape = 1/((lw + grad*length/2*stark)*np.sqrt(2*np.pi)) * np.exp(-1/2 * (delta + grad*stark*x)**2/(lw + grad*length/2*stark)**2)
    return OD*line_shape/max(1/(lw*np.sqrt(2*np.pi)) * np.exp(-1/2 * delta**2/lw**2))

def grad_switch(t,t_switch):
    if t > t_switch:
        return 1
    else:
        return -1
    
def gauss(t,t_peak,amp,sigma):    
    line_shape = amp*np.exp(-(t - t_peak)**2 / (2 * sigma**2))
    return line_shape

def MBE_pol_z0(t, alpha, t0, A, sigma, delta, x, grad, stark, t_switch, g, gamma, w):

    E = gauss(t,t0,A,sigma)

    dalpha_dt = -(gamma/2)*alpha - (1j) * (delta + grad * stark * grad_switch(t_switch,t) * x)*alpha - (1j) * g * w * E

    return dalpha_dt

def interpolate_E_t(E,t_arr,tn):

    interp = np.interp(tn,t_arr,E)
    return interp

def interpolate_E_z(E,x_arr,xn):

    interp = np.interp(xn,x_arr,E)
    return interp

def MBE_pol_z(t, alpha, E, t_arr, delta, x_n, grad, stark, t_switch, g, gamma, w):
    
    Ez = interpolate_E_t(E,t_arr,t)

    dalpha_dt = -(gamma/2)*alpha - (1j) * (delta + grad * stark * grad_switch(t_switch,t) * x_n)*alpha + (1j) * g * w * Ez

    return dalpha_dt

def MBE_field(x, E, delta, t0, t1, alpha0, E0, t_arr, x_arr, grad, stark, t_switch, g, gamma, w0, OD, lw, L):

    # Calculate the optical field over dz
    # Since the the function calculates the spatial derivative of 
    # E and the eqn does not include any terms with E then then initial value of E is not needed
    # rhs_int is the initial value

    E1 = E0#-E0*OD*x

    # print(np.shape(E[:,0]))

    sol_alpha_z = solve_ivp(MBE_pol_z,
                            (t0, t1),
                            alpha0,
                            # args=(E[:,0], t_arr, delta, x, grad, stark, t_switch, g, gamma, w0),
                            args=(E1, t_arr, delta, x, grad, stark, t_switch, g, gamma, w0),
                            t_eval=t_arr,
                            dense_output=True,
                            vectorized=False,
                            method='RK23')
        
    alpha = sol_alpha_z.y.T

    dE_dx = -1 * OD * simps(lorentz(x_arr,delta,lw,grad,stark,L)*alpha,x=delta,axis=1)

    print(x)

    return dE_dx