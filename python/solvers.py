from ngsolve.la import InnerProduct
from math import sqrt


def CG(mat, rhs, pre=None, sol=None, precision=1e-12, maxsteps = 100, printrates = True, initialize = True):
    """preconditioned conjugate gradient method"""

    u = sol if sol else rhs.CreateVector()
    d = rhs.CreateVector()
    w = rhs.CreateVector()
    s = rhs.CreateVector()

    if initialize: u[:] = 0.0
    d.data = rhs - mat * u
    w.data = pre * d if pre else d
    err0 = sqrt(InnerProduct(d,w))
    s.data = w
    wdn = InnerProduct (w,d)
    
    for it in range(maxsteps):
        w.data = mat * s
        wd = wdn
        as_s = InnerProduct (s, w)
        alpha = wd / as_s
        u.data += alpha * s
        d.data += (-alpha) * w

        w.data = pre*d if pre else d
        
        wdn = InnerProduct (w, d)
        beta = wdn / wd

        s *= beta
        s.data += w

        err = sqrt(wd)
        if err < precision*err0: break
            
        if printrates:
            print ("it = ", it, " err = ", err)

    return u




def QMR(matrix, rhs, fdofs, pre1=None, pre2=None, sol=None, maxsteps = 100, printrates = True, initialize = True, ep = 1.0, tol = 1e-7):
	
	u = sol if sol else rhs.CreateVector()
	
	r = rhs.CreateVector()
	v = rhs.CreateVector()	
	v_tld = rhs.CreateVector()
	w = rhs.CreateVector()
	w_tld = rhs.CreateVector()
	y = rhs.CreateVector()
	y_tld = rhs.CreateVector()
	z = rhs.CreateVector()
	z_tld = rhs.CreateVector()	
	p = rhs.CreateVector()
	p_tld = rhs.CreateVector()
	q = rhs.CreateVector()	
	d = rhs.CreateVector()
	s = rhs.CreateVector()

	if (initialize): u[:] = 0.0

	r.data = rhs - matrix * u
	v_tld.data = r
	y.data = pre1 * v_tld if pre1 else v_tld
	
	rho = InnerProduct(y,y)
	rho = sqrt(rho)
	
	w_tld.data = r 
	z.data = pre2.T * w_tld if pre2 else w_tld
	
	xi = InnerProduct(z,z)
	xi = sqrt(xi)
	
	gamma = 1.0
	eta = -1.0
	theta = 0.0

	
	for i in range(1,maxsteps+1):
		if (rho == 0.0):
			print('Breakdown in rho')
			return
		if (xi == 0.0):
			print('Breakdown in xi')
			return
		v.data = (1.0/rho) * v_tld
		y.data = (1.0/rho) * y

		w.data = (1.0/xi) * w_tld
		z.data = (1.0/xi) * z

		delta = InnerProduct(z,y)
		if (delta == 0.0):
			print('Breakdown in delta')
			return
		
		y_tld.data = pre2 * y if pre2 else y
		z_tld.data = pre1.T * z if pre1 else z

		if (i > 1):
			p.data = (-xi*delta / ep) * p
			p.data += y_tld

			q.data = (-rho * delta / ep) * q
			q.data += z_tld
		else:
			p.data = y_tld
			q.data = z_tld
		
		p_tld.data = matrix * p
		ep = InnerProduct(q, p_tld)
		if (ep == 0.0):
			print('Breakdown in epsilon')
			return

		beta = ep/delta
		if (beta == 0.0):
			print('Breakdown in beta')
			return
		
		v_tld.data = p_tld - beta * v;		

		y.data = pre1 * v_tld if pre1 else v_tld
		
		rho_1 = rho
		rho = InnerProduct(y,y)
		rho = sqrt(rho)		
		
		
		w_tld.data = matrix.T * q
		w_tld.data -= beta * w
		
		z.data = pre2.T * w_tld if pre2 else w_tld		
		
		xi = InnerProduct(z,z)
		xi = sqrt(xi)

		gamma_1 = gamma
		theta_1 = theta

		theta = rho/(gamma_1 * abs(beta)) 
		gamma = 1.0 / sqrt(1.0 + theta * theta)
		if (gamma == 0.0):
			print('Breakdown in gamma')
			return
	
		eta = -eta * rho_1 * gamma * gamma / (beta * gamma_1 * gamma_1);		

		if (i > 1):
			d.data = (theta_1 * theta_1 * gamma * gamma) * d
			d.data += eta * p

			s.data = (theta_1 * theta_1 * gamma * gamma) * s	
			s.data += eta * p_tld

		else: 
			d.data = eta * p
			s.data = eta * p_tld
		
	    
		u.data += d
		r.data -= s

		#Projected residuum: Better terminating condition necessary?
		ResNorm = sqrt( np.dot(r.FV().NumPy()[fdofs],r.FV().NumPy()[fdofs]))
		#ResNorm = sqrt(InnerProduct(r,r))
		
		if (ResNorm <= tol):
			break
				
		if (printrates):
			print ("it = ", i, " Residuennorm = ", ResNorm)




#Compare: Iterative Krylov methods for large linear systems - Henk A. van der Vorst
def MinRes(matrix, rhs, pre=None, sol=None, maxsteps = 100, printrates = True, initialize = True, tol = 1e-7):

	u = sol if sol else rhs.CreateVector()

	v_new = rhs.CreateVector()
	v = rhs.CreateVector()	
	v_old = rhs.CreateVector()
	w = rhs.CreateVector()
	w_old = rhs.CreateVector()
	w_2old = rhs.CreateVector()
	hv = rhs.CreateVector()
	
	if (initialize):
		u[:] = 0.0
		hv.data = rhs
	else:
		hv.data = rhs - matrix * u
		
	v.data = pre * hv if pre else hv
	
	#First Step
	beta = sqrt(InnerProduct(v,v))
	ResNorm = beta
	ResNorm_old = beta
	
	if (printrates):
		print("it = ", 0, " Residuennorm = ", ResNorm)		
		
	eta = beta
	gamma = 1
	gamma_old = 1
	sigma = 0
	sigma_old = 0
	
	j = 1
	while (j < maxsteps+1 and ResNorm > tol):

		hv.data = 1.0 / beta * v
		v.data = hv

		v_new.data = matrix * v
		hv.data = pre * v_new if pre else v_new
		
		alpha = InnerProduct(v,hv)
			
		v_new.data = hv - alpha * v - beta * v_old	
		
		beta_new = sqrt(InnerProduct(v_new, v_new))

		delta = gamma * alpha - gamma_old * sigma * beta
		rho1 = sqrt(delta * delta + beta_new * beta_new)
		rho2 = sigma * alpha + gamma_old * gamma * beta
		rho3 = sigma_old * beta
	
		gamma_new = delta / rho1
		sigma_new = beta_new / rho1

		w.data = 1.0 / rho1 * (v - rho3 * w_2old - rho2 * w_old)
		u.data += gamma_new * eta * w

		tmp = -sigma_new * eta
		eta = tmp

		#update of residuum
		ResNorm = abs(sigma_new) * ResNorm_old
		
		print("it = ", j, " Residuennorm = ", ResNorm)		
	
		j += 1

		#updates
		beta = beta_new
		gamma_old = gamma
		gamma = gamma_new
		sigma_old = sigma
		sigma = sigma_new
		ResNorm_old = ResNorm
		
		v_old.data = v
		v.data = v_new

		w_2old.data = w_old
		w_old.data = w


