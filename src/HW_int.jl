
module HW_int


	# question 1 b) 

	using FastGaussQuadrature
	using Roots
	using Sobol
	using PyPlot
	using Distributions

	# demand function
	q(p) = 2*p.^-0.5

	# gauss-legendre adjustment factors for map change
	ba2(ub) = (ub-0)/2
	ab(ub) = (0 + ub)/2

	# eqm condition for question 2
	# this is the equilibrium condition: total demand = supply, 
	# i.e. domestic + export demand = 2
	function dd(p,t1,t2)
	    exp(t1)/p .+ exp(t2) .* p^-0.5 - 2
	end
	# this is the integration function
	# it's the weighted sum from the slides.
	function intfun(px,g::Matrix,w::Vector)
		sum( w .* dd(px,g[:,1],g[:,2])) 
	end

	function question_1b(n)

		gl = gausslegendre(n);

		# bounds of integration
		# a = 0
		# b = pstar

		# equilibrium quantity
		# qstar = q(pstar)

		# CS1: p=4
		# ========

		pstar=4.0
		pqstar = pstar * q(pstar)
		figure()
		plot(gl[1],q(ba2(pstar).*gl[1] .+ ab(pstar)),"o")
		title("Gauss Legendre")

		CS1 = ba2(pstar) * (gl[2]' * q(ba2(pstar).*gl[1] .+ ab(pstar)))  - pqstar

		# CS2: p=1
		# ========
		pstar=1.0
		pqstar = pstar * q(pstar)
		CS2 = ba2(pstar) * (gl[2]' * q(ba2(pstar).*gl[1] .+ ab(pstar)))  - pqstar

		d = CS1-CS2

		println("estimated change in CS using $n gauss legendre nodes is $(round(d[1],5))")
		println("i.e. an error of $(round(abs(100*(2-d[1]))/2,5)) percent")
		println("")
	end

	# question 1 c)

	function question_1c(n)

		# CS1: p=4
		# ========

		pstar=4.0
		pqstar = pstar * q(pstar)

		pts = rand(n)*pstar
		figure()
		plot(pts,q(pts),"o")
		title("Monte Carlo")
		Integ = mean(q(pts)) 

		CS1 = Integ - pqstar

		# CS2: p=1
		# ========
		pstar=1.0
		pqstar = pstar * q(pstar)

		pts = rand(n)*pstar
		Integ = mean(q(pts)) 

		CS2 = Integ - pqstar

		d = CS1-CS2

		println("estimated change in CS using $n monte carlo nodes is $(round(d[1],5))")
		println("i.e. an error of $(round(abs(100*(2-d[1]))/2,5)) percent")
		println("")

	end

	function question_1d(n)

		# CS1: p=4
		# ========

		pstar=4.0
		pqstar = pstar * q(pstar)

		s = SobolSeq(1,[0],[pstar])  # 1-dimensional sobol sequence in [0,pstar]
		pts = zeros(n)
		for i in 1:n
			pts[i] = next(s)[1]
		end
		figure()
		plot(pts,q(pts),"o")
		title("Quasi Monte Carlo")
		Integ = mean(q(pts)) 

		CS1 = Integ - pqstar

		# CS2: p=1
		# ========
		pstar=1.0
		pqstar = pstar * q(pstar)

		s = SobolSeq(1,[0],[pstar])  # 1-dimensional sobol sequence
		pts = zeros(n)
		for i in 1:n
			pts[i] = next(s)[1]
		end
		pts = rand(n)*pstar
		Integ = mean(q(pts)) 

		CS2 = Integ - pqstar

		d = CS1-CS2

		println("estimated change in CS using $n quasi-monte carlo nodes is $(round(d[1],5))")
		println("i.e. an error of $(round(abs(100*(2-d[1]))/2,5)) percent")
		println("")

	end

	# question 2

	function question_2a(n)

		# expected price in equilibrium is such that expected demand equals supply.
		# pfun(p,theta1,theta2) = 
		# \int \int (theta1 1/p + theta2 * 1/sqrt(p) ) w(theta1,theta2) dthetea1 dtheta2 - 2
		# find root of pfun

		gh = gausshermite(n)

		Sigma = hcat([0.02, 0.01],[0.01,0.01])
		Omega = chol(Sigma,Val{:L})

		mu = [0.0;0.0]

		# kronecker product of grids and weights
		gr = hcat(kron(ones(n),gh[1]),kron(gh[1],ones(n)))
		wt = kron(gh[2],gh[2]) / sqrt(pi)	# watch out for the pi!

		# make adjustment for correlation in shocks
		grids = gr * Omega + zeros(n*n,2)   # zeros here would be a matrix with mu
		
		EP = fzero(x->intfun(x,grids,wt),0.01,10)

		VAR = fzero(x-> intfun(x^2,grids,wt) - EP^2,0.01,10) 

		return Dict("E[p]"=>EP, "Var[p]"=>VAR)

	end

	function question_2b(n)

		# expected price in equilibrium is such that expected demand equals supply.
		# pfun(p,theta1,theta2) = 
		# \int \int (theta1 1/p + theta2 * 1/sqrt(p) ) w(theta1,theta2) dthetea1 dtheta2 - 2
		# find root of pfun

		ln1 = LogNormal(0.0,0.02)
		ln2 = LogNormal(0.0,0.01)
		pt1 = rand(ln1,n)
		pt2 = rand(ln2,n)

		Sigma = hcat([0.02, 0.01],[0.01,0.01])
		Omega = chol(Sigma,Val{:L})

		mu = [0.0;0.0]

		# kronecker product of grids 
		gr = hcat(kron(ones(n),pt1),kron(pt2,ones(n)))
		wt = ones(n^2) ./ n^2

		# make adjustment for correlation in shocks
		grids = gr * Omega + zeros(n*n,2)   # zeros here would be a matrix with mu
		
		EP = fzero(x->mean(dd(x,grids,wt)),0.01,10)

		VAR = fzero(x-> mean(dd(x^2,grids,wt)) - EP^2,0.01,10) 

		return Dict("E[p]"=>EP, "Var[p]"=>VAR)

	end

	# function to run all questions
	function runall(n=10)
		println("running all questions of HW-integration:")
		println("results of question 1:")
		question_1b(n)	# make sure your function prints some kind of result!
		question_1c(n)
		question_1d(n)
		println("")
		println("results of question 2:")
		q2 = question_2a(n)
		println(q2)
		q2b = question_2b(n)
		println(q2b)
		println("end of HW-integration")
	end

end

