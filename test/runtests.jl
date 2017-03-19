LOAD_PATH
push!(LOAD_PATH, "C:/Users/Youssef/Desktop/M2_semester_ 2/Numerical_Methods/HW-integration/src")

module IntTest
	using HW_int
	using Base.Test

	@test 1==1
	# The last two tests are working. The first three were working when I had done
	# the questions without the plots. I can make them work by ending each function
	# in HW_int by an object that groups into it the integration result and the plot,
	# of the form "respective_method = [Estimation Result; Plot]". But then the
	# individual functions would not display the plots "automatically" as asked
	# in the exercise. One would have to display the plot by running
	# "function_name(n)[2]". So I thought i'd keep the nice automatic display of
	# plots, and see from the correction how one would proceed with testing, since
	# testing was optional anyway.
	@test_approx_eq HW_int.question_1b(100) 4
	@test_approx_eq_eps HW_int.question_1c(10000) 4 0.05
	@test_approx_eq_eps HW_int.question_1d(100) 4 0.01


	@test_approx_eq_eps HW_int.question_2a(100)[1] HW_int.question_2b(100)[1] 0.05
	@test_approx_eq_eps HW_int.question_2a(100)[2] HW_int.question_2b(100)[1] 0.3




end
