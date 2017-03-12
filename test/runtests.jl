LOAD_PATH
push!(LOAD_PATH, "C:/Users/Youssef/Desktop/M2_semester_ 2/Numerical_Methods/HW-integration/src")

module IntTest
	using HW_int
	using Base.Test

	@test 1==1
	@test_approx_eq HW_int.question_1b(100) 4
	@test_approx_eq_eps HW_int.question_1c(10000) 4 0.05
	@test_approx_eq_eps HW_int.question_1d(100) 4 0.01

	@test_approx_eq_eps HW_int.question_2a(100)[1] HW_int.question_2b(100)[1] 0.05
	@test_approx_eq_eps HW_int.question_2a(100)[2] HW_int.question_2b(100)[1] 0.3




end
