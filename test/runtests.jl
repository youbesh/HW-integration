

module IntTest
	using HW_int
	using FactCheck

	facts("check demand") do
		@fact HW_int.q(1) --> 2
		@fact HW_int.q(4) --> 1
		@fact HW_int.q(9) --> 2/3
	end


	facts("check gauss adjustment") do
		@fact HW_int.ba2(4) --> 2
	end

	facts("eqm condition for Q2") do 
		@fact HW_int.dd(1,0,0) --> 0

	end
end 