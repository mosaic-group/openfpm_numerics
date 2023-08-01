//
// Created by jstark on 28.05.21.
//

#ifndef OPENFPM_PDATA_TESTHELPERS_HPP
#define OPENFPM_PDATA_TESTHELPERS_HPP

// modified, based on https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
template <typename T>
bool AlmostEqualRelativeAndAbs(T A, T B,
                               T maxDiff, T maxRelDiff = std::numeric_limits<T>::epsilon())
{
	// Check if the numbers are really close -- needed
	// when comparing numbers near zero.
	T diff = fabs(A - B);
	if (diff <= maxDiff)
		return true;
	
	A = fabs(A);
	B = fabs(B);
	T largest = (B > A) ? B : A;
	
	if (diff <= largest * maxRelDiff)
		return true;
	return false;
}


template <size_t Field, typename grid_type_1, typename grid_type_2, std::string Field_type>
bool are_equal_grid(grid_type_1 & grid_1, grid_type_2 & grid_2)
{
	bool equal = true;
	assert(grid_type_1::dims == grid_type_2::dims);
	assert(grid_1.size() == grid_2.size());

	for (int d = 0; d < grid_dest_type::dims; ++d)
	{
		auto dom_1 = grid_1.getDomainIterator();
		auto dom_2 = grid_2.getDomainIterator();
		while (dom_sc.isNext())
		{
			auto key_1 = dom_1.get();
			auto key_2 = dom_2.get();
			
			auto v1 = grid_1.template get<attr_ds>(key_1)[d];
			auto v2 = grid_2.template get<attr_sc>(key_2)[d];
			
			if (!AlmostEqualRelativeAndAbs(v1, v2, std::numeric_limits<Field_type>::epsilon()))
			{
				equal = false;
				break;
			}
			
			++dom_1;
			++dom_2;
		}
	}
	return equal;
}





#endif //OPENFPM_PDATA_TESTHELPERS_HPP
