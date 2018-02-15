/*
 * EMatrix_unit_tests.cpp
 *
 *  Created on: Feb 13, 2018
 *      Author: i-bird
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "DMatrix/EMatrix.hpp"
#include "memory/HeapMemory.hpp"

BOOST_AUTO_TEST_SUITE (EMatrix_test)

BOOST_AUTO_TEST_CASE( EMatrix_test_use)
{
	{
	EMatrixXd em;

	em.resize(8,5);

	for (size_t i = 0 ; i < 8 ; i++)
	{
		for (size_t j = 0 ; j < 5 ; j++)
		{em(i,j) = i*8+j;}
	}

	size_t pr = 0;
	em.packRequest(pr);

	// allocate the memory
	HeapMemory pmem;
	pmem.allocate(pr);
	ExtPreAlloc<HeapMemory> & mem = *(new ExtPreAlloc<HeapMemory>(pr,pmem));
	mem.incRef();

	BOOST_REQUIRE_EQUAL(pr,8*5*sizeof(double) + 2*sizeof(size_t));

	Pack_stat sts;
	em.pack(mem,sts);

	// Reset to zero
	for (size_t i = 0 ; i < 8 ; i++)
	{
		for (size_t j = 0 ; j < 5 ; j++)
		{em(i,j) = 0;}
	}

	Unpack_stat ps;
	em.unpack(mem,ps);

	for (size_t i = 0 ; i < 8 ; i++)
	{
		for (size_t j = 0 ; j < 5 ; j++)
		{BOOST_REQUIRE_EQUAL(em(i,j),i*8+j);}
	}
	}


	{
	EMatrix3d em;

	em.resize(3,3);

	for (size_t i = 0 ; i < 3 ; i++)
	{
		for (size_t j = 0 ; j < 3 ; j++)
		{em(i,j) = i*8+j;}
	}

	size_t pr = 0;
	em.packRequest(pr);

	// allocate the memory
	HeapMemory pmem;
	pmem.allocate(pr);
	ExtPreAlloc<HeapMemory> & mem = *(new ExtPreAlloc<HeapMemory>(pr,pmem));
	mem.incRef();

	BOOST_REQUIRE_EQUAL(pr,3*3*sizeof(double) + 2*sizeof(size_t));

	Pack_stat sts;
	em.pack(mem,sts);

	// Reset to zero
	for (size_t i = 0 ; i < 3 ; i++)
	{
		for (size_t j = 0 ; j < 3 ; j++)
		{em(i,j) = 0;}
	}

	Unpack_stat ps;
	em.unpack(mem,ps);

	for (size_t i = 0 ; i < 3 ; i++)
	{
		for (size_t j = 0 ; j < 3 ; j++)
		{BOOST_REQUIRE_EQUAL(em(i,j),i*8+j);}
	}
	}
}

BOOST_AUTO_TEST_SUITE_END()

