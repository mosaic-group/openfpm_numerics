/*
 * PETSC_util.hpp
 *
 *  Created on: Jul 7, 2015
 *      Author: Pietro Incardona
 */

#ifndef PETSC_UTIL_HPP_
#define PETSC_UTIL_HPP_

#include <iostream>

#define PETSC_SAFE_CALL(call) {\
	PetscErrorCode err = call;\
	if (err != 0) {\
		std::string msg("Petsc error: ");\
		msg += std::string(__FILE__) + std::string(" ") + std::to_string(__LINE__);\
		PetscError(MPI_COMM_WORLD,__LINE__,__FUNCT__,__FILE__,err,PETSC_ERROR_INITIAL,"Error petsc");\
	}\
}


#endif /* MPI_UTIL_HPP_ */
