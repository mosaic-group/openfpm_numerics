//
// Created by jstark on 01.06.21.
//

#ifndef OPENFPM_NUMERICS_BOUNDARYCONDITIONS_METHODOFIMAGES_HPP
#define OPENFPM_NUMERICS_BOUNDARYCONDITIONS_METHODOFIMAGES_HPP

#include <cmath>
// Include OpenFPM header files
#include "Vector/vector_dist_subset.hpp"

#define PID_VECTOR_TYPE openfpm::vector<aggregate<int>>
#define KEY_VECTOR_TYPE openfpm::vector<vect_dist_key_dx>


/**@brief Class for getting mirror particles to impose Neumann BCs
 *
 * @details Source refers to the particles close to the boundary that will we mirrored at the boundary. The source
 * particles do not form an individual subset because subsets mustn't overlap and source particles belong to the real
 * particles. Instead, store the keys of the desired source particles in a vector and pass them to the constructor.
 *
 * @tparam SurfaceNormal Index of property to which the surface normal should be written to.
 * @tparam vd_type Template type of input particle vector_dist.
 */
template <size_t SurfaceNormal, typename vd_type>
class MethodOfImages {

public:
	typedef vector_dist_subset<vd_type::dims, typename vd_type::stype, typename vd_type::value_type> vd_subset_type;
	typedef Point<vd_type::dims, typename vd_type::stype> point_type;

	/**@brief Constructor
	 *
	 * @param vd Input particle vector_dist of type vd_type.
	 * @param keys_source Vector containing keys of source particles (key-type: openfpm::vector<vect_dist_key_dx>).
	 * @param subset_id_mirror ID of subset containing the mirror particles (default=1).
	 */
	MethodOfImages(
			vd_type & vd,
			const KEY_VECTOR_TYPE & keys_source,
			const size_t subset_id_real = 0,
			const size_t subset_id_mirror = 1)
			: keys_source(keys_source)
			, subset_id_real(subset_id_real)
			, subset_id_mirror(subset_id_mirror)
			, Real(vd, subset_id_real)
			, Mirror(vd, subset_id_mirror)

	{
		#ifdef SE_CLASS1
		check_if_ghost_isometric(vd);
		#endif // SE_CLASS1
	}

	//	Member variables
	size_t subset_id_real; ///< ID of subset containing the real particles (default=0).
	size_t subset_id_mirror; ///< ID of subset containing the mirror particles (default=1).
	KEY_VECTOR_TYPE keys_source; ///< Vector containing keys of source particles.
	PID_VECTOR_TYPE pid_mirror; ///< Vector containing indices of mirror particles.
	vd_subset_type Mirror; ///< Subset containing the mirror particles.
	vd_subset_type Real;

	/**@brief Place mirror particles along the surface normal.
	 *
	 * @param vd Input particle vector_dist of type vd_type.
	 */
	void get_mirror_particles(vd_type & vd)
	{
		std::cout << "Ghost size before placing mirror particles = "
				<< vd.size_local_with_ghost() - vd.size_local() << std::endl;

		for (int i = 0; i < keys_source.size(); i++)
		{
			auto key        = keys_source.get(i);
			point_type xp   = vd.getPos(key);
			point_type n    = vd.template getProp<SurfaceNormal>(key);
			
			point_type xm   = xp + 2 * n;

			#ifdef SE_CLASS1
			if(!point_lies_on_this_processor(vd, xm))
			{
				std::cerr << __FILE__ << ":" << __LINE__ << " Error: Ghost layer is too small. Source and mirror"
															" particles that belong together must lie on the same"
															" processor."
															" Create a bigger ghost layer which is bigger than the"
															" mirror particle layer. Aborting..." << std::endl;
				abort();
			}
			#endif // SE_CLASS1
			
			vd.add();
			for (size_t d = 0; d < vd_type::dims; d++)
			{
				vd.getLastPos()[d] = xm[d];
			}
			vd.getLastSubset(subset_id_mirror);
		}
// No vd.map() here, because we want to keep the source and the mirror particles on the same node
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		vd.template ghost_get();
		std::cout << "Ghost size after placing mirror particles = "
				<< vd.size_local_with_ghost() - vd.size_local() << std::endl;
		Mirror.update();
		pid_mirror = Mirror.getIds();
		#ifdef SE_CLASS1
		check_size_mirror_source_equal();
		#endif // SE_CLASS1
	}

	
	
	/**@brief Copies the values stored in PropToMirror from each source particle to its respective mirror particles
	 *
	 * @tparam PropToMirror Index of property storing the values that should be mirrored.
	 * @param vd Input particle vector_dist of type vd_type.
	 */
	template <size_t PropToMirror>
	void apply_reflection(vd_type & vd)
	{
		check_size_mirror_source_equal();
		vd.template ghost_get<PropToMirror>(KEEP_PROPERTIES); // Update Ghost layer.
		for (int i = 0; i < keys_source.size(); ++i)
		{
			auto key_source = keys_source.get<0>(i); // Get key of one source particle
			auto key_mirror = pid_mirror.get<0>(i); // Get key of corresponding mirror particle to that source particle
			vd.template getProp<PropToMirror>(key_mirror) = vd.template getProp<PropToMirror>(key_source);
		}
	}


private:
	void check_if_ghost_isometric(vd_type & vd)
	{
		for(int d=0; d<vd_type::dims; d++)
		{
			if (vd.getDecomposition().getGhost().getLow(0) != vd.getDecomposition().getGhost().getLow(d))
			{
				std::cerr << __FILE__ << ":" << __LINE__ << "Ghost layer doesn't have the same size in all dimensions"
															". Use an isometric ghost layer. Aborting..."
															<< std::endl;
				abort();
			}
		}
	}

	/**@brief Checks if the size of the ghost layer is bigger or equal the size of the mirror layer. This is needed s
	 * .t. added mirror particles can be accessed by the same processor on which the corresponding source lies.
	 * @param vd Input particle vector_dist of type vd_type.
	 * @param p Point for which we want to know if it is covered by the current processor (incl. its ghost).
	 * @return True, if processor has access to point. False, if point lays too far away.
	 */
	bool point_lies_on_this_processor(vd_type & vd, point_type p)
	{
		double g_width = fabs(vd.getDecomposition().getGhost().getLow(0));
		Ghost<vd_type::dims, typename vd_type::stype> g(g_width);
		
		auto & subs = vd.getDecomposition().getSubDomains();
		
		bool is_inside = false;
		
		for (int i = 0 ; i < subs.size() ; i++)
		{
			SpaceBox<vd_type::dims, typename vd_type::stype> sub = subs.get(i);
			sub.enlarge(g);
			is_inside |= sub.isInside(p);
		}

		if (!is_inside)  {std::cout << "Processor does not have the point" << std::endl;}

		return is_inside;
	}

	/**@brief Checks if local vector containing source particle ids and vector containing mirror particle ids match
	 * in size. Necessary, because in apply_mirror, source ids and mirror ids are iterated in same loop on same
	 * processor.
	 *
	 */
	void check_size_mirror_source_equal()
	{
		if (pid_mirror.size() != keys_source.size())
		{
			std::cout << "pid_mirror.size() = " << pid_mirror.size() << ", keys_source.size() = " << keys_source.size()
					<< std::endl;
			std::cerr << __FILE__ << ":" << __LINE__
					<< " Error: Local vector of source-IDs has different size than local vector of mirror-IDs. Matching "
					   "source and mirror particle IDs must be stored on same processor." << std::endl;
			abort();
		}
	}



};




#endif //OPENFPM_NUMERICS_BOUNDARYCONDITIONS_METHODOFIMAGES_HPP
