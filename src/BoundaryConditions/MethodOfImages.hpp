//
// Created by jstark on 01.06.21.
//

#ifndef OPENFPM_NUMERICS_BOUNDARYCONDITIONS_METHODOFIMAGES_HPP
#define OPENFPM_NUMERICS_BOUNDARYCONDITIONS_METHODOFIMAGES_HPP

#include <cmath>

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
	typedef Point<vd_type::dims, double> point_type;
	
	typedef PID_VECTOR_TYPE openfpm::vector<aggregate<int>>
	typedef KEY_VECTOR_TYPE openfpm::vector<vect_dist_key_dx>
	
	/**@brief Constructor
	 *
	 * @param vd Input particle vector_dist of type vd_type.
	 * @param keys_source Vector containing keys of source particles (key-type: openfpm::vector<vect_dist_key_dx>).
	 * @param subset_id_mirror ID of subset containing the mirror particles (default=1).
	 */
	MethodOfImages(
			vd_type & vd,
			const KEY_VECTOR_TYPE & keys_source,
			const size_t subset_id_mirror = 1)
			: keys_source(keys_source)
			, subset_id_mirror(subset_id_mirror)
			, Mirror(vd, subset_id_mirror)
	
	{
		check_ghost_thick_enough(vd);
	}
	
	//	Member variables
	size_t subset_id_mirror; ///< ID of subset containing the mirror particles (default=2).
	KEY_VECTOR_TYPE keys_source; ///< Vector containing keys of source particles.
	PID_VECTOR_TYPE pid_mirror; ///< Vector containing indices of mirror particles.
	vd_subset_type Mirror; ///< Subset containing the mirror particles.
	
	/**@brief Place mirror particles along the surface normal.
	 *
	 * @param vd Input particle vector_dist of type vd_type.
	 */
	void get_mirror_particles(vd_type & vd)
	{
		for (int i = 0; i < keys_source.size(); i++)
		{
			auto key        = keys_source.get(i);
			point_type xp   = vd.getPos(key);
			point_type n    = vd.template getProp<SurfaceNormal>(key);
			double distance = n.norm() * 2.0;
			
			point_type xm   = xp + n * distance;
			
			vd.add();
			for (size_t d = 0; d < vd_type::dims; d++)
			{
				vd.getLastPos()[d] = xm[d];
			}
			vd.getLastSubset(subset_id_mirror);
		}
// No vd.map() here, because we want to keep the source and the mirror particles on the same node
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		Mirror.update();
		pid_mirror = Mirror.getIds();
		check_size_mirror_source_equal();
	}
	
	/**@brief Copies the values stored in PropToMirror from each source particle to its respective mirror particles
	 *
	 * @tparam PropToMirror Index of property storing the values that should be mirrored.
	 * @param vd Input particle vector_dist of type vd_type.
	 */
	template <size_t PropToMirror>
	void apply_reflection(vd_type & vd)
	{
		vd.template ghost_get<PropToMirror>(KEEP_PROPERTIES); // Update Ghost layer.
		for (int i = 0; i < keys_source.size(); ++i)
		{
			auto key_source = keys_source.get<0>(i); // Get key of one source particle
			auto key_mirror = pid_mirror.get<0>(i); // Get key of corresponding mirror particle to that source particle
			vd.template getProp<PropToMirror>(key_mirror) = vd.template getProp<PropToMirror>(key_source);
		}
	}

	
private:
	/**@brief Checks if the size of the ghost layer is bigger or equal the size of the mirror layer. This is needed s
	 * .t. added mirror particles can be accessed by the same processor on which the corresponding source lies.
	 *
	 * @param vd
	 */
	void check_ghost_thick_enough(vd_type & vd)
	{
		double mirror_thickness = _b_up - _b_low;
		for (size_t i = 0 ; i < vd_type::dims ; i++)
		{
			if (fabs(vd.getDecomposition().getGhost().getLow(i)) < mirror_thickness)
			{
				std::cerr << __FILE__ << ":" << __LINE__ << " Error: thickness of the mirror layer (" <<
						mirror_thickness <<	") is bigger than the ghost layer on the dimension " << i << " which is "
						<< fabs(vd.getDecomposition().getGhost().getLow(i)) << ". Create a bigger ghost layer." << std::endl;
				abort();
			}
		}
	}
	
	/**@brief Checks if local vector containing source particle ids and vector containing mirror particle ids match
	 * in size. Necessary, because in apply_mirror, source ids and mirror ids are iterated in same loop on same
	 * processor.
	 *
	 */
	void check_size_mirror_source_equal()
	{
		std::cout << "pid_mirror.size() = " << pid_mirror.size() << ", keys_source.size() = " << keys_source.size()
				<< std::endl;
		if (pid_mirror.size() != keys_source.size())
		{
			std::cerr << __FILE__ << ":" << __LINE__
					<< " Error: Local vector of source-IDs has different size than local vector of mirror-IDs. Matching "
					   "source and mirror particle IDs must be stored on same processor." << std::endl;
			abort();
		}
	}
	
	
	
};




#endif //OPENFPM_NUMERICS_BOUNDARYCONDITIONS_METHODOFIMAGES_HPP
