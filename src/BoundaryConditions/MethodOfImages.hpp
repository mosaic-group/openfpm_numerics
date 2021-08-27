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
	openfpm::vector<openfpm::vector<size_t>> key_map_source_mirror;
	
	/**@brief Place mirror particles along the surface normal.
	 *
	 * @param vd Input particle vector_dist of type vd_type.
	 */
	void get_mirror_particles(vd_type & vd)
	{
		auto lastReal = vd.size_local();
		for (int i = 0; i < keys_source.size(); i++)
		{
			auto key_source  = keys_source.get(i);
			
			size_t id_source = key_source.getKey();
			size_t id_mirror = lastReal + i;
			
			point_type xp   = vd.getPos(key_source);
			point_type n    = vd.template getProp<SurfaceNormal>(key_source);
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
			openfpm::vector<size_t> pair = {id_source, id_mirror};
			key_map_source_mirror.add(pair);
		}
		vd.map(); // We redistribute the mirror particles to the processors they belong to. Source and mirror will be
		// still seen by the same processor via the ghost layer.
		vd.template ghost_get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
	void apply_noflux(vd_type & vd)
	{
		vd.template ghost_get<PropToMirror>(KEEP_PROPERTIES); // Update Ghost layer.
		
		for(int i = 0; i < key_map_source_mirror.size(); ++i)
		{
			openfpm::vector<size_t> row = key_map_source_mirror.get(i);
			vect_dist_key_dx key_source, key_mirror;
			key_source.setKey(row.get(0));
			key_mirror.setKey(row.get(1));
			vd.template getProp<PropToMirror>(key_mirror) = vd.template getProp<PropToMirror>(key_source);
		}
	}


private:
	/**@brief Checks if the ghost layer has the same size in all dimensions. This is required to ensure that the
	 * ghost layer is bigger than the mirror layer in all dimensions. This is needed s.t. added mirror particles can
	 * be accessed by the same processor on which the corresponding source lies.
	 *
	 * @param vd Input particle vector_dist of type vd_type.
	 */
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
#ifdef SE_CLASS1
	/**@brief Checks if a point, i.e. the added mirror particle, lies on the same processor as the source particle, i
	 * .e. the ghost is big enough to cover the space where the mirror particle was placed. Mirror and Source must
	 * lie on same processor in order to avoid communication overhead when mirror (= copy from source to mirror
	 * particle) is applied.
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

	/**@brief Checks if local vector containing source particle ids and mirror particle subset match in size, i.e.
	 * for each source particle a mirror particle exists and vice versa.
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
#endif // SE_CLASS1

};




#endif //OPENFPM_NUMERICS_BOUNDARYCONDITIONS_METHODOFIMAGES_HPP
