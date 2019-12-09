/*
 * FiniteDifferences.hpp
 *
 *  Created on: Sep 17, 2015
 *      Author: i-bird
 *  Modified on: Dec 09, 2019
 *      Author: amfoggia
 */

#ifndef OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_FDSCHEME_HPP_
#define OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_FDSCHEME_HPP_

#include "../Matrix/SparseMatrix.hpp"
#include "Grid/grid_dist_id.hpp"
#include "Grid/Iterators/grid_dist_id_iterator_sub.hpp"
#include "eq.hpp"
#include "NN/CellList/CellDecomposer.hpp"
#include "Grid/staggered_dist_grid_util.hpp"
#include "Grid/grid_dist_id.hpp"
#include "Vector/Vector_util.hpp"
#include "Grid/staggered_dist_grid.hpp"

/*! \brief Finite Differences
 *
 * This class is able to discretize on a Matrix any system of equations producing a linear system of type \f$Ax=b\f$. In order to create a consistent
 * Matrix it is required that each processor must contain a contiguous range on grid points without
 * holes. In order to ensure this, each processor produce a contiguous local labeling of its local
 * points. Each processor also add an offset equal to the number of local
 * points of the processors with id smaller than him, to produce a global and non overlapping
 * labeling. An example is shown in the figures down, here we have
 * a grid 8x6 divided across four processors each processor label locally its grid points
 *
 * \verbatim
 *
 +--------------------------+
 | 1   2   3   4| 1  2  3  4|
 |              |           |
 | 5   6   7   8| 5  6  7  8|
 |              |           |
 | 9  10  11  12| 9 10 11 12|
 +--------------------------+
 |13  14  15| 13 14 15 16 17|
 |          |               |
 |16  17  18| 18 19 20 21 22|
 |          |               |
 |19  20  21| 23 24 25 26 27|
 +--------------------------+

 *
 *
 * \endverbatim
 *
 * To the local relabelling is added an offset to make the local id global and non overlapping
 *
 *
 * \verbatim
 *
 +--------------------------+
 | 1   2   3   4|23 24 25 26|
 |              |           |
 | 5   6   7   8|27 28 29 30|
 |              |           |
 | 9  10  12  13|31 32 33 34|
 +--------------------------+
 |14  15  16| 35 36 37 38 39|
 |          |               |
 |17  18  19| 40 41 42 43 44|
 |          |               |
 |20  21  22| 45 46 47 48 49|
 +--------------------------+
 *
 *
 * \endverbatim
 *
 * \tparam Sys_eqs Definition of the system of equations
 *
 * # Examples
 *
 * ## Solve lid-driven cavity 2D for incompressible fluid (inertia=0 --> Re=0)
 *
 * In this case the system of equation to solve is
 *
 * \f$
 \left\{
 \begin{array}{c}
 \eta\nabla v_x + \partial_x P = 0 \quad Eq1 \\
 \eta\nabla v_y + \partial_y P = 0 \quad Eq2 \\
 \partial_x v_x + \partial_y v_y = 0 \quad Eq3
 \end{array}
 \right.  \f$

 and  boundary conditions

 * \f$
 \left\{
 \begin{array}{c}
 v_x = 0, v_y = 0 \quad x = 0 \quad B1\\
 v_x = 0, v_y = 1.0 \quad x = L \quad B2\\
 v_x = 0, v_y = 0 \quad y = 0 \quad B3\\
 v_x = 0, v_y = 0 \quad y = L \quad B4\\
 \end{array}
 \right.  \f$

 *
 * with \f$v_x\f$ and \f$v_y\f$ the velocity in x and y and \f$P\f$ Pressure
 *
 * In order to solve such system first we define the general properties of the system
 *
 *	\snippet eq_unit_test.hpp Definition of the system
 *
 * ## Define the equations of the system
 *
 * \snippet eq_unit_test.hpp Definition of the equation of the system in the bulk and at the boundary
 *
 * ## Define the domain and impose the equations
 *
 * \snippet eq_unit_test.hpp lid-driven cavity 2D
 *
 * # 3D
 *
 * A 3D case is given in the examples
 *
 */

template<typename Sys_eqs>
class FDScheme
{
public:

  //! Distributed grid map
  typedef grid_dist_id<Sys_eqs::dims,typename Sys_eqs::stype,aggregate<size_t>,typename Sys_eqs::b_grid::decomposition::extended_type> g_map_type;

  //! Type that specify the properties of the system of equations
  typedef Sys_eqs Sys_eqs_typ;

private:

  //! Encapsulation of the b term as constant
  struct constant_b
  {
    //! scalar
    typename Sys_eqs::stype scal;

    /*! \brief Constrictor from a scalar
     *
     * \param scal scalar
     *
     */
    constant_b(typename Sys_eqs::stype scal)
    {
      this->scal = scal;
    }

    /*! \brief Get the b term on a grid point
     *
     * \note It does not matter the grid point it is a scalar
     *
     * \param  key grid position (unused because it is a constant)
     *
     * \return the scalar
     *
     */
    typename Sys_eqs::stype get(grid_dist_key_dx<Sys_eqs::dims> & key)
    {
      return scal;
    }
  };

  //! Encapsulation of the b term as grid
  template<typename grid, unsigned int prp>
  struct grid_b
  {
    //! b term fo the grid
    grid & gr;

    /*! \brief gr grid that encapsulate
     *
     * \param gr grid
     *
     */
    grid_b(grid & gr)
      :gr(gr)
    {}

    /*! \brief Get the value of the b term on a grid point
     *
     * \param key grid point
     *
     * \return the value
     *
     */
    typename Sys_eqs::stype get(grid_dist_key_dx<Sys_eqs::dims> & key)
    {
      return gr.template get<prp>(key);
    }
  };

  //! Padding
  Padding<Sys_eqs::dims> pd;

  //! Sparse matrix triplet type
  typedef typename Sys_eqs::SparseMatrix_type::triplet_type triplet;

  //! Vector b
  typename Sys_eqs::Vector_type b;

  //! Domain Grid informations
  const grid_sm<Sys_eqs::dims,void> & gs;

  //! Get the grid spacing
  typename Sys_eqs::stype spacing[Sys_eqs::dims];

  //! mapping grid
  g_map_type g_map;

  //! row of the matrix
  size_t row;

  //! row on b
  size_t row_b;

  //! Grid points that has each processor
  openfpm::vector<size_t> pnt;

  //! Staggered position for each property
  comb<Sys_eqs::dims> s_pos[Sys_eqs::nvar];

  //! Each point in the grid has a global id, to decompose correctly the Matrix each processor contain a
  //! contiguos range of global id, example processor 0 can have from 0 to 234 and processor 1 from 235 to 512
  //! no processors can have holes in the sequence, this number indicate where the sequence start for this
  //! processor
  size_t s_pnt;

  /*! \brief Equation id + space position
   *
   */
  struct key_and_eq
  {
    //! space position
    grid_key_dx<Sys_eqs::dims> key;

    //! equation id
    size_t eq;
  };

  /*! \brief From the row Matrix position to the spatial position
   *
   * \param row Matrix row
   *
   * \return spatial position + equation id
   *
   */
  inline key_and_eq from_row_to_key(size_t row)
  {
    key_and_eq ke;
    auto it = g_map.getDomainIterator();

    while (it.isNext())
      {
	size_t row_low = g_map.template get<0>(it.get());

	if (row >= row_low * Sys_eqs::nvar && row < row_low * Sys_eqs::nvar + Sys_eqs::nvar)
	  {
	    ke.eq = row - row_low * Sys_eqs::nvar;
	    ke.key = g_map.getGKey(it.get());
	    return ke;
	  }

	++it;
      }
    std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " the row does not map to any position" << "\n";

    return ke;
  }

  /*! \brief calculate the mapping grid size with padding
   *
   * \param sz original grid size
   * \param pd padding
   *
   * \return padded grid size
   *
   */
  inline const std::vector<size_t> padding( const size_t (& sz)[Sys_eqs::dims], Padding<Sys_eqs::dims> & pd)
  {
    std::vector<size_t> g_sz_pad(Sys_eqs::dims);

    for (size_t i = 0 ; i < Sys_eqs::dims ; i++)
      g_sz_pad[i] = sz[i] + pd.getLow(i) + pd.getHigh(i);

    return g_sz_pad;
  }

  /*! \brief Check if the Matrix is consistent
   *
   */
  void consistency()
  {
    openfpm::vector<triplet> & trpl = A.getMatrixTriplets();

    // A and B must have the same rows
    if (row != row_b)
      std::cerr << "Error " << __FILE__ << ":" << __LINE__ << "the term B and the Matrix A for Ax=B must contain the same number of rows\n";

    // Indicate all the non zero rows
    openfpm::vector<unsigned char> nz_rows;
    nz_rows.resize(row_b);

    for (size_t i = 0 ; i < trpl.size() ; i++)
      nz_rows.get(trpl.get(i).row() - s_pnt*Sys_eqs::nvar) = true;

    // Indicate all the non zero colums
    // This check can be done only on single processor

    Vcluster<> & v_cl = create_vcluster();
    if (v_cl.getProcessingUnits() == 1)
      {
	openfpm::vector<unsigned> nz_cols;
	nz_cols.resize(row_b);

	for (size_t i = 0 ; i < trpl.size() ; i++)
	  nz_cols.get(trpl.get(i).col()) = true;

	// all the rows must have a non zero element
	for (size_t i = 0 ; i < nz_rows.size() ; i++)
	  {
	    if (nz_rows.get(i) == false)
	      {
		key_and_eq ke = from_row_to_key(i);
		std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " Ill posed matrix row " << i <<  " is not filled, position " << ke.key.to_string() << " equation: " << ke.eq << "\n";
	      }
	  }

	// all the colums must have a non zero element
	for (size_t i = 0 ; i < nz_cols.size() ; i++)
	  {
	    if (nz_cols.get(i) == false)
	      std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " Ill posed matrix colum " << i << " is not filled\n";
	  }
      }
  }

  /*! \brief Copy a given solution vector in a staggered grid
   *
   * \tparam Vct Vector type
   * \tparam Grid_dst target grid
   * \tparam pos set of properties
   *
   * \param v Vector
   * \param g_dst target staggered grid
   *
   */
  template<typename Vct, typename Grid_dst ,unsigned int ... pos> void copy_staggered(Vct & v, Grid_dst & g_dst)
  {
    // check that g_dst is staggered
    if (g_dst.is_staggered() == false)
      std::cerr << __FILE__ << ":" << __LINE__ << " The destination grid must be staggered " << std::endl;

#ifdef SE_CLASS1

    if (g_map.getLocalDomainSize() != g_dst.getLocalDomainSize())
      std::cerr << __FILE__ << ":" << __LINE__ << " The staggered and destination grid in size does not match " << std::endl;
#endif

    // sub-grid iterator over the grid map
    auto g_map_it = g_map.getDomainIterator();

    // Iterator over the destination grid
    auto g_dst_it = g_dst.getDomainIterator();

    while (g_map_it.isNext() == true)
      {
	typedef typename to_boost_vmpl<pos...>::type vid;
	typedef boost::mpl::size<vid> v_size;

	auto key_src = g_map_it.get();
	size_t lin_id = g_map.template get<0>(key_src);

	// destination point
	auto key_dst = g_dst_it.get();

	// Transform this id into an id for the Eigen vector

	copy_ele<Sys_eqs_typ, Grid_dst,Vct> cp(key_dst,g_dst,v,lin_id,g_map.size());

	boost::mpl::for_each_ref<boost::mpl::range_c<int,0,v_size::value>>(cp);

	++g_map_it;
	++g_dst_it;
      }
  }

  /*! \brief Copy a given solution vector in a normal grid
   *
   * \tparam Vct Vector type
   * \tparam Grid_dst target grid
   * \tparam pos set of property
   *
   * \param v Vector
   * \param g_dst target normal grid
   *
   */
  template<typename Vct, typename Grid_dst ,unsigned int ... pos> void copy_normal(Vct & v, Grid_dst & g_dst)
  {
    // check that g_dst is staggered
    if (g_dst.is_staggered() == true)
      std::cerr << __FILE__ << ":" << __LINE__ << " The destination grid must be normal " << std::endl;

    grid_key_dx<Grid_dst::dims> start;
    grid_key_dx<Grid_dst::dims> stop;

    for (size_t i = 0 ; i < Grid_dst::dims ; i++)
      {
	start.set_d(i,pd.getLow(i));
	stop.set_d(i,g_map.size(i) - pd.getHigh(i));
      }

    // sub-grid iterator over the grid map
    auto g_map_it = g_map.getSubDomainIterator(start,stop);

    // Iterator over the destination grid
    auto g_dst_it = g_dst.getDomainIterator();

    while (g_dst_it.isNext() == true)
      {
	typedef typename to_boost_vmpl<pos...>::type vid;
	typedef boost::mpl::size<vid> v_size;

	auto key_src = g_map_it.get();
	size_t lin_id = g_map.template get<0>(key_src);

	// destination point
	auto key_dst = g_dst_it.get();

	// Transform this id into an id for the vector

	copy_ele<Sys_eqs_typ, Grid_dst,Vct> cp(key_dst,g_dst,v,lin_id,g_map.size());

	boost::mpl::for_each_ref<boost::mpl::range_c<int,0,v_size::value>>(cp);

	++g_map_it;
	++g_dst_it;
      }
  }

  /*! \brief Impose an operator
   *
   * This function impose an operator on a particular grid region to produce the system
   *
   * Ax = b
   *
   * ## Stokes equation 2D, lid driven cavity with one splipping wall
   * \snippet eq_unit_test.hpp Copy the solution to grid
   *
   * \param op Operator to impose (A term)
   * \param num right hand side of the term (b term)
   * \param id Equation id in the system that we are imposing
   * \param it_d iterator that define where you want to impose
   *
   */
  template<typename T, typename bop, typename iterator> void impose_git(const T & op,
									bop num,
									long int id,
									const iterator & it_d)
  {
    openfpm::vector<triplet> & trpl = A.getMatrixTriplets();

    auto it = it_d;
    grid_sm<Sys_eqs::dims,void> gs = g_map.getGridInfoVoid();

    std::unordered_map<long int,float> cols;

    // iterate all the grid points
    while (it.isNext())
      {
	// get the position
	auto key = it.get();

	// Calculate the non-zero colums
	op.value(g_map,key,gs,spacing,cols,1.0);

	// indicate if the diagonal has been set
	bool is_diag = false;

	// create the triplet
	for ( auto it = cols.begin(); it != cols.end(); ++it )
	  {
	    trpl.add();
	    trpl.last().row() = g_map.template get<0>(key)*Sys_eqs::nvar + id;
	    trpl.last().col() = it->first;
	    trpl.last().value() = it->second;

	    if (trpl.last().row() == trpl.last().col())
	      is_diag = true;

	    //				std::cout << "(" << trpl.last().row() << "," << trpl.last().col() << "," << trpl.last().value() << ")" << "\n";
	  }

	// If does not have a diagonal entry put it to zero
	if (is_diag == false)
	  {
	    trpl.add();
	    trpl.last().row() = g_map.template get<0>(key)*Sys_eqs::nvar + id;
	    trpl.last().col() = g_map.template get<0>(key)*Sys_eqs::nvar + id;
	    trpl.last().value() = 0.0;
	  }

	b(g_map.template get<0>(key)*Sys_eqs::nvar + id) = num.get(key);

	cols.clear();

	// if SE_CLASS1 is defined check the position
#ifdef SE_CLASS1
	//			T::position(key,gs,s_pos);
#endif

	++row;
	++row_b;
	++it;
      }
  }

  /*! \brief Construct the gmap structure
   *
   */
  void construct_gmap()
  {
    Vcluster<> & v_cl = create_vcluster();

    // Calculate the size of the local domain
    size_t sz = g_map.getLocalDomainSize();

    // Get the total size of the local grids on each processors
    v_cl.allGather(sz,pnt);
    v_cl.execute();
    s_pnt = 0;

    // calculate the starting point for this processor
    for (size_t i = 0 ; i < v_cl.getProcessUnitID() ; i++)
      s_pnt += pnt.get(i);

    // resize b if needed
    b.resize(Sys_eqs::nvar * g_map.size(),Sys_eqs::nvar * sz);

    // Calculate the starting point

    // Counter
    size_t cnt = 0;

    // Create the re-mapping grid
    auto it = g_map.getDomainIterator();

    while (it.isNext())
      {
	auto key = it.get();

	g_map.template get<0>(key) = cnt + s_pnt;

	++cnt;
	++it;
      }

    // sync the ghost
    g_map.template ghost_get<0>();
  }

  /*! \initialize the object FDScheme
   *
   * \param domain simulation domain
   *
   */
  void Initialize(const Box<Sys_eqs::dims,typename Sys_eqs::stype> & domain)
  {
    construct_gmap();

    // Create a CellDecomposer and calculate the spacing

    size_t sz_g[Sys_eqs::dims];
    for (size_t i = 0 ; i < Sys_eqs::dims ; i++)
      {
	if (Sys_eqs::boundary[i] == NON_PERIODIC)
	  sz_g[i] = gs.getSize()[i] - 1;
	else
	  sz_g[i] = gs.getSize()[i];
      }

    CellDecomposer_sm<Sys_eqs::dims,typename Sys_eqs::stype> cd(domain,sz_g,0);

    for (size_t i = 0 ; i < Sys_eqs::dims ; i++)
      spacing[i] = cd.getCellBox().getHigh(i);
  }

public:

  /*! \brief set the staggered position for each property
   *
   * \param sp vector containing the staggered position for each property
   *
   */
  void setStagPos(comb<Sys_eqs::dims> (& sp)[Sys_eqs::nvar])
  {
    for (size_t i = 0 ; i < Sys_eqs::nvar ; i++)
      s_pos[i] = sp[i];
  }

  /*! \brief compute the staggered position for each property
   *
   * This is compute from the value_type stored by Sys_eqs::b_grid::value_type
   * the position of the staggered properties
   *
   *
   */
  void computeStag()
  {
    typedef typename Sys_eqs::b_grid::value_type prp_type;

    openfpm::vector<comb<Sys_eqs::dims>> c_prp[prp_type::max_prop];

    stag_set_position<Sys_eqs::dims,prp_type> ssp(c_prp);

    boost::mpl::for_each_ref< boost::mpl::range_c<int,0,prp_type::max_prop> >(ssp);
  }

  /*! \brief Get the specified padding
   *
   * \return the padding specified
   *
   */
  const Padding<Sys_eqs::dims> & getPadding()
  {
    return pd;
  }

  /*! \brief Return the map between the grid index position and the position in the distributed vector
   *
   * It is the map explained in the intro of the FDScheme
   *
   * \return the map
   *
   */
  const g_map_type & getMap()
  {
    return g_map;
  }

  /*! \brief Constructor
   *
   * The periodicity is given by the grid b_g
   *
   * \param stencil maximum extension of the stencil on each directions
   * \param domain the domain
   * \param b_g object grid that will store the solution
   *
   */
  FDScheme(const Ghost<Sys_eqs::dims,long int> & stencil,
	   const Box<Sys_eqs::dims,typename Sys_eqs::stype> & domain,
	   const typename Sys_eqs::b_grid & b_g)
    :pd({0,0,0},{0,0,0}),gs(b_g.getGridInfoVoid()),g_map(b_g,stencil,pd),row(0),row_b(0)
  {
    Initialize(domain);
  }

  /*! \brief Constructor
   *
   * The periodicity is given by the grid b_g
   *
   * \param pd Padding, how many points out of boundary are present
   * \param stencil maximum extension of the stencil on each directions
   * \param domain the domain
   * \param b_g object grid that will store the solution
   *
   */
  FDScheme(Padding<Sys_eqs::dims> & pd,
	   const Ghost<Sys_eqs::dims,long int> & stencil,
	   const Box<Sys_eqs::dims,typename Sys_eqs::stype> & domain,
	   const typename Sys_eqs::b_grid & b_g)
    :pd(pd),gs(b_g.getGridInfoVoid()),g_map(b_g,stencil,pd),row(0),row_b(0)
  {
    Initialize(domain);
  }

  /*! \brief Impose an operator
   *
   * This function impose an operator on a box region to produce the system
   *
   * Ax = b
   *
   * ## Stokes equation, lid driven cavity with one splipping wall
   * \snippet eq_unit_test.hpp lid-driven cavity 2D
   *
   * \param op Operator to impose (A term)
   * \param num right hand side of the term (b term)
   * \param id Equation id in the system that we are imposing
   * \param start starting point of the box
   * \param stop stop point of the box
   * \param skip_first skip the first point
   *
   */
  template<typename T> void impose(const T & op,
				   typename Sys_eqs::stype num,
				   long int id,
				   const long int (& start)[Sys_eqs::dims],
				   const long int (& stop)[Sys_eqs::dims],
				   bool skip_first = false)
  {
    grid_key_dx<Sys_eqs::dims> start_k;
    grid_key_dx<Sys_eqs::dims> stop_k;

    bool increment = false;
    if (skip_first == true)
      {
	start_k = grid_key_dx<Sys_eqs::dims>(start);
	stop_k = grid_key_dx<Sys_eqs::dims>(start);

	auto it = g_map.getSubDomainIterator(start_k,stop_k);

	if (it.isNext() == true)
	  increment = true;
      }

    // add padding to start and stop
    start_k = grid_key_dx<Sys_eqs::dims>(start);
    stop_k = grid_key_dx<Sys_eqs::dims>(stop);

    auto it = g_map.getSubDomainIterator(start_k,stop_k);

    if (increment == true)
      ++it;

    constant_b b(num);

    impose_git(op,b,id,it);

  }

  /*! \brief In case we want to impose a new b re-using FDScheme we have to call
   *         This function
   */
  void new_b()
  {row_b = 0;}

  /*! \brief In case we want to impose a new A re-using FDScheme we have to call
   *         This function
   *
   */
  void new_A()
  {row = 0;}

  //! type of the sparse matrix
  typename Sys_eqs::SparseMatrix_type A;

  /*! \brief produce the Matrix
   *
   *  \return the Sparse matrix produced
   *
   */
  typename Sys_eqs::SparseMatrix_type & getA()
  {
#ifdef SE_CLASS1
    consistency();
#endif
    A.resize(g_map.size()*Sys_eqs::nvar,g_map.size()*Sys_eqs::nvar,
	     g_map.getLocalDomainSize()*Sys_eqs::nvar,
	     g_map.getLocalDomainSize()*Sys_eqs::nvar);

    return A;

  }

  /*! \brief produce the B vector
   *
   *  \return the vector produced
   *
   */
  typename Sys_eqs::Vector_type & getB()
  {
#ifdef SE_CLASS1
    consistency();
#endif

    return b;
  }

  /*! \brief Copy the vector into the grid
   *
   * ## Copy the solution into the grid
   * \snippet eq_unit_test.hpp Copy the solution to grid
   *
   * \tparam scheme Discretization scheme
   * \tparam Grid_dst type of the target grid
   * \tparam pos target properties
   *
   * \param v Vector that contain the solution of the system
   * \param start point
   * \param stop point
   * \param g_dst Destination grid
   *
   */
  template<unsigned int ... pos, typename Vct, typename Grid_dst> void copy(Vct & v,const long int (& start)[Sys_eqs_typ::dims], const long int (& stop)[Sys_eqs_typ::dims], Grid_dst & g_dst)
  {
    if (is_grid_staggered<Sys_eqs>::value())
      {
	if (g_dst.is_staggered() == true)
	  copy_staggered<Vct,Grid_dst,pos...>(v,g_dst);
	else
	  {
	    // Create a temporal staggered grid and copy the data there
	    auto & g_map = this->getMap();

	    // Convert the ghost in grid units

	    Ghost<Grid_dst::dims,long int> g_int;
	    for (size_t i = 0 ; i < Grid_dst::dims ; i++)
	      {
		g_int.setLow(i,g_map.getDecomposition().getGhost().getLow(i) / g_map.spacing(i));
		g_int.setHigh(i,g_map.getDecomposition().getGhost().getHigh(i) / g_map.spacing(i));
	      }

	    staggered_grid_dist<Grid_dst::dims,typename Grid_dst::stype,typename Grid_dst::value_type,typename Grid_dst::decomposition::extended_type, typename Grid_dst::memory_type, typename Grid_dst::device_grid_type> stg(g_dst,g_int,this->getPadding());
	    stg.setDefaultStagPosition();
	    copy_staggered<Vct,decltype(stg),pos...>(v,stg);

	    // sync the ghost and interpolate to the normal grid
	    stg.template ghost_get<pos...>();
	    stg.template to_normal<Grid_dst,pos...>(g_dst,this->getPadding(),start,stop);
	  }
      }
    else
      {
	copy_normal<Vct,Grid_dst,pos...>(v,g_dst);
      }
  }

  /*! \brief Copy the vector into the grid
   *
   * ## Copy the solution into the grid
   * \snippet eq_unit_test.hpp Copy the solution to grid
   *
   * \tparam scheme Discretization scheme
   * \tparam Grid_dst type of the target grid
   * \tparam pos target properties
   *
   * \param v Vector that contain the solution of the system
   * \param g_dst Destination grid
   *
   */
  template<unsigned int ... pos, typename Vct, typename Grid_dst> void copy(Vct & v, Grid_dst & g_dst)
  {
    long int start[Sys_eqs_typ::dims];
    long int stop[Sys_eqs_typ::dims];

    for (size_t i = 0 ; i < Sys_eqs_typ::dims ; i++)
      {
	start[i] = 0;
	stop[i] = g_dst.size(i);
      }

    this->copy<pos...>(v,start,stop,g_dst);
  }
};

#define EQS_FIELDS 0
#define EQS_SPACE 1


#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_FDSCHEME_HPP_ */
