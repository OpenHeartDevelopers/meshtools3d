//  -*- c++ -*-
//  Personalised labeled mesh domain for images    Oct/21/2016
//
//
//  This library  is distributed in the  hope that it will  be useful, but
//  WITHOUT   ANY  WARRANTY;   without  even   the  implied   warranty  of
//  MERCHANTABILITY  or FITNESS  FOR A  PARTICULAR PURPOSE.   See  the GNU
//  Lesser General Public License for more details.
//
//  To Use MyWrapper 
//
//  Developer: Cesare Corrado cesare.corrado@kcl.ac.uk
//==========================================================================

#ifndef CGAL_MY_LABELED_IMAGE_MESH_DOMAIN_3_H
#define CGAL_MY_LABELED_IMAGE_MESH_DOMAIN_3_H


#include <CGAL/Random.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Bbox_3.h>
#include "MyImageWrapper.hpp"
#include "INRreader.hpp"

namespace CGAL {

/**
 * @class My_Labeled_image_mesh_domain_3
 *
 *
 */
template<class BGT,
         typename Image_word_type = unsigned char,typename Subdomain_index = int,
         class Wrapper = MyImageWrapper<BGT, Image_word_type, Subdomain_index> >

class My_Labeled_image_mesh_domain_3
: public Labeled_mesh_domain_3<Wrapper, BGT>
{
public:
  typedef Labeled_mesh_domain_3<Wrapper, BGT> Base;

  typedef typename Base::Sphere_3 Sphere_3;
  typedef typename Base::FT FT;
  typedef BGT Geom_traits;
  typedef CGAL::Bbox_3 Bbox_3;

  /// Constructor
  My_Labeled_image_mesh_domain_3(const INRreader& image,
                              const FT& error_bound = FT(1e-3),
                              CGAL::Random* p_rng = NULL)
    :Base(Wrapper(image),compute_bounding_box(image), error_bound, p_rng)
  {}

  My_Labeled_image_mesh_domain_3(const INRreader& image,
                              const CGAL::Bbox_3& bbox,
                              const FT& error_bound = FT(1e-3),
                              CGAL::Random* p_rng = NULL)
    :Base(Wrapper(image), bbox, error_bound, p_rng)
  {}

  /// Destructor
  virtual ~My_Labeled_image_mesh_domain_3() {}

  using Base::bbox;

private:
  /// Returns a box enclosing image \c im
  Bbox_3 compute_bounding_box(const INRreader& im) const
  {
    return Bbox_3(-1,-1,-1,
                  im.xdim()*im.vx()+1, im.ydim()*im.vy()+1, im.zdim()*im.vz()+1);
  }

private:
  // Disabled copy constructor & assignment operator
  typedef My_Labeled_image_mesh_domain_3<BGT> Self;
  My_Labeled_image_mesh_domain_3(const Self& src);
  Self& operator=(const Self& src);

};  // end class My_Labeled_image_mesh_domain_3



}  // end namespace CGAL



#endif // CGAL_LABELED_IMAGE_MESH_DOMAIN_3_H
