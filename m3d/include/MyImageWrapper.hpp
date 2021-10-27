//  -*- c++ -*-
//  Wrapper for INR image                                    Oct/21/2016
//
//
//  This library  is distributed in the  hope that it will  be useful, but
//  WITHOUT   ANY  WARRANTY;   without  even   the  implied   warranty  of
//  MERCHANTABILITY  or FITNESS  FOR A  PARTICULAR PURPOSE.   See  the GNU
//  Lesser General Public License for more details.
//
//  Wrapper for INR image: this library is used to wrap INRreader class into a
//  suitable function for CGAL domain definition; This Class is used to override
//  the label numbering during mesh generation step
//
//  Developer: Cesare Corrado cesare.corrado@kcl.ac.uk
//==========================================================================



#ifndef MY_IMAGE_WRAPPER_H
#define MY_IMAGE_WRAPPER_H

#include "INRreader.hpp"
#include <CGAL/function_objects.h>

/*
namespace CGAL {

namespace Mesh_3 {
*/


/**
 * @class MyImageWrapper
 *
 * Wraps a labeled image into a labeled function which takes his values into
 * N. Uses trilinear interpolation.
 * Note: Image has to be labeled with unsigned char
 * returns only 0 and 1
 */
template<class BGT,
         typename Image_word_type = unsigned char,
         typename Return_type = int,
         typename Transform = CGAL::Identity<Image_word_type> >
class MyImageWrapper
{
public:
  // Types
  typedef Return_type return_type;
  typedef Image_word_type word_type;
  typedef typename BGT::Point_3   Point_3;

  /// Constructor
  MyImageWrapper(const INRreader& image, 
                                const Transform& transform = Transform(),
                                const Image_word_type value_outside = 0)
    : r_im_(image)
    , transform(transform)
    , value_outside(value_outside)
  {
    CGAL_assertion(transform(value_outside) == return_type());
  }

  // Default copy constructor and assignment operator are ok

  /// Destructor
  ~MyImageWrapper() {}

  /**
   * Returns an int corresponding to the label at point \c p
   * @param p the input point
   * @return the label at point \c p
   */
  return_type operator()(const Point_3& p, const bool = true) const
  {

    bool isin = r_im_.isPointInsideSegmentation(p.x(),p.y(),p.z());
    //bool isin = r_im_.isPointInsideSegmentationTrilinear(p.x(),p.y(),p.z());
    if(isin)
    {
        return 1;
    }
    else
    {
        return 0;
    }
  }

private:
  /// Labeled image to wrap
  const INRreader & r_im_;
  const Transform transform;
  const Image_word_type value_outside;

};  // end class Image_to_labeled_function_wrapper

/*
}  // end namespace Mesh_3

}  // end namespace CGAL
*/
#endif // CGAL_MESH_3_IMAGE_TO_LABELED_FUNCTION_WRAPPER_H
