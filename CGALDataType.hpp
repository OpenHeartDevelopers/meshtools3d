#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Labeled_image_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>


// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Labeled_image_mesh_domain_3<CGAL::Image_3,K> Mesh_domain;

// Triangulation
#ifdef CGAL_LINKED_WITH_TBB
  typedef CGAL::Mesh_triangulation_3< Mesh_domain,
    CGAL::Kernel_traits<Mesh_domain>::Kernel, // Same as sequential
    CGAL::Parallel_tag                        // Tag to activate parallelism
  >::type Tr;
#else
  typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
#endif
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef CGAL::Mesh_facet_criteria_3<Tr> Facet_criteria;
typedef CGAL::Mesh_cell_criteria_3<Tr> Cell_criteria;

typedef C3t3::Cells_in_complex_iterator Cell_iterator;

typedef Tr::Geom_traits::FT 	FaceNumericalType;
typedef Tr::Geom_traits::FT 	CellNumericalType; // in the guide was Tr::FT but compiler says does not name a type
