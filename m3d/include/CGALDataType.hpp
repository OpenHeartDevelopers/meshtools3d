#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_constant_domain_field_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>
#include "My_Labeled_image_mesh_domain_3.hpp"
#include "MyImageWrapper.hpp"

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;

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



// Personalised data types for domain with manual segmentation
// Domain
typedef CGAL::My_Labeled_image_mesh_domain_3<K> Mesh_domain_manualseg;
// Triangulation
#ifdef CGAL_LINKED_WITH_TBB
  typedef CGAL::Mesh_triangulation_3< Mesh_domain_manualseg,
    CGAL::Kernel_traits<Mesh_domain_manualseg>::Kernel, // Same as sequential
    CGAL::Parallel_tag                        // Tag to activate parallelism
  >::type Tr_manualseg;
#else
  typedef CGAL::Mesh_triangulation_3<Mesh_domain_manualseg>::type Tr_manualseg;
#endif
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr_manualseg> C3t3_manualseg;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr_manualseg> Mesh_criteria_manualseg;
typedef CGAL::Mesh_facet_criteria_3<Tr_manualseg> Facet_criteria_manualseg;
typedef CGAL::Mesh_cell_criteria_3<Tr_manualseg> Cell_criteria_manualseg;

typedef C3t3_manualseg::Cells_in_complex_iterator Cell_iterator_manualseg;

struct MeshingParams {
    FaceNumericalType facet_angle         = 30;
    FaceNumericalType facet_size          = 0.8;
    FaceNumericalType facet_distance      = 4;
    CellNumericalType cell_rad_edge_ratio = 2.0;
    CellNumericalType cell_size           = 1.0;
    double            rescale_factor      = 1.0;
};

// Traits bundles for the templated meshing dispatcher in MeshingPipeline.hpp.
// Each traits struct groups the CGAL types that vary between the labeled-image
// branch and the manual-segmentation branch of main.cpp.
struct LabeledMeshingTraits {
    using Domain        = Mesh_domain;
    using C3t3          = ::C3t3;
    using FacetCriteria = Facet_criteria;
    using CellCriteria  = Cell_criteria;
    using MeshCriteria  = Mesh_criteria;
};

struct ManualSegMeshingTraits {
    using Domain        = Mesh_domain_manualseg;
    using C3t3          = C3t3_manualseg;
    using FacetCriteria = Facet_criteria_manualseg;
    using CellCriteria  = Cell_criteria_manualseg;
    using MeshCriteria  = Mesh_criteria_manualseg;
};










