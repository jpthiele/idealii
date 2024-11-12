// ---------------------------------------------------------------------
//
// Copyright (C) 2022 - 2024 by the ideal.II authors
//
// This file is part of the ideal.II library.
//
// The ideal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of ideal.II.
//
// ---------------------------------------------------------------------

#ifndef INCLUDE_IDEAL_II_FE_FE_TOOLS_HH_
#define INCLUDE_IDEAL_II_FE_FE_TOOLS_HH_

#include <deal.II/lac/vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/fe/fe_tools.h>

#include <ideal.II/dofs/slab_dof_handler.hh>   
namespace idealii::slab::FETools{

    template <int dim,class VectorType>
    void interpolate_space(
        const slab::DoFHandler<dim> & dof1,
        const VectorType &            u1,
        const slab::DoFHandler<dim> &       dof2,
        VectorType & u2
    );

    template <int dim,class VectorType>
    void interpolate_space(
        const slab::DoFHandler<dim> & dof1,
        const VectorType &            u1,
        const slab::DoFHandler<dim> &       dof2,
        const dealii::AffineConstraints
            <typename VectorType::value_type> & constraints,
        VectorType & u2
    );


    template <int dim,class VectorType>
    void interpolate_time(
        const slab::DoFHandler<dim> & dof1,
        const VectorType &            u1,
        const slab::DoFHandler<dim> &       dof2,
        VectorType & u2
    );
    
    template <int dim,class VectorType>
    void interpolate_time(
        const slab::DoFHandler<dim> & dof1,
        const VectorType &            u1,
        const slab::DoFHandler<dim> &       dof2,
        const dealii::AffineConstraints
            <typename VectorType::value_type> & constraints,
        VectorType & u2
    );



    template <int dim,class VectorType>
    void interpolate(
        const slab::DoFHandler<dim> & dof1,
        const VectorType &            u1,
        const slab::DoFHandler<dim> &       dof2,
        VectorType & u2
    );
    
    template <int dim,class VectorType>
    void interpolate(
        const slab::DoFHandler<dim> & dof1,
        const VectorType &            u1,
        const slab::DoFHandler<dim> &       dof2,
        const dealii::AffineConstraints
            <typename VectorType::value_type> & constraints,
        VectorType & u2
    );
}

#endif //INCLUDE_IDEAL_II_FE_FE_TOOLS_HH_