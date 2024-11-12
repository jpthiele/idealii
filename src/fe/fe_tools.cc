// ---------------------------------------------------------------------
//
// Copyright (C) 2022 - 2023 by the ideal.II authors
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

#include <deal.II/base/exceptions.h>
#include <deal.II/fe/fe_tools.h>
#include <ideal.II/fe/fe_tools.hh>
#include <ideal.II/numerics/vector_tools.hh>
#include <deal.II/lac/affine_constraints.h>

namespace idealii::slab::FETools
{
    
    template<int dim, class VectorType>
    void interpolate_time(
        const slab::DoFHandler<dim> & dof1,
        const VectorType &            u1,
        const slab::DoFHandler<dim> & dof2,
        VectorType &                  u2
    ){
        dealii::AffineConstraints<typename VectorType::value_type> dummy;
        dummy.close();
        interpolate_time(dof1, u1, dof2, dummy, u2);
    }

    template<int dim, class VectorType>
    void interpolate_time(
        const slab::DoFHandler<dim> & dof1,
        const VectorType &            u1,
        const slab::DoFHandler<dim> & dof2,
        const dealii::AffineConstraints
            <typename VectorType::value_type> & constraints,
        VectorType &                  u2
    ){
        Assert(u1.size()==dof1.n_dofs_spacetime(),
            dealii::ExcDimensionMismatch(u1.size(),dof1.n_dofs_spacetime()))
        Assert(u2.size()==dof2.n_dofs_spacetime(),
            dealii::ExcDimensionMismatch(u2.size(),dof2.n_dofs_spacetime()))
        
        //Todo: Assert space FEs are equal
        //Todo: More asserts?

        const auto n_dofs_t1 = dof1.dofs_per_cell_time();
        const auto n_dofs_t2 = dof2.dofs_per_cell_time();
        const auto n_dofs_x = dof1.n_dofs_space();
        const auto n_dofs_xt1 = n_dofs_t1*n_dofs_x;
        const auto n_dofs_xt2 = n_dofs_t2*n_dofs_x;

        dealii::FullMatrix<typename VectorType::value_type> ip;
        ip.reinit(n_dofs_t2,n_dofs_t1);

        dealii::FETools::get_interpolation_matrix(
            dof1.temporal()->get_fe(), 
            dof2.temporal()->get_fe(), 
            ip);

        u2=0;
        //Go over spatial DoFs
        unsigned int l = 0; //temporal elements
        for (unsigned int k = 0 ; k < n_dofs_x; ++k)
        {
            //loop over temporal elements
            for (const auto &cell_time :
                dof1.temporal()->active_cell_iterators())
            {
                l = cell_time->index();
                for (unsigned int j = 0 ; j < n_dofs_t2; ++j)
                for (unsigned int i = 0 ; i < n_dofs_t1; ++i)
                {
                    u2[k+j*n_dofs_x+l*n_dofs_xt2]+=
                        ip(j,i)*u1[k+i*n_dofs_x+l*n_dofs_xt1];
                }
            }
        }
        constraints.distribute(u2);
    }

    template<int dim, class VectorType>
    void interpolate_space(
        const slab::DoFHandler<dim> & dof1,
        const VectorType &            u1,
        const slab::DoFHandler<dim> & dof2,
        VectorType &                  u2
    ){
        dealii::AffineConstraints<typename VectorType::value_type> dummy;
        dummy.close();
        interpolate_space(dof1, u1, dof2, dummy, u2);
    }


    template<int dim, class VectorType>
    void interpolate_space(
        const slab::DoFHandler<dim> & dof1,
        const VectorType &            u1,
        const slab::DoFHandler<dim> & dof2,
        const dealii::AffineConstraints
            <typename VectorType::value_type> & constraints,
        VectorType &                  u2
    ){
        VectorType spatial1;
        VectorType spatial2;
        spatial1.reinit(dof1.n_dofs_space());
        spatial2.reinit(dof2.n_dofs_space());

        //loop over temporal elements
        unsigned int offset;
        for (const auto &cell_time :
            dof1.temporal()->active_cell_iterators())
        {
            offset = cell_time->index()*dof1.dofs_per_cell_time();
            for (unsigned int i = 0 ; i < dof1.dofs_per_cell_time();++i)
            {
                VectorTools::extract_subvector_at_time_dof(
                    u1,
                    spatial1,
                    i+offset
                );
                dealii::FETools::interpolate(
                    *dof1.spatial(),
                    spatial1,
                    *dof2.spatial(),
                    spatial2
                );
                for (unsigned int j = 0 ; j < dof2.n_dofs_space(); ++j){
                    u2[j+(offset+i)*dof2.n_dofs_space()]
                        = spatial2[j];
                }
            }
        }

        constraints.distribute(u2);
    }
}

#include "fe_tools.inst"

