/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
 */


#include <aspect/simulator.h>
#include <aspect/free_surface.h>
#include <aspect/global.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1_eulerian.h>

#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/numerics/vector_tools.h>
//#include <deal.II/numerics/matrix_tools.h>
//#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>


//from step-26
#include <deal.II/base/utilities.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/matrix_tools.h>

namespace aspect
{
  namespace Assemblers
  {
    template <int dim>
    void
    ApplyStabilization<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>       &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim>      &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>& > (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>& > (data_base);

      if (!this->get_parameters().free_surface_enabled)
        return;

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();

      const typename DoFHandler<dim>::active_cell_iterator cell (&this->get_triangulation(),
                                                                 scratch.finite_element_values.get_cell()->level(),
                                                                 scratch.finite_element_values.get_cell()->index(),
                                                                 &this->get_dof_handler());

      const unsigned int n_face_q_points = scratch.face_finite_element_values.n_quadrature_points;
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();

      // only apply on free surface faces
      if (cell->at_boundary() && cell->is_locally_owned())
        for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
          if (cell->face(face_no)->at_boundary())
            {
              const types::boundary_id boundary_indicator
                = cell->face(face_no)->boundary_id();

              if (this->get_parameters().free_surface_boundary_indicators.find(boundary_indicator)
                  == this->get_parameters().free_surface_boundary_indicators.end())
                continue;

              scratch.face_finite_element_values.reinit(cell, face_no);

              this->compute_material_model_input_values (this->get_solution(),
                                                         scratch.face_finite_element_values,
                                                         cell,
                                                         false,
                                                         scratch.face_material_model_inputs);

              this->get_material_model().evaluate(scratch.face_material_model_inputs, scratch.face_material_model_outputs);

              for (unsigned int q_point = 0; q_point < n_face_q_points; ++q_point)
                {
                  for (unsigned int i = 0, i_stokes = 0; i_stokes < stokes_dofs_per_cell; /*increment at end of loop*/)
                    {
                      if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                        {
                          scratch.phi_u[i_stokes] = scratch.face_finite_element_values[introspection.extractors.velocities].value(i, q_point);
                          ++i_stokes;
                        }
                      ++i;
                    }

                  const Tensor<1,dim>
                  gravity = this->get_gravity_model().gravity_vector(scratch.face_finite_element_values.quadrature_point(q_point));
                  double g_norm = gravity.norm();

                  // construct the relevant vectors
                  const Tensor<1,dim> n_hat = scratch.face_finite_element_values.normal_vector(q_point);
                  const Tensor<1,dim> g_hat = (g_norm == 0.0 ? Tensor<1,dim>() : gravity/g_norm);

                  const double pressure_perturbation = scratch.face_material_model_outputs.densities[q_point] *
                                                       this->get_timestep() *
                                                       this->get_free_surface_handler().get_stabilization_term() *
                                                       g_norm;

                  // see Kaus et al 2010 for details of the stabilization term
                  for (unsigned int i=0; i< stokes_dofs_per_cell; ++i)
                    for (unsigned int j=0; j< stokes_dofs_per_cell; ++j)
                      {
                        // The fictive stabilization stress is (phi_u[i].g)*(phi_u[j].n)
                        const double stress_value = -pressure_perturbation*
                                                    (scratch.phi_u[i]*g_hat) * (scratch.phi_u[j]*n_hat)
                                                    *scratch.face_finite_element_values.JxW(q_point);

                        data.local_matrix(i,j) += stress_value;
                      }
                }
            }
    }
  }

  template <int dim>
  FreeSurfaceHandler<dim>::FreeSurfaceHandler (Simulator<dim> &simulator,
                                               ParameterHandler &prm)
    : sim(simulator),  // reference to the simulator that owns the FreeSurfaceHandler
      free_surface_fe (FE_Q<dim>(1),dim), // Q1 elements which describe the mesh geometry
      free_surface_dof_handler (sim.triangulation)
  {
    parse_parameters(prm);
  }

  template <int dim>
  FreeSurfaceHandler<dim>::~FreeSurfaceHandler ()
  {
    // Free the Simulator's mapping object, otherwise
    // when the FreeSurfaceHandler gets destroyed,
    // the mapping's reference to the mesh displacement
    // vector will be invalid.
    sim.mapping.reset();
  }



  template <int dim>
  void FreeSurfaceHandler<dim>::set_assemblers()
  {
    aspect::Assemblers::ApplyStabilization<dim> *surface_stabilization
      = new aspect::Assemblers::ApplyStabilization<dim>();

    sim.assemblers->stokes_system.push_back(
      std_cxx11::unique_ptr<aspect::Assemblers::ApplyStabilization<dim> > (surface_stabilization));

    // Note that we do not want face_material_model_data, because we do not
    // connect to a face assembler. We instead connect to a normal assembler,
    // and compute our own material_model_inputs in apply_stabilization
    // (because we want to use the solution instead of the current_linearization_point
    // to compute the material properties).
    sim.assemblers->stokes_system_assembler_on_boundary_face_properties.needed_update_flags |= (update_values  |
        update_gradients |
        update_quadrature_points |
        update_normal_vectors |
        update_JxW_values);
  }

  template <int dim>
  void FreeSurfaceHandler<dim>::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection ("Free surface");
    {
      prm.declare_entry("Free surface stabilization theta", "0.5",
                        Patterns::Double(0,1),
                        "Theta parameter described in Kaus et. al. 2010. "
                        "An unstabilized free surface can overshoot its "
                        "equilibrium position quite easily and generate "
                        "unphysical results.  One solution is to use a "
                        "quasi-implicit correction term to the forces near the "
                        "free surface.  This parameter describes how much "
                        "the free surface is stabilized with this term, "
                        "where zero is no stabilization, and one is fully "
                        "implicit.");
      prm.declare_entry("Surface velocity projection", "normal",
                        Patterns::Selection("normal|vertical"),
                        "After each time step the free surface must be "
                        "advected in the direction of the velocity field. "
                        "Mass conservation requires that the mesh velocity "
                        "is in the normal direction of the surface. However, "
                        "for steep topography or large curvature, advection "
                        "in the normal direction can become ill-conditioned, "
                        "and instabilities in the mesh can form. Projection "
                        "of the mesh velocity onto the local vertical direction "
                        "can preserve the mesh quality better, but at the "
                        "cost of slightly poorer mass conservation of the "
                        "domain.");
      prm.declare_entry ("Additional tangential mesh velocity boundary indicators", "",
                         Patterns::List (Patterns::Anything()),
                         "A comma separated list of names denoting those boundaries "
                         "where there the mesh is allowed to move tangential to the "
                         "boundary. All tangential mesh movements along "
                         "those boundaries that have tangential material velocity "
                         "boundary conditions are allowed by default, this parameters "
                         "allows to generate mesh movements along other boundaries that are "
                         "open, or have prescribed material velocities or tractions."
                         "\n\n"
                         "The names of the boundaries listed here can either be "
                         "numbers (in which case they correspond to the numerical "
                         "boundary indicators assigned by the geometry object), or they "
                         "can correspond to any of the symbolic names the geometry object "
                         "may have provided for each part of the boundary. You may want "
                         "to compare this with the documentation of the geometry model you "
                         "use in your model.");
      prm.declare_entry("Diffusivity", "1e-4",
                        Patterns::Double(0,1),
                        "Diffusivity constant");						 
    }
    prm.leave_subsection ();
  }

  template <int dim>
  void FreeSurfaceHandler<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection ("Free surface");
    {
      free_surface_theta = prm.get_double("Free surface stabilization theta");
	  std::string advection_dir = prm.get("Surface velocity projection");
	  diffusivity = prm.get_double("Diffusivity");
      if ( advection_dir == "normal")
        advection_direction = SurfaceAdvection::normal;
      else if ( advection_dir == "vertical")
        advection_direction = SurfaceAdvection::vertical;
      else
        AssertThrow(false, ExcMessage("The surface velocity projection must be ``normal'' or ``vertical''."));


      // Create the list of tangential mesh movement boundary indicators
      try
        {
          const std::vector<types::boundary_id> x_additional_tangential_mesh_boundary_indicators
            = sim.geometry_model->translate_symbolic_boundary_names_to_ids(Utilities::split_string_list
                                                                           (prm.get ("Additional tangential mesh velocity boundary indicators")));

          tangential_mesh_boundary_indicators = sim.parameters.tangential_velocity_boundary_indicators;
          tangential_mesh_boundary_indicators.insert(x_additional_tangential_mesh_boundary_indicators.begin(),
                                                     x_additional_tangential_mesh_boundary_indicators.end());
        }
      catch (const std::string &error)
        {
          AssertThrow (false, ExcMessage ("While parsing the entry <Free surface/Additional tangential "
                                          "mesh velocity boundary indicators>, there was an error. Specifically, "
                                          "the conversion function complained as follows: "
                                          + error));
        }
    }
    prm.leave_subsection ();
  }



  template <int dim>
  void FreeSurfaceHandler<dim>::execute()
  {
    if (!sim.parameters.free_surface_enabled)
      return;
    sim.computing_timer.enter_section("Free surface");


    //Make the constraints for the elliptic problem.  On the free surface, we
    //constrain mesh velocity to be v.n, on free slip it is constrainted to
    //be tangential, and on no slip boundaries it is zero.
    // std::cout<<"Entered execute()\n";
		make_constraints();
		// std::cout<<"make_constraints complete\n";


    // Assemble and solve the vector Laplace problem which determines
    // the mesh displacements in the interior of the domain
    compute_mesh_displacements();

    // Interpolate the mesh velocity into the same
    // finite element space as used in the Stokes solve, which
    // is needed for the ALE corrections.
    interpolate_mesh_velocity();

    // After changing the mesh we need to rebuild things
    sim.rebuild_stokes_matrix = sim.rebuild_stokes_preconditioner = true;
    sim.computing_timer.exit_section("Free surface");
  }



  template <int dim>
  void FreeSurfaceHandler<dim>::make_constraints()
  {
    if (!sim.parameters.free_surface_enabled)
      return;

    // Now construct the mesh displacement constraints
    mesh_displacement_constraints.clear();
    mesh_displacement_constraints.reinit(mesh_locally_relevant);

    // mesh_displacement_constraints can use the same hanging node
    // information that was used for mesh_vertex constraints.
    mesh_displacement_constraints.merge(mesh_vertex_constraints);

    // Add the vanilla periodic boundary constraints
    typedef std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int> > periodic_boundary_pairs;
    periodic_boundary_pairs pbp = sim.geometry_model->get_periodic_boundary_pairs();
    for (periodic_boundary_pairs::iterator p = pbp.begin(); p != pbp.end(); ++p)
      DoFTools::make_periodicity_constraints(free_surface_dof_handler, (*p).first.first, (*p).first.second, (*p).second, mesh_displacement_constraints);

    //Zero out the displacement for the zero-velocity boundary indicators
    for (std::set<types::boundary_id>::const_iterator p = sim.parameters.zero_velocity_boundary_indicators.begin();
         p != sim.parameters.zero_velocity_boundary_indicators.end(); ++p)
      VectorTools::interpolate_boundary_values (free_surface_dof_handler, *p,
                                                ZeroFunction<dim>(dim), mesh_displacement_constraints);

    // Zero out the displacement for the prescribed velocity boundaries
    // if the boundary is not in the set of tangential mesh boundaries
    for (std::map<types::boundary_id, std::pair<std::string, std::string> >::const_iterator p = sim.parameters.prescribed_velocity_boundary_indicators.begin();
         p != sim.parameters.prescribed_velocity_boundary_indicators.end(); ++p)
		{
			if (tangential_mesh_boundary_indicators.find(p->first) == tangential_mesh_boundary_indicators.end())
			  {
				VectorTools::interpolate_boundary_values (free_surface_dof_handler, p->first,
														  ZeroFunction<dim>(dim), mesh_displacement_constraints);
			  } 
		  
		  
			// std::vector<bool> mask(sim.introspection.component_masks.velocities.size(), false);
			// std::cout<<"Length of mask = "<<mask.size()<<std::endl;
			// const std::string &comp = sim.parameters.prescribed_velocity_boundary_indicators[p->first].first;

			// if (comp.length()>0)
			// {
			  // for (std::string::const_iterator direction=comp.begin(); direction!=comp.end(); ++direction)
				// {
				  // switch (*direction)
					// {
					  // case 'x':
						// mask[sim.introspection.component_indices.velocities[0]] = true;
						// break;
					  // case 'y':
						// mask[sim.introspection.component_indices.velocities[1]] = true;
						// break;
					  // case 'z':
						// // we must be in 3d, or 'z' should never have gotten through
						// Assert (dim==3, ExcInternalError());
						// if (dim==3)
						  // mask[sim.introspection.component_indices.velocities[2]] = true;
						// break;
					  // default:
						// Assert (false, ExcInternalError());
					// }
				// }
			// }
			// else
			// {
			  // // no mask given -- take all velocities
			  // for (unsigned int i=0; i<sim.introspection.component_masks.velocities.size(); ++i)
				// mask[i]=sim.introspection.component_masks.velocities[i];
			// }

			// std::cout<<"Length of mask = "<<mask.size()<<std::endl;
			// VectorTools::interpolate_boundary_values (free_surface_dof_handler,//sim.dof_handler,
														// p->first,
														// ZeroFunction<dim>(dim),
														// mesh_displacement_constraints,
														// mask); 
		  
		}
	  

    // Zero out the displacement for the traction boundaries
    // if the boundary is not in the set of tangential mesh boundaries
    for (std::map<types::boundary_id, std::pair<std::string, std::string> >::const_iterator p = sim.parameters.prescribed_traction_boundary_indicators.begin();
         p != sim.parameters.prescribed_traction_boundary_indicators.end(); ++p)
      {
        if (tangential_mesh_boundary_indicators.find(p->first) == tangential_mesh_boundary_indicators.end())
          {
            VectorTools::interpolate_boundary_values (free_surface_dof_handler, p->first,
                                                      ZeroFunction<dim>(dim), mesh_displacement_constraints);
          }
      }

    sim.signals.pre_compute_no_normal_flux_constraints(sim.triangulation);
    // Make the no flux boundary constraints for boundaries with tangential mesh boundaries
    VectorTools::compute_no_normal_flux_constraints (free_surface_dof_handler,
                                                     /* first_vector_component= */
                                                     0,
                                                     tangential_mesh_boundary_indicators,
                                                     mesh_displacement_constraints, *sim.mapping);

    // make the periodic boundary indicators no displacement normal to the boundary
    std::set< types::boundary_id > periodic_boundaries;
    for (periodic_boundary_pairs::iterator p = pbp.begin(); p != pbp.end(); ++p)
      {
        periodic_boundaries.insert((*p).first.first);
        periodic_boundaries.insert((*p).first.second);
      }
    VectorTools::compute_no_normal_flux_constraints (free_surface_dof_handler,
                                                     /* first_vector_component= */
                                                     0,
                                                     periodic_boundaries,
                                                     mesh_displacement_constraints, *sim.mapping);
    sim.signals.post_compute_no_normal_flux_constraints(sim.triangulation);

    // For the free surface indicators we constrain the displacement to be v.n
    LinearAlgebra::Vector boundary_velocity;
    boundary_velocity.reinit(mesh_locally_owned, mesh_locally_relevant, sim.mpi_communicator);
    
		//Project velocity from stokes solution onto boundary
		project_velocity_onto_boundary( boundary_velocity );
		//Apply hillslope diffusion
		if (diffusivity) diffuse_surface(boundary_velocity);
		

    // now insert the relevant part of the solution into the mesh constraints
    IndexSet constrained_dofs;
    DoFTools::extract_boundary_dofs(free_surface_dof_handler, ComponentMask(dim, true),
                                    constrained_dofs, sim.parameters.free_surface_boundary_indicators);
    for ( unsigned int i = 0; i < constrained_dofs.n_elements();  ++i)
      {
        types::global_dof_index index = constrained_dofs.nth_index_in_set(i);
        if (mesh_displacement_constraints.can_store_line(index))
          if (mesh_displacement_constraints.is_constrained(index)==false)
            {
              mesh_displacement_constraints.add_line(index);
              mesh_displacement_constraints.set_inhomogeneity(index, boundary_velocity[index]);
            }
      }

    mesh_displacement_constraints.close();
  }

  
  
	struct xu
	{
		
		//xu(): x(0),y(0), jxw(0) { }
		double x = 0;
		double y = 0;
		double jxw = 0;
		
	};
	bool xuCompare(const xu& firstElem, const xu& secondElem) 
	{
		return firstElem.x < secondElem.x;
	}  
	
	double f(double x, double L=1)
	{
		
		//^-shaped
		//return -std::abs(L/2-x)*0.2+0.1;
		
		//quadratic
		return 0.4*x-0.4*x*x;
	}
	double F_def(double L=1, double n=1)
	{
		double pn = M_PI*n;
		
		//quadratic
		return -(2*pn*std::sin(pn)+4*std::cos(pn)-4)/(5*pn*pn*pn);
	}
	
	template <int dim>
  void FreeSurfaceHandler<dim>::diffuse_surface(LinearAlgebra::Vector &output)
	{
		
		
		
		
		// std::cout<<"entered diffuse_surface\n";
		const unsigned int ts = sim.timestep_number;
		
		// std::cout <<"Time = " <<sim.time <<std::endl;
		LinearAlgebra::SparseMatrix system_matrix;
		LinearAlgebra::Vector system_rhs, solution;
    system_rhs.reinit(mesh_locally_owned, sim.mpi_communicator);
    solution.reinit(mesh_locally_owned, sim.mpi_communicator);		
		
		
		//set up constraints
		
    ConstraintMatrix constraints(mesh_locally_relevant);
    DoFTools::make_hanging_node_constraints(free_surface_dof_handler, constraints);
	//make an empty oconstraint matrix
	ConstraintMatrix blank_constraints(mesh_locally_relevant);

	
    // typedef std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int> > periodic_boundary_pairs;
    // periodic_boundary_pairs pbp = sim.geometry_model->get_periodic_boundary_pairs();
    // for (periodic_boundary_pairs::iterator p = pbp.begin(); p != pbp.end(); ++p)
      // DoFTools::make_periodicity_constraints(free_surface_dof_handler,
                                             // (*p).first.first, (*p).first.second, (*p).second, mass_matrix_constraints);
    
		// //Zero out the displacement for the zero-velocity boundary indicators
    // VectorTools::interpolate_boundary_values (free_surface_dof_handler, 0,
                                                // ZeroFunction<dim>(dim), constraints);		
    // VectorTools::interpolate_boundary_values (free_surface_dof_handler, 1,
                                                // ZeroFunction<dim>(dim), constraints);																									
		
		
	// std::cout<<"completed periodic boundaries loop\n";
    constraints.close();
		
#ifdef ASPECT_USE_PETSC
    LinearAlgebra::DynamicSparsityPattern sparsity_pattern(mesh_locally_relevant);
#else
    TrilinosWrappers::SparsityPattern sparsity_pattern (mesh_locally_owned,
                                          mesh_locally_owned,
                                          mesh_locally_relevant,
                                          sim.mpi_communicator);
#endif
    DoFTools::make_sparsity_pattern (free_surface_dof_handler, sparsity_pattern, constraints, false,
                                     Utilities::MPI::this_mpi_process(sim.mpi_communicator));
#ifdef ASPECT_USE_PETSC
    SparsityTools::distribute_sparsity_pattern(sparsity_pattern,
                                               free_surface_dof_handler.n_locally_owned_dofs_per_processor(),
                                               sim.mpi_communicator, mesh_locally_relevant);

    sparsity_pattern.compress();
    system_matrix.reinit (mesh_locally_owned, mesh_locally_owned, sparsity_pattern, sim.mpi_communicator);
#else
    sparsity_pattern.compress();
    system_matrix.reinit (sparsity_pattern);
#endif
	

	// std::cout<<"completed sparsity patterns setup\n";
	
		
		
		
		
		
		QGauss<dim-1>  face_quadrature(free_surface_fe.degree);
		QGauss<dim>  quadrature(free_surface_fe.degree+1);
		// QTrapez<dim-1>  face_quadrature;
		// QTrapez<dim>  quadrature;
    UpdateFlags update_flags = UpdateFlags(update_values | update_gradients | update_JxW_values | update_quadrature_points | update_normal_vectors);
	UpdateFlags update_flags2 = UpdateFlags(update_values | update_gradients | update_JxW_values | update_quadrature_points );
    FEFaceValues<dim> fs_fe_face_values (*sim.mapping, free_surface_fe, face_quadrature, update_flags);
	FEValues<dim> fs_fe_values (*sim.mapping, free_surface_fe, quadrature, update_flags2);
	FEValues<dim> fe_values (*sim.mapping, sim.finite_element, quadrature, update_flags2);
	FEFaceValues<dim> fe_face_values (*sim.mapping, sim.finite_element, face_quadrature, update_flags);
		LinearAlgebra::Vector displacements = this->mesh_displacements;//free_surface_dof_handler.mesh_displacements;

		
		// std::cout<<"completed FEValues setup\n";
		
		double theta = 1;
		
		
		//Using deal.II tools
		// SparseMatrix<double> mass_matrix;
		// SparseMatrix<double> laplace_matrix;	
		// SparsityPattern      dealii_sparsity_pattern;
		// ConstraintMatrix     dealii_constraints;
		// SparseMatrix<double> dealii_system_matrix;
		
		// dealii_constraints.clear ();
		// DoFTools::make_hanging_node_constraints (free_surface_dof_handler,
												 // dealii_constraints);
		// dealii_constraints.close();

		// DynamicSparsityPattern dsp(free_surface_dof_handler.n_dofs());
		// DoFTools::make_sparsity_pattern(free_surface_dof_handler,
										// dsp,
										// dealii_constraints,
										// /*keep_constrained_dofs = */ true);
		// dealii_sparsity_pattern.copy_from(dsp);		
		
		
		// mass_matrix.reinit(sparsity_pattern);
		// laplace_matrix.reinit(sparsity_pattern);		
		
		// MatrixCreator::create_mass_matrix(free_surface_dof_handler,
								  // quadrature,
								  // mass_matrix);
		// MatrixCreator::create_laplace_matrix(free_surface_dof_handler,
											 // quadrature,
											 // laplace_matrix);
											 
		// dealii_system_matrix.copy_from(mass_matrix);
        // dealii_system_matrix.add(theta * sim.time_step * diffusivity, laplace_matrix);		

		
		//Output resulting matrix
		// if (sim.timestep_number == 0)	
		// {
			// std::ofstream file_dealii_system_matrix("dealii_system_matrix", std::ios::app);
			// for (int i =0; i<dealii_system_matrix.m();i++)
				// file_dealii_system_matrix << dealii_system_matrix[i]<<" ";
			// file_dealii_system_matrix << std::endl;
			// file_dealii_system_matrix.close();	
		// }
		
		// std::ofstream file_dealii_system_matrix("dealii_system_matrix", std::ios::app);
		// dealii_system_matrix.print	(file_dealii_system_matrix);			
		
		
		
		
		//Shortcuts
		const unsigned int dofs_per_cell = free_surface_fe.dofs_per_cell;
		const unsigned int dofs_per_cell_reg = sim.finite_element.dofs_per_cell;
		const unsigned int dofs_per_face = free_surface_fe.dofs_per_face;
		
		//const unsigned int dofs_per_cell = free_surface_fe.dofs_per_cell;
		// std::cout<<"dofs_per_cell="<<dofs_per_cell<<std::endl;
		// std::cout<<"dofs_per_cell_reg="<<dofs_per_cell_reg<<std::endl;
		// std::cout<<"dofs_per_face="<<dofs_per_face<<std::endl;
		
		const unsigned int n_q_points = fs_fe_values.n_quadrature_points;
		const unsigned int n_q_points_reg = fe_values.n_quadrature_points;
		const unsigned int n_face_q_points = fs_fe_face_values.n_quadrature_points;
		// std::cout<<"n_q_points = "<<n_q_points<<std::endl;
		// std::cout<<"n_q_points_reg = "<<n_q_points_reg<<std::endl;
		// std::cout<<"n_face_q_points = "<<n_face_q_points<<std::endl;
		//Create local matrices
		FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
		// FullMatrix<double>   cell_mass (dofs_per_cell, dofs_per_cell);
		// FullMatrix<double>   cell_laplace (dofs_per_cell, dofs_per_cell);
		Vector<double>       cell_rhs (dofs_per_cell);
		//Vector to store global numbers of local DOFs
		std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
		
		
		//Create cell iterator
		typename DoFHandler<dim>::active_cell_iterator
		cell = sim.dof_handler.begin_active(),
		endc = sim.dof_handler.end();
		typename DoFHandler<dim>::active_cell_iterator
		fscell = free_surface_dof_handler.begin_active(),
		endfscell = free_surface_dof_handler.end();
		
		FEValuesExtractors::Scalar  vertical_displacement(dim-1);
		std::vector<double> local_displacements (n_face_q_points);
		
		
		
		//std::cout<<"Reached loop over cells\n";
		
		
		typedef std::vector<xu> displacements_vector;
		typedef std::vector<displacements_vector> displacements_matrix;
		static displacements_matrix u;
		u.push_back(displacements_vector());	


		std::vector<double> disp;
		double max_rhs = -999;
		double norm = 0;
		
                int cell_counter = 0;
                // std::ofstream assembly("assembly", std::ios::app);
		
                // assembly <<std::endl;
                // assembly << "=========== Timestep "<<sim.timestep_number<<" ===========";
                // assembly <<std::endl;
                
		//Loop over cells
		for (; cell!=endc; ++cell, ++fscell)
		{
			//We're only interested in boundary cells
			if (cell->at_boundary() /* && cell->boundary_id()==dim-1 */ && cell->is_locally_owned())
				for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
					//...and specifially, in faces lying on the boundary
					{
					//if (cell->face(face_no)->at_boundary() && cell->face(face_no)->boundary_id()==dim*2-1)
					const types::boundary_id boundary_indicator = cell->face(face_no)->boundary_id();
					if (sim.parameters.free_surface_boundary_indicators.find(boundary_indicator)
						== sim.parameters.free_surface_boundary_indicators.end())
						continue;	
                                                
                                        
                                                // assembly <<std::endl;
                                                // assembly << "=========== Cell no. "<<cell_counter++<<" ===========";
                                                // assembly <<std::endl;
				
						//recompute values, gradients of the shape functions and determinants of the Jacobian matrices for the current cell
						fs_fe_face_values.reinit (fscell, face_no);
						//Reset the local cell's contributions to global matrix
						cell_matrix = 0;
						// cell_mass = 0;
						// cell_laplace = 0;
						cell_rhs = 0;
						//Fill global numbers of local DOFs for current cell
						fscell->get_dof_indices (local_dof_indices);
						
						fs_fe_face_values[vertical_displacement].get_function_values (displacements,
																						local_displacements);
						
						
						
						
						//Integrate over the fscell by looping over quadrature points
						for (unsigned int q_index=0; q_index<n_face_q_points; ++q_index)
						{
                                                        const Point<dim> p = fs_fe_face_values.quadrature_point(q_index);
							double displacement = local_displacements[q_index];
							disp.push_back(displacement);
							
							//if (sim.timestep_number==0) std::cout<<fs_fe_face_values.JxW (q_index)<<" ";
							// if (sim.timestep_number==0) std::cout<<fs_fe_face_values.JxW (q_index)<<" ";
							
							
							
							
							for (unsigned int i=0; i<dofs_per_cell; ++i)
							{
								if (free_surface_fe.system_to_component_index(i).first != dim-1)
									continue;
								//if (sim.timestep_number==0) std::cout<<"dof.no.="<<i<<" ";
								//if (sim.timestep_number==0) std::cout<<"c.i="<<free_surface_fe.system_to_component_index(i).first<<" ";
								
								//Set local matrix elements					
								Tensor<2, dim, double> rot;
								outer_product(rot, fs_fe_face_values.normal_vector(q_index), fs_fe_face_values.normal_vector(q_index));
								Tensor<2, dim, double> I;
								for (int k = 0; k<dim; k++)
									I[k][k]=1;
								rot =  I - rot;
								
								//Set right hand side -- do we need rot vector? 
								double rhs_contribution = (fs_fe_face_values.shape_value (i, q_index) 
												* displacement 
												- sim.time_step*diffusivity*(1-theta)*fs_fe_face_values.shape_value (i, q_index)*displacement)
												* fs_fe_face_values.JxW (q_index);
								cell_rhs(i) += rhs_contribution;
								
								// assembly << "q_index="<<q_index<<" q_x="<<p[0]<<" i="<<i<<" shape_value="<<fs_fe_face_values.shape_value(i, q_index)<<" shape_grad="<<fs_fe_face_values.shape_grad(i, q_index)<<" rhs_contribution="<<rhs_contribution<<" ";
								// assembly <<std::endl;
                                                                
								// if (sim.timestep_number==0) std::cout<<fs_fe_face_values.shape_value (i, q_index)<<" ";

								

								// if (fs_fe_face_values.normal_vector(q_index)[0] || fs_fe_face_values.normal_vector(q_index)[1])
									// std::cout<<fs_fe_face_values.normal_vector(q_index)[0]<<" "<<fs_fe_face_values.normal_vector(q_index)[1]<<std::endl;
								//normal vectors are fine

								//if (rot.norm())
									//std::cout<<rot[0][0]<<" "<<rot[0][1]<<rot[1][0]<<" "<<rot[1][1]<<std::endl;

								
								//rot =  dealii::identity_tensor<dim> () - rot;
								
								
                                                                
                                                                
								for (unsigned int j=0; j<dofs_per_cell; ++j)
								{
									if (free_surface_fe.system_to_component_index(j).first != dim-1)
										continue;
									/*fs_fe_face_values.shape_grad (i|j, q_index) - gradients of shape functions*/
									
									double matrix_contribution = (sim.time_step*diffusivity * theta * 
									(rot*fs_fe_face_values.shape_grad (i, q_index)  * rot*fs_fe_face_values.shape_grad (j, q_index)) 
										+  (fs_fe_face_values.shape_value (i, q_index) * fs_fe_face_values.shape_value (j, q_index)))
										* fs_fe_face_values.JxW (q_index); 
									
									cell_matrix(i,j) += matrix_contribution;
                                                                        
									// assembly << "q_index="<<q_index<<" q_x="<<p[0]<<" i="<<i<<" j="<<j<<" matrix_contribution="<<matrix_contribution<<" ";
									// assembly <<std::endl;
																		
									// cell_laplace(i,j) += (fe_values.shape_grad (i, q_index)  * fe_values.shape_grad (j, q_index)) 
															// * fe_values.JxW (q_index); 
										
									// cell_mass(i,j) += (fe_values.shape_value (i, q_index) * fe_values.shape_value (j, q_index))
															// * fe_values.JxW (q_index); 																		
                                                                        
										/*
									cell_rhs(i,j) += (displacement * sim.time_step * diffusivity * (1-theta) * (rot*fs_fe_face_values.shape_grad (i, q_index)  * rot*fs_fe_face_values.shape_grad (j, q_index)) 
										+  (fs_fe_face_values.shape_value (i, q_index) * fs_fe_face_values.shape_value (j, q_index)))
										* fs_fe_face_values.JxW (q_index); */

										// need SURFACE gradient
									// if (sim.time_step*diffusivity * (fs_fe_face_values.shape_grad (i, q_index)  * fs_fe_face_values.shape_grad (j, q_index))> 1e-13)
										// std::cout<<sim.time_step*diffusivity * (fs_fe_face_values.shape_grad (i, q_index)  * fs_fe_face_values.shape_grad (j, q_index)) <<" ";
									//around 1e-15
									//if (fs_fe_face_values.shape_value (i, q_index) * fs_fe_face_values.shape_value (j, q_index) > 1e-13)
										//std::cout<<fs_fe_face_values.shape_value (i, q_index) * fs_fe_face_values.shape_value (j, q_index)<<" ";
									//around 0.1~0.6
									// if (fs_fe_face_values.shape_grad (i, q_index).norm()  > 1e-13)
										// std::cout<<" q_index="<<q_index<<" i="<<i<<" grad="<<fs_fe_face_values.shape_grad (i, q_index);
									//around 1
									
									//no diffusion
									// cell_matrix(i,j) += (fs_fe_face_values.shape_value(i, q_index) * fs_fe_face_values.shape_value(j, q_index))
										// * fs_fe_face_values.JxW (q_index);            // need SURFACE gradient				
									
								}	
							}
                                                        
                                                        
                                                        	
							
							
							xu temp;
							temp.x=p[0];
							temp.y=displacement;//p[dim-1];
							temp.jxw = fs_fe_face_values.JxW (q_index);
							u[ts].push_back(temp);							
							
						}

						
																																	
																																
						// if (sim.timestep_number == 0)	
						// {							
							// std::ofstream myfile_matrix("cell_matrix", std::ios::app);
							
							// for (unsigned int i=0; i<cell_matrix.m(); i++) 
							// {	
								// for (unsigned int j=0; j<cell_matrix.n(); j++)
									// myfile_matrix << cell_matrix[i][j]<<" ";
								// myfile_matrix <<std::endl;
							// }
							// myfile_matrix.close();		

							// std::ofstream myfile_rhs("rhs", std::ios::app);
							// for (unsigned int i=0; i<cell_rhs.size(); ++i) 
								// myfile_rhs << cell_rhs[i]<<" ";
							// myfile_rhs.close();							
						// }
						
						constraints.distribute_local_to_global (cell_matrix, cell_rhs,
                                                                  local_dof_indices, system_matrix, system_rhs, false);
						
						//std::cout<<"system_matrix.m()="<<system_matrix.m()<<"\n";
						//std::cout<<"system_matrix.n()="<<system_matrix.n()<<"\n";
						// for (unsigned int i=0; i<dofs_per_cell; ++i)
						// {
							// //transfer the local rhs to the global matrix
							
							// system_rhs(local_dof_indices[i]) += cell_rhs(i);
							// for (unsigned int j=0; j<dofs_per_cell; ++j)
							// {
								// //transfer the local matrix elements to the global matrix
								// std::cout<<"i="<<i<<"\n";
								// std::cout<<"j="<<j<<"\n";
								// std::cout<<"local_dof_indices[i]="<<local_dof_indices[i]<<"\n";
								// std::cout<<"local_dof_indices[j]="<<local_dof_indices[j]<<"\n";
								// system_matrix.add (local_dof_indices[i],
																	 // local_dof_indices[j],
																	 // cell_matrix(i,j));						
							// }
						// }
					}
		}
                
                //assembly.close();
                
		// std::cout<<"ts = "<<ts<<"sim.timestep_number = "<<sim.timestep_number<<std::endl;
		// std::cout<<"len (u[ts]) = "<<u[ts].size()<<std::endl;
		std::sort(u[ts].begin(), u[ts].end(), xuCompare);
		// for (unsigned int j=0; j<u[ts].size(); ++j) 
			// std::cout<<u[ts][j].x<<" ";
		
		
		
		std::ofstream myfile("u", std::ios::app);
		// std::cout<<"Max disp = "<<*max_element(std::begin(disp), std::end(disp))<<std::endl;
		
		
		
		// for (int i = 0; i<system_rhs.size();i++)
			// if (max_rhs<system_rhs(i)) max_rhs=system_rhs(i);
		// std::cout<<"Max rhs = "<<max_rhs<<std::endl;
		// //std::cout<<"Max rhs = "<<*max_element(std::begin(system_rhs), std::end(system_rhs))<<std::endl;
		
		for (unsigned int j=0; j<u[ts].size(); ++j) 
			{
				myfile << u[ts][j].x<<" ";//std::ostream_iterator<int>(std::cout, " "));
				myfile << u[ts][j].y<<std::endl;
			}
		myfile << std::endl;
		myfile.close();		
		
		
		
		//std::cout<<"Exit loop over cells\n";
		//Generate a list of pairs of global degree of freedom numbers (on the boundary) 
		//and their boundary values (which are zero here for all entries)
		
		/*
		std::map<types::global_dof_index, double> boundary_values;
		switch (dim)
		{
			case 2:
			{
				
				VectorTools::interpolate_boundary_values (*sim.mapping,
																									free_surface_dof_handler,
																									3,
																									ZeroFunction<dim>(2),
																									boundary_values);				
			}
			case 3:
			{
				VectorTools::interpolate_boundary_values (*sim.mapping,
																									free_surface_dof_handler,//sim.dof_handler,
																									5,
																									ZeroFunction<dim>(),
																									boundary_values);				
			}
		}*/
		// Apply b.c. to the system of equations
		// MatrixTools::apply_boundary_values (constraints,
																				// system_matrix,
																				// solution,
																				// system_rhs);		
		
		system_rhs.compress (VectorOperation::add);
		system_matrix.compress(VectorOperation::add);	
		
		//Set parameters for the solver: (max_iterations, tolerance); 
		//std::cout<<"rhs size="<<system_rhs.size()<<"\n";
		SolverControl solver_control (1000000, 1e-7);
		//SolverControl solver_control (5*system_rhs.size(), sim.parameters.linear_stokes_solver_tolerance*system_rhs.l2_norm());
		//Create solver itself
		SolverCG<LinearAlgebra::Vector> solver (solver_control);
		//Solve
		
		solver.solve (system_matrix, solution, system_rhs,
									PreconditionIdentity());
									
		//std::cout<<"Diffusion equation solved in "<<solver_control.last_step() << " iterations."<< std::endl;							
		
		// Do I need this?? ************
		//constraints.distribute (solution);
		
		LinearAlgebra::Vector output_temp(mesh_locally_owned, mesh_locally_relevant, sim.mpi_communicator);
		// std::cout<<"output_temp size:" <<output_temp.size()<<"\n";
		// std::cout<<"solution size:" <<solution.size()<<"\n";
		// std::cout<<"output size:" <<output.size()<<"\n";
		// std::cout<<"displacements size:" <<displacements.size()<<"\n";
		output_temp = solution;
		//std::cout<<"solution size ="<<solution.size()<<std::endl;
		//std::cout<<"rhs size ="<<system_rhs.size()<<std::endl;
		
		// if (sim.timestep_number==0)
		// {
			// std::ofstream file_system_matrix("system_matrix", std::ios::app);
			// for (int i =0; i<system_matrix.size();i++)
				// file_system_matrix << system_matrix[i]<<" ";
			// file_system_matrix << std::endl;
			// file_system_matrix.close();		
		// }
		
		
		// std::ofstream file_sol("sol", std::ios::app);
		// for (int i =0; i<solution.size();i++)
			// file_sol << solution[i]<<" ";
		// file_sol << std::endl;
		// file_sol.close();
			
		// std::ofstream file_system_rhs("system_rhs", std::ios::app);
		// for (int i =0; i<system_rhs.size();i++)
			// file_system_rhs << system_rhs[i]<<" ";
		// file_system_rhs << std::endl;
		// file_system_rhs.close();			
		
		output_temp -= displacements; //resulting elevation - initial elevation
		//UNCOMMENT//output_temp.sadd(1,displacements);
		//output_temp.sadd(1/sim.time_step,);
		
		//displacement[dim-1] = p[dim-1]*0.1*std::sin(p[0]*M_PI);
		
		//=====
		
		//Analytical solution
		//if (MPI_COMM_RANK==0)
		//{
			//Checking cooling from constant initial temperature
			// double u0 = 0.1;
			// for (int i=0; i<u[0].size(); i++)
				// u[0][i].y=u0;
			
			
			double time_ua = std::max(sim.time-sim.time_step/*sim.old_time_step*/,0.0);
			
			const double L = 1; //x extent of the domain
			//const double T = L*L / diffusivity; //charachteristic temperature
			//const double U = 0.1; //initial displacement amplitude
			
			static displacements_matrix ua;
			std::vector<double> u_diff(u[ts].size(), 0);
			std::vector<double> ua_diff(u[ts].size(), 0);
			ua.push_back(displacements_vector());
			
			// std::cout<<"u[0]size="<<u[0].size()<<" ";
			
			double eps = 1e-10;
			
			//double max_x = 1;//sim.geometry_model->get_extents()[0];
			
			for (int i = 0; i<u[0].size(); ++i)
			{
				xu temp;
				temp.x=u[ts][i].x;
				
				int counter = 0;
				double series = 0;
				for (int n=1; n<=1000; n++)
				{
					
					//cooling of a rod from constant initial temperature
					//double term = exp(-(2*n-1)*(2*n-1)*M_PI*M_PI*time_ua) * std::sin((2*n-1)*M_PI*temp.x)/(2*n-1);
					
					// if (n==1)
					// {
					// //explicit Bn
					// double Bn = 4 * u[0][i].y / M_PI;
					// std::cout << "Bn expl = "<<Bn<<" ";
					
					// //numerical Bn calculation
					// double Bnn = 0; double y = 0;
					// double num_segments = 1000;
					// double step = L/num_segments;
					// double np = n*M_PI;
					// for (int k = 0; k<num_segments; k++)
					// {
						// //std::cout << "k = "<<k <<" ";
						// y = std::sin(np*(k+0.5)*step) * 0.1;
						// Bnn = Bnn + 2* y * step;
					// }
					// std::cout << "Bnn = "<<Bnn<<" ";
					// }
					//Diffusing a hat shaped initial topography
					double pn = n*M_PI;
					
					//explicit Bn calculation
					//double Bn = (-0.2/np)*(L/2 - (2*std::sin(np*L/2)/np) - L/2*std::cos (np*L) + std::sin(np*L)/np);
					//online calculated formula
					//double Bn1 = -113*std::sin(pn)/(1775*M_PI*n*n) + (std::cos(3927*n/2500)/(10*pn) - (35969*std::sin(3927*n/2500)/(1000*n)
						//- std::cos(3927*n/2500)/7569)/(1775*n));
					//if (n==1) std::cout << "Bn1 = "<<Bn1 <<" ";
					//double Bn=4/(5*M_PI*M_PI);
					//numerical Bn calculation
					// double Bn = 0; double y = 0;
					// double num_segments = 1000;
					// double step = L/num_segments;
					// for (int k = 0; k<num_segments; k++)
					// {
						// //std::cout << "k = "<<k <<" ";
						// y = std::sin(np*k*step) * ( -0.2 * std::abs(L/2 - k*step) + 0.1);
						// Bn = Bn + y * step;
					// }
					
					double Bn = 2*F_def(1,n);
					double term = Bn * std::sin (pn*temp.x) * std::exp(-pn*pn*diffusivity*time_ua);
					series += term;
					//if (n==1)
					//{
						//std::cout << "Bn = "<<Bn<<" ";
						// std::cout << "exp = "<<std::exp(-np*np*diffusivity*time_ua)<<" ";
						//std::cout << "series = "<<series <<" ";
						// std::cout << "term = "<<term <<" ";					
					//}				
					
					if (std::abs(term) < eps)
						counter++;
					else 
						counter=0;
					if (counter>3)
					{
						// std::cout << "n reached= "<<n<<" ";
						break;
					}
				}
				//if (i==0)
					//std::cout << std::endl<<"n (series)= " <<n<<" series = "<< series <<std::endl;
				
				//if (ts==0)
					//temp.y = u[ts][i].y;
				//else
					//temp.y=4 * u[0][i].y * series / M_PI ;
					temp.y = series;
				ua[ts].push_back(temp);				
				
				
			}
		//}
		
		
		std::ofstream myfile2("ua", std::ios::app);
		std::ofstream myfile3("u_diff", std::ios::app);
		std::ofstream myfile4("ua_diff", std::ios::app);
		std::ofstream myfile5("u-ua", std::ios::app);
		
		
		
		
		for (unsigned int j=0; j<ua[ts].size(); ++j) 
		{
			myfile2 << ua[ts][j].x<<" ";
			myfile2 << ua[ts][j].y<<std::endl;
			myfile5 << ua[ts][j].x<<" ";
			myfile5 << u[ts][j].y-ua[ts][j].y<<std::endl;
			
			if (ts)
			{
				if (j==15)
				{
					myfile3 << ts<<" ";
					myfile3 << u[ts][j].y-u[ts-1][j].y<<std::endl;
					myfile4 << ts<<" ";
					myfile4 << ua[ts][j].y-ua[ts-1][j].y<<std::endl;
					
				}
			}
			
		}
			
			
		// calculating norm of error against resolution
		
		std::ofstream myfile6("u-ua_norm", std::ios::app);
		if (sim.timestep_number==3)
		{
			for (unsigned int j=0; j<ua[ts].size(); ++j) 
				norm += pow(u[ts][j].y-ua[ts][j].y, 2)*u[ts][j].jxw;
			norm = pow(norm, 0.5);
			myfile6 << "Refinement "<<sim.parameters.initial_global_refinement<<"+"<<sim.parameters.initial_adaptive_refinement<<" ";
			myfile6 << norm<<std::endl;
		}
		
			
		myfile2 << std::endl;
		myfile5 << std::endl;
		// myfile3 << std::endl;
		// myfile4 << std::endl;
		myfile2.close();			
		myfile3.close();			
		myfile4.close();			
		myfile5.close();
		myfile6.close();
		//Compare obtained solution with analytical solution
		
		
		

		
		//std::ofstream out ("sparsity_pattern2.svg");
		//sparsity_pattern.print_svg (out);		
		
		
		
		
		//std::cout<<"Reached last for loop\n";
		//for (unsigned int i = 0; i < solution.size(); i++)
		//{
			//output_temp[i] = solution[i];// - displacements[i];
		if (sim.time_step)
			output_temp /= sim.time_step;
		//}
		//std::cout<<"Exit last for loop\n";
		output=output_temp;
		
		
		
		
		
		

/*  cg.solve (mesh_matrix, velocity_solution, rhs, preconditioner_stiffness);
    sim.pcout << "   Solving mesh velocity system... " << solver_control.last_step() <<" iterations."<< std::endl;

    mass_matrix_constraints.distribute (solution);

    //Update the free surface mesh velocity vector
    fs_mesh_velocity = velocity_solution;

    //Update the mesh displacement vector
    LinearAlgebra::Vector distributed_mesh_displacements(mesh_locally_owned, sim.mpi_communicator);
    distributed_mesh_displacements = mesh_displacements;
    distributed_mesh_displacements.add(sim.time_step, velocity_solution);
		
    mesh_displacements = distributed_mesh_displacements; */		
		
		
		//std::cout<<"Reached end of important part\n";
		
		
		
		//Visual output
		// //std::cout<<"sim.time_step = "<<sim.time_step<<"\n";
		// Vector<double> diff;
		// diff = solution;
		// diff -= displacements;
		
		// {
			// DataOut<dim> data_out;
			// data_out.attach_dof_handler (free_surface_dof_handler);
			// data_out.add_data_vector (solution, "solution");
			// data_out.add_data_vector (diff, "diff");
			// data_out.add_data_vector (displacements, "displacements");
			// data_out.build_patches ();
			
			// std::ofstream out ("solutions.gpl");
			// data_out.write_gnuplot (out);
		// }
	}
	

  template <int dim>
  void FreeSurfaceHandler<dim>::project_velocity_onto_boundary(LinearAlgebra::Vector &output)
  {
    // TODO: should we use the extrapolated solution?

    // stuff for iterating over the mesh
    QGauss<dim-1> face_quadrature(free_surface_fe.degree+1);
    UpdateFlags update_flags = UpdateFlags(update_values | update_quadrature_points
                                           | update_normal_vectors | update_JxW_values);
    FEFaceValues<dim> fs_fe_face_values (*sim.mapping, free_surface_fe, face_quadrature, update_flags);
    FEFaceValues<dim> fe_face_values (*sim.mapping, sim.finite_element, face_quadrature, update_flags);
    const unsigned int n_face_q_points = fe_face_values.n_quadrature_points,
                       dofs_per_cell = fs_fe_face_values.dofs_per_cell;

    // stuff for assembling system
    std::vector<types::global_dof_index> cell_dof_indices (dofs_per_cell);
    Vector<double> cell_vector (dofs_per_cell);
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    // stuff for getting the velocity values
    std::vector<Tensor<1,dim> > velocity_values(n_face_q_points);

    // set up constraints
    ConstraintMatrix mass_matrix_constraints(mesh_locally_relevant);
    DoFTools::make_hanging_node_constraints(free_surface_dof_handler, mass_matrix_constraints);

    typedef std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int> > periodic_boundary_pairs;
    periodic_boundary_pairs pbp = sim.geometry_model->get_periodic_boundary_pairs();
    for (periodic_boundary_pairs::iterator p = pbp.begin(); p != pbp.end(); ++p)
      DoFTools::make_periodicity_constraints(free_surface_dof_handler,
                                             (*p).first.first, (*p).first.second, (*p).second, mass_matrix_constraints);

    mass_matrix_constraints.close();

    // set up the matrix
    LinearAlgebra::SparseMatrix mass_matrix;
#ifdef ASPECT_USE_PETSC
    LinearAlgebra::DynamicSparsityPattern sp(mesh_locally_relevant);

#else
    TrilinosWrappers::SparsityPattern sp (mesh_locally_owned,
                                          mesh_locally_owned,
                                          mesh_locally_relevant,
                                          sim.mpi_communicator);
#endif
    DoFTools::make_sparsity_pattern (free_surface_dof_handler, sp, mass_matrix_constraints, false,
                                     Utilities::MPI::this_mpi_process(sim.mpi_communicator));
#ifdef ASPECT_USE_PETSC
    SparsityTools::distribute_sparsity_pattern(sp,
                                               free_surface_dof_handler.n_locally_owned_dofs_per_processor(),
                                               sim.mpi_communicator, mesh_locally_relevant);

    sp.compress();
    mass_matrix.reinit (mesh_locally_owned, mesh_locally_owned, sp, sim.mpi_communicator);
#else
    sp.compress();
    mass_matrix.reinit (sp);
#endif

    FEValuesExtractors::Vector extract_vel(0);

    // make distributed vectors.
    LinearAlgebra::Vector rhs, dist_solution;
    rhs.reinit(mesh_locally_owned, sim.mpi_communicator);
    dist_solution.reinit(mesh_locally_owned, sim.mpi_communicator);

    typename DoFHandler<dim>::active_cell_iterator
    cell = sim.dof_handler.begin_active(), endc= sim.dof_handler.end();
    typename DoFHandler<dim>::active_cell_iterator
    fscell = free_surface_dof_handler.begin_active();

    for (; cell!=endc; ++cell, ++fscell)
      if (cell->at_boundary() && cell->is_locally_owned())
        for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
          if (cell->face(face_no)->at_boundary())
            {
              const types::boundary_id boundary_indicator
                = cell->face(face_no)->boundary_id();
              if (sim.parameters.free_surface_boundary_indicators.find(boundary_indicator)
                  == sim.parameters.free_surface_boundary_indicators.end())
                continue;

              fscell->get_dof_indices (cell_dof_indices);
              fs_fe_face_values.reinit (fscell, face_no);
              fe_face_values.reinit (cell, face_no);
              fe_face_values[sim.introspection.extractors.velocities].get_function_values(sim.solution, velocity_values);

              cell_vector = 0;
              cell_matrix = 0;
              for (unsigned int point=0; point<n_face_q_points; ++point)
                {
                  // Select the direction onto which to project the velocity solution
                  Tensor<1,dim> direction;
                  if ( advection_direction == SurfaceAdvection::normal ) // project onto normal vector
                    direction = fs_fe_face_values.normal_vector(point);
                  else if ( advection_direction == SurfaceAdvection::vertical ) // project onto local gravity
                    direction = sim.gravity_model->gravity_vector( fs_fe_face_values.quadrature_point(point) );
                  else AssertThrow(false, ExcInternalError());
                  direction *= ( direction.norm() > 0.0 ? 1./direction.norm() : 0.0 );

                  for (unsigned int i=0; i<dofs_per_cell; ++i)
                    {
                      for (unsigned int j=0; j<dofs_per_cell; ++j)
                        {
                          cell_matrix(i,j) += (fs_fe_face_values[extract_vel].value(j,point) *
                                               fs_fe_face_values[extract_vel].value(i,point) ) *
                                              fs_fe_face_values.JxW(point);
                        }

                      cell_vector(i) += (fs_fe_face_values[extract_vel].value(i,point) * direction)
                                        * (velocity_values[point] * direction)
                                        * fs_fe_face_values.JxW(point);
                    }
                }

              mass_matrix_constraints.distribute_local_to_global (cell_matrix, cell_vector,
                                                                  cell_dof_indices, mass_matrix, rhs, false);
            }

    rhs.compress (VectorOperation::add);
    mass_matrix.compress(VectorOperation::add);

    // Jacobi seems to be fine here.  Other preconditioners (ILU, IC) run into troubles
    // because the matrix is mostly empty, since we don't touch internal vertices.
    LinearAlgebra::PreconditionJacobi preconditioner_mass;
    preconditioner_mass.initialize(mass_matrix);

    SolverControl solver_control(5*rhs.size(), sim.parameters.linear_stokes_solver_tolerance*rhs.l2_norm());
    SolverCG<LinearAlgebra::Vector> cg(solver_control);
    cg.solve (mass_matrix, dist_solution, rhs, preconditioner_mass);

    mass_matrix_constraints.distribute (dist_solution);
    output = dist_solution;
	
	// std::ofstream file_sol("stokes_sol", std::ios::app);
	// for (int i =0; i<output.size();i++)
		// file_sol << output[i]<<" ";
	// file_sol << std::endl;
	// file_sol.close();	
	
	
	
  }


  template <int dim>
  void FreeSurfaceHandler<dim>::compute_mesh_displacements()
  {
    QGauss<dim> quadrature(free_surface_fe.degree + 1);
    UpdateFlags update_flags = UpdateFlags(update_values | update_JxW_values | update_gradients);
    FEValues<dim> fe_values (*sim.mapping, free_surface_fe, quadrature, update_flags);

    const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                       dofs_per_face = sim.finite_element.dofs_per_face,
                       n_q_points    = fe_values.n_quadrature_points;

    std::vector<types::global_dof_index> cell_dof_indices (dofs_per_cell);
    std::vector<unsigned int> face_dof_indices (dofs_per_face);
    Vector<double> cell_vector (dofs_per_cell);
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    // We are just solving a Laplacian in each spatial direction, so
    // the degrees of freedom for different dimensions do not couple.
    Table<2,DoFTools::Coupling> coupling (dim, dim);
    coupling.fill(DoFTools::none);

    for (unsigned int c=0; c<dim; ++c)
      coupling[c][c] = DoFTools::always;

    LinearAlgebra::SparseMatrix mesh_matrix;
#ifdef ASPECT_USE_PETSC
    LinearAlgebra::DynamicSparsityPattern sp(mesh_locally_relevant);
#else
    TrilinosWrappers::SparsityPattern sp (mesh_locally_owned,
                                          mesh_locally_owned,
                                          mesh_locally_relevant,
                                          sim.mpi_communicator);
#endif
    DoFTools::make_sparsity_pattern (free_surface_dof_handler,
                                     coupling, sp,
                                     mesh_displacement_constraints, false,
                                     Utilities::MPI::
                                     this_mpi_process(sim.mpi_communicator));
#ifdef ASPECT_USE_PETSC
    SparsityTools::distribute_sparsity_pattern(sp,
                                               free_surface_dof_handler.n_locally_owned_dofs_per_processor(),
                                               sim.mpi_communicator, mesh_locally_relevant);
    sp.compress();
    mesh_matrix.reinit (mesh_locally_owned, mesh_locally_owned, sp, sim.mpi_communicator);
#else
    sp.compress();
    mesh_matrix.reinit (sp);
#endif

    // carry out the solution
    FEValuesExtractors::Vector extract_vel(0);

    LinearAlgebra::Vector rhs, velocity_solution;
    rhs.reinit(mesh_locally_owned, sim.mpi_communicator);
    velocity_solution.reinit(mesh_locally_owned, sim.mpi_communicator);

    typename DoFHandler<dim>::active_cell_iterator cell = free_surface_dof_handler.begin_active(),
                                                   endc= free_surface_dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          cell->get_dof_indices (cell_dof_indices);
          fe_values.reinit (cell);

          cell_vector = 0;
          cell_matrix = 0;
          for (unsigned int point=0; point<n_q_points; ++point)
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              {
                for (unsigned int j=0; j<dofs_per_cell; ++j)
                  cell_matrix(i,j) += scalar_product( fe_values[extract_vel].gradient(j,point),
                                                      fe_values[extract_vel].gradient(i,point) ) *
                                      fe_values.JxW(point);
              }

          mesh_displacement_constraints.distribute_local_to_global (cell_matrix, cell_vector,
                                                                    cell_dof_indices, mesh_matrix, rhs, false);
        }

    rhs.compress (VectorOperation::add);
    mesh_matrix.compress (VectorOperation::add);

    // Make the AMG preconditioner
    std::vector<std::vector<bool> > constant_modes;
    DoFTools::extract_constant_modes (free_surface_dof_handler,
                                      ComponentMask(dim, true),
                                      constant_modes);
    // TODO: think about keeping object between time steps
    LinearAlgebra::PreconditionAMG preconditioner_stiffness;
    LinearAlgebra::PreconditionAMG::AdditionalData Amg_data;
#ifdef ASPECT_USE_PETSC
    Amg_data.symmetric_operator = false;
#else
    Amg_data.constant_modes = constant_modes;
    Amg_data.elliptic = true;
    Amg_data.higher_order_elements = false;
    Amg_data.smoother_sweeps = 2;
    Amg_data.aggregation_threshold = 0.02;
#endif
    preconditioner_stiffness.initialize(mesh_matrix);

    SolverControl solver_control(5*rhs.size(), sim.parameters.linear_stokes_solver_tolerance*rhs.l2_norm());
    SolverCG<LinearAlgebra::Vector> cg(solver_control);

    cg.solve (mesh_matrix, velocity_solution, rhs, preconditioner_stiffness);
    sim.pcout << "   Solving mesh velocity system... " << solver_control.last_step() <<" iterations."<< std::endl;

    mesh_displacement_constraints.distribute (velocity_solution);

    // Update the free surface mesh velocity vector
    fs_mesh_velocity = velocity_solution;

    // Update the mesh displacement vector
    LinearAlgebra::Vector distributed_mesh_displacements(mesh_locally_owned, sim.mpi_communicator);
    distributed_mesh_displacements = mesh_displacements;
    distributed_mesh_displacements.add(sim.time_step, velocity_solution);
		
    mesh_displacements = distributed_mesh_displacements;

  }


  template <int dim>
  void FreeSurfaceHandler<dim>::interpolate_mesh_velocity()
  {
    // Interpolate the mesh vertex velocity onto the Stokes velocity system for use in ALE corrections
    LinearAlgebra::BlockVector distributed_mesh_velocity;
    distributed_mesh_velocity.reinit(sim.introspection.index_sets.system_partitioning, sim.mpi_communicator);

    const std::vector<Point<dim> > support_points
      = sim.finite_element.base_element(sim.introspection.component_indices.velocities[0]).get_unit_support_points();

    Quadrature<dim> quad(support_points);
    UpdateFlags update_flags = UpdateFlags(update_values | update_JxW_values);
    FEValues<dim> fs_fe_values (*sim.mapping, free_surface_fe, quad, update_flags);
    FEValues<dim> fe_values (*sim.mapping, sim.finite_element, quad, update_flags);
    const unsigned int n_q_points = fe_values.n_quadrature_points,
                       dofs_per_cell = fe_values.dofs_per_cell;

    std::vector<types::global_dof_index> cell_dof_indices (dofs_per_cell);
    FEValuesExtractors::Vector extract_vel(0);
    std::vector<Tensor<1,dim> > velocity_values(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator
    cell = sim.dof_handler.begin_active(), endc= sim.dof_handler.end();
    typename DoFHandler<dim>::active_cell_iterator
    fscell = free_surface_dof_handler.begin_active();

    for (; cell!=endc; ++cell, ++fscell)
      if (cell->is_locally_owned())
        {
          cell->get_dof_indices (cell_dof_indices);

          fe_values.reinit (cell);
          fs_fe_values.reinit (fscell);
          fs_fe_values[extract_vel].get_function_values(fs_mesh_velocity, velocity_values);
          for (unsigned int j=0; j<n_q_points; ++j)
            for (unsigned int dir=0; dir<dim; ++dir)
              {
                unsigned int support_point_index
                  = sim.finite_element.component_to_system_index(/*velocity component=*/ sim.introspection.component_indices.velocities[dir],
                                                                                         /*dof index within component=*/ j);
                distributed_mesh_velocity[cell_dof_indices[support_point_index]] = velocity_values[j][dir];
              }
        }

    distributed_mesh_velocity.compress(VectorOperation::insert);
    mesh_velocity = distributed_mesh_velocity;
  }



	
  template <int dim>
  void FreeSurfaceHandler<dim>::setup_dofs() //called from core every mesh repartitioning
  {
    if (!sim.parameters.free_surface_enabled)
      return;

    // these live in the same FE as the velocity variable:
    mesh_velocity.reinit(sim.introspection.index_sets.system_partitioning,
                         sim.introspection.index_sets.system_relevant_partitioning,
                         sim.mpi_communicator);


    free_surface_dof_handler.distribute_dofs(free_surface_fe);

    sim.pcout << "Number of free surface degrees of freedom: "
              << free_surface_dof_handler.n_dofs()
              << std::endl;

    // Renumber the DoFs hierarchical so that we get the
    // same numbering if we resume the computation. This
    // is because the numbering depends on the order the
    // cells are created.
    DoFRenumbering::hierarchical (free_surface_dof_handler);

    mesh_locally_owned = free_surface_dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs (free_surface_dof_handler,
                                             mesh_locally_relevant);

    mesh_displacements.reinit(mesh_locally_owned, mesh_locally_relevant, sim.mpi_communicator);
    fs_mesh_velocity.reinit(mesh_locally_owned, mesh_locally_relevant, sim.mpi_communicator);


    //if we are just starting, we need to initialize the mesh displacement vector.
    
	//Here we can set up initial topography
	if (sim.timestep_number == 0)
	//mesh_displacements = 0.;
	{
		//Interpolate the mesh vertex velocity onto the Stokes velocity system for use in ALE corrections
		LinearAlgebra::Vector distributed_initial_displacements;
		distributed_initial_displacements.reinit(mesh_locally_owned, sim.mpi_communicator);

		const std::vector<Point<dim> > support_points
			= free_surface_fe.base_element(0).get_unit_support_points();

		Quadrature<dim> quad(support_points);
		UpdateFlags update_flags = UpdateFlags(update_quadrature_points);
		FEValues<dim> fe_values (free_surface_fe, quad, update_flags);
		//FEFaceValues<dim> fe_face_values (*sim.mapping, sim.finite_element, quadrature_formula, update_values | update_gradients | update_JxW_values)
		const unsigned int n_q_points = fe_values.n_quadrature_points,
											 dofs_per_cell = fe_values.dofs_per_cell;
		// std::cout<<"n_q_points="<<n_q_points<<"\n";
		// std::cout<<"dim="<<dim<<"\n";
		// std::cout<<"length of support_points="<<support_points.size()<<"\n";

		std::vector<types::global_dof_index> cell_dof_indices (dofs_per_cell);
		FEValuesExtractors::Vector extract_displacement(0);
		std::vector<Tensor<1,dim> > displacement_values(n_q_points);
		
		typename DoFHandler<dim>::active_cell_iterator
		cell = free_surface_dof_handler.begin_active(), 	endc= free_surface_dof_handler.end();
		
		//===
		
			  //added by lev for benchmarking purposes

		// typedef std::vector<xu> displacements_vector;
		// typedef std::vector<displacements_vector> displacements_matrix;
		// static displacements_matrix u;
		// u.push_back(displacements_vector());
		
		
		
		for (; cell!=endc; ++cell)
			if (cell->is_locally_owned())
				{
					cell->get_dof_indices (cell_dof_indices);

					fe_values.reinit (cell);
					for (unsigned int j=0; j<n_q_points; ++j)
					{
						const Point<dim> p = fe_values.quadrature_point(j);
						Point<dim> displacement;
						displacement[0] = 0.;
						displacement[1] = 0.;
						//displacement[dim-1] = p[dim-1]*0.1*std::sin(p[0]*M_PI);
						double L=1; double amp=0.1;
						
						//Zero topography
						displacement[dim-1] = 0;
						
						//^-shaped topography
						//displacement[dim-1] = -p[dim-1]*2*amp*(std::abs(L/2-p[0]))+amp;
						//constant topography
						//displacement[dim-1] = 0;
						//linear topography
						//displacement[dim-1] = p[dim-1]*amp*p[0]; 
						//quadratic topography
						//displacement[dim-1] = p[dim-1]*amp*(4*p[0]-4*p[0]*p[0]); 
						
						//Hammocky topography
						// if (dim==3)
							// displacement[dim-1] = amp+amp*p[2]*std::sin(p[0])*std::sin(p[1]); 
						// else
							// displacement[dim-1] = amp*p[1]*std::sin(p[0]); 
							
						unsigned int support_point_index;
						for (unsigned int dir=0; dir<dim; ++dir)
							{
								
								support_point_index = free_surface_fe.component_to_system_index(/*velocity component=*/ dir,
																															/*dof index within component=*/ j);
								distributed_initial_displacements[cell_dof_indices[support_point_index]] = displacement[dir];
								
							}
						// xu temp;
						// temp.x=p[0];
						// temp.y=displacement[dim-1];
						// u[0].push_back(temp);
					}
				}
		
			
		
		//std::cout<<"indx = "<<indx<<std::endl;
		// std::cout<<"len (u[0]) = "<<u[0].size()<<std::endl;
		
		// std::sort(u[0].begin(), u[0].end(), xuCompare);
		
		// for (unsigned int j=0; j<u[0].size(); ++j) 
			// std::cout<<u[0][j].x<<" ";
		// std::ofstream myfile("displacements.txt", std::ios::app);
		
		// /* fancy way but can't assign xu to string
		// std::ostream_iterator<int> output_iterator(myfile, "\t");
		// myfile << "Displacements, timestep "<<sim.time_step<<std::endl;	
		
		
		// for(typename displacements_matrix::const_iterator it=u.begin(); it!=u.end(); ++it) 
		// {
			// std::copy(it->begin(), it->end(), output_iterator);//std::ostream_iterator<int>(std::cout, " "));
			// std::cout << std::endl;
		// }
		// */
		// //for(unsigned int i=0; i<u.size(); ++i) 
			// for (unsigned int j=0; j<u[0].size(); ++j) 
				// {
					// myfile << u[0][j].x<<" ";//std::ostream_iterator<int>(std::cout, " "));
					// myfile << u[0][j].y<<std::endl;
				// }
		// myfile.close();

		distributed_initial_displacements.compress(VectorOperation::insert);
		mesh_displacements = distributed_initial_displacements;
		
		
		// std::ofstream file_initial("initial", std::ios::app);
		// for (int i =0; i<distributed_initial_displacements.size();i++)
			// file_initial << distributed_initial_displacements[i]<<" ";
		// file_initial << std::endl;
		// file_initial.close();		
		//std::cout<<"Setup_DOFs complete\n";
	}

    // We would like to make sure that the mesh stays conforming upon
    // redistribution, so we construct mesh_vertex_constraints, which
    // keeps track of hanging node constraints.
    // Note: this would be a more natural fit in make_constraints(),
    // but we would like to be able to apply vertex constraints directly
    // after setup_dofs(), as is done, for instance, during mesh
    // refinement.
    mesh_vertex_constraints.clear();
    mesh_vertex_constraints.reinit(mesh_locally_relevant);

    DoFTools::make_hanging_node_constraints(free_surface_dof_handler, mesh_vertex_constraints);

    // We can safely close this now
    mesh_vertex_constraints.close();

    // Now reset the mapping of the simulator to be something that captures mesh deformation in time.
    sim.mapping.reset (new MappingQ1Eulerian<dim, LinearAlgebra::Vector> (free_surface_dof_handler,
                                                                          mesh_displacements));
  }

  template <int dim>
  double FreeSurfaceHandler<dim>::get_stabilization_term() const
  {
    return free_surface_theta;
  }

	/*
  template <int dim>
  void
  FreeSurfaceHandler<dim>::detach_manifolds()
  {
    //remove manifold ids for all cells
    for (typename Triangulation<dim>::active_cell_iterator
         cell = sim.triangulation.begin_active();
         cell != sim.triangulation.end(); ++cell)
      cell->set_all_manifold_ids (numbers::invalid_manifold_id);

    //also remove boundary description for the free surface indicators
    std::set<types::boundary_id>::iterator id = sim.parameters.free_surface_boundary_indicators.begin();
    for (; id != sim.parameters.free_surface_boundary_indicators.end(); ++id)
      sim.triangulation.set_boundary( *id );
  }
	*/
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template class FreeSurfaceHandler<dim>;

  ASPECT_INSTANTIATE(INSTANTIATE)
}
