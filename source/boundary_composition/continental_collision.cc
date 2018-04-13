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


#include <aspect/boundary_composition/continental_collision.h>

#include <utility>
#include <limits>


namespace aspect
{
  namespace BoundaryComposition
  {

    template <int dim>
    double
    ContinentalCollision<dim>::
    boundary_composition (const Point<dim> &p/*position*/,
                          const unsigned int compositional_field) const
    {
    // #set Function constants  = total_depth=512e3, total_width=1024e3, total_length = 1024e3, crustal_depth=20e3, lower_crustal_depth=30e3, lithospheric_depth=120e3, wz_width=20.e3, wz_height=120.e3, wz_depth=60e3, wz_angle =  60, pi=3.14159265359 #30*3.141592654/180 # angle from vertical
    // set Function constants  = total_depth=512e3, total_width=512e3, total_length = 512e3, oceanic_depth=15e3, continental_depth=30e3, continental_length = 300e3, oceanic_width = 112e3, continental_width = 300e3, lithospheric_depth=120e3, wz_width=20.e3, wz_height=40.e3, wz_depth=60e3
    // #Width : x, Depth: z, length: y



    // #Variable angle weak zone
    // #set Function expression = if (z > total_depth-crustal_depth && (z < total_depth-wz_depth-wz_height/2 || z > total_depth-wz_depth+wz_height/2 || (abs(x-total_width/2-(total_depth-wz_depth-z)/tan(pi/2-wz_angle*pi/180*(total_length/2-y)/(total_length/2)))>wz_width/2)), 1, 0);
    // if (z > total_depth-lower_crustal_depth && z <= total_depth-crustal_depth && (z < total_depth-wz_depth-wz_height/2 || z > total_depth-wz_depth+wz_height/2 || (abs(x-total_width/2-(total_depth-wz_depth-z)/tan(pi/2-wz_angle*pi/180*(total_length/2-y)/(total_length/2)))>wz_width/2)), 1, 0);
    // if (z >= total_depth-wz_depth-wz_height/2 && z <= total_depth-wz_depth+wz_height/2 && (abs(x-total_width/2-(total_depth-wz_depth-z)/tan(pi/2-wz_angle*pi/180*(total_length/2-y)/(total_length/2)))<=wz_width/2), 1, 0);
    // if (z > total_depth-lithospheric_depth && z <= total_depth-lower_crustal_depth && (z < total_depth-wz_depth-wz_height/2 || z > total_depth-wz_depth+wz_height/2 || (abs(x-total_width/2-(total_depth-wz_depth-z)/tan(pi/2-wz_angle*pi/180*(total_length/2-y)/(total_length/2)))>wz_width/2)), 1, 0);    if (z <= total_depth-lithospheric_depth, 1, 0)

      const double x=p(0);
      const double y=p(1);
      //z is always the vertical direction, regardless of the number of dimensions
      const double z=p(dim-1);

      //This will be the number of our compositional field at point p. Initialize as Sublithospheric mantle.
      //Each following c.f. cuts through the previous ones
      unsigned int resulting_cf = 4;

      double total_depth = 512e3;
      double total_width = 512e3;
      double total_length = 512e3;
      double lithospheric_depth = 120e3;
      double oceanic_depth = 30e3;
      double continental_depth = 60e3;
      //y-direction
      double continental_length = 300e3;
      //x-direction
      double continental_width = 300e3;

      double wz_width=40e3;
      double wz_height=80e3;
      double wz_depth=80e3;
      double wz_angle = 60;

      //Lithospheric mantle
      if (z > total_depth-lithospheric_depth) resulting_cf = 3;

      //Oceanic crust
      if (z > total_depth-oceanic_depth) resulting_cf = 1;

      //Continental crust
      if (dim == 3)
        if (z > total_depth-continental_depth
            && abs(x - (total_width/2)) <= continental_width/2
            && abs(y - (total_length/2)) <= continental_length/2)
          resulting_cf = 0;
      else
        if (z > total_depth-continental_depth
            && abs(x - (total_width/2)) <= continental_width/2)
          resulting_cf = 0;

      //Weak zone
      if (z >= total_depth-wz_depth-wz_height/2
          && z <= total_depth-wz_depth+wz_height/2
          && (abs(x-total_width/2-(total_depth-wz_depth-z)/std::tan(M_PI/2-wz_angle*M_PI/180*(total_length/2-y)/(total_length/2)))<=wz_width/2))
          resulting_cf = 2;

      if (resulting_cf == compositional_field) return 1;
      else return 0;



    } //boundary_composition


    template <int dim>
    void
    ContinentalCollision<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary composition model");
      {
        prm.enter_subsection("Continental collision");
        {
          prm.declare_entry ("Left composition", "",
                             Patterns::List(Patterns::Double ()),
                             "A comma separated list of composition boundary values "
                             "at the left boundary (at minimal x-value). This list must have as many "
                             "entries as there are compositional fields. Units: none.");
          prm.declare_entry ("Right composition", "",
                             Patterns::List(Patterns::Double ()),
                             "A comma separated list of composition boundary values "
                             "at the right boundary (at maximal x-value). This list must have as many "
                             "entries as there are compositional fields. Units: none.");
          prm.declare_entry ("Left composition lithosphere", "",
                             Patterns::List(Patterns::Double ()),
                             "A comma separated list of composition boundary values "
                             "at the left boundary (at minimal x-value). This list must have as many "
                             "entries as there are compositional fields. Units: none.");
          prm.declare_entry ("Right composition lithosphere", "",
                             Patterns::List(Patterns::Double ()),
                             "A comma separated list of composition boundary values "
                             "at the right boundary (at maximal x-value). This list must have as many "
                             "entries as there are compositional fields. Units: none.");
          prm.declare_entry ("Bottom composition", "",
                             Patterns::List(Patterns::Double ()),
                             "A comma separated list of composition boundary values "
                             "at the bottom boundary (at minimal y-value in 2d, or minimal "
                             "z-value in 3d). This list must have as many "
                             "entries as there are compositional fields. Units: none.");
          prm.declare_entry ("Top composition", "",
                             Patterns::List(Patterns::Double ()),
                             "A comma separated list of composition boundary values "
                             "at the top boundary (at maximal y-value in 2d, or maximal "
                             "z-value in 3d). This list must have as many "
                             "entries as there are compositional fields. Units: none.");

        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    ContinentalCollision<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary composition model");
      {
        prm.enter_subsection("Box with lithosphere boundary indicators");
        {
          switch (dim)
            {
              case 2:
                composition_values[0] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Left composition")));
                composition_values[1] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Right composition")));
                composition_values[2] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Bottom composition")));
                composition_values[3] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Top composition")));
                composition_values[4] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Left composition lithosphere")));
                composition_values[5] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Right composition lithosphere")));
                break;

              case 3:
                composition_values[0] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Left composition")));
                composition_values[1] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Right composition")));
                composition_values[2] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Front composition")));
                composition_values[3] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Back composition")));
                composition_values[4] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Bottom composition")));
                composition_values[5] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Top composition")));
                composition_values[6] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Left composition lithosphere")));
                composition_values[7] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Right composition lithosphere")));
                composition_values[8] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Front composition lithosphere")));
                composition_values[9] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Back composition lithosphere")));
                break;

              default:
                Assert (false, ExcNotImplemented());
            }
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }



    template <int dim>
    void
    ContinentalCollision<dim>::initialize()
    {
      // verify that each of the lists for boundary values
      // has the requisite number of elements
      for (unsigned int f=0; f<2*dim+2*(dim-1); ++f)
        AssertThrow (composition_values[f].size() == this->n_compositional_fields(),
                     ExcMessage (std::string("The specification of boundary composition values for the `box with "
                                             "lithosphere boundary indicators' model "
                                             "requires as many values for each boundary indicator the box as there "
                                             "are compositional "
                                             "fields. However, for boundary indicator ")
                                 +
                                 Utilities::int_to_string(f)
                                 +
                                 ", the input file specifies "
                                 +
                                 Utilities::int_to_string(composition_values[f].size())
                                 +
                                 " values even though there are "
                                 +
                                 Utilities::int_to_string(this->n_compositional_fields())
                                 +
                                 " compositional fields."));
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryComposition
  {
    ASPECT_REGISTER_BOUNDARY_COMPOSITION_MODEL(ContinentalCollision,
                                               "Continental collision",
                                               "A model in which the composition is chosen constant on "
                                               "all the sides of a box. Additional boundary indicators "
                                               "are added to the lithospheric parts of the vertical boundaries. "
                                               "This model is to be used with the 'Two Merged Boxes' Geometry Model.")
  }
}
