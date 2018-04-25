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


#include <aspect/initial_composition/continental_collision.h>

#include <utility>
#include <limits>


namespace aspect
{
  namespace InitialComposition
  {

    template <int dim>
    double
    ContinentalCollision<dim>::
    initial_composition (const Point<dim> &p/*position*/,
                          const unsigned int compositional_field) const
    {
      if (compositional_field<0 || compositional_field>4) std::cout << "compositional_field = " << compositional_field;

      const double x=p(0);
      const double y=p(1);
      //z is always the vertical direction, regardless of the number of dimensions
      const double z=p(dim-1);

      //This will be the number of our compositional field at point p. Initialize as Sublithospheric mantle.
      //Each following c.f. cuts through the previous ones
      unsigned int resulting_cf = 4;


      //Lithospheric mantle
      if (z >= total_depth-lithospheric_depth) resulting_cf = 3;

      //Oceanic crust
      if (z >= total_depth-oceanic_depth) resulting_cf = 1;

      //Continental crust
      if (dim == 3)
      {
        if (configuration == "continent centered")
        {
          if (z >= total_depth-continental_depth
              && abs(x - (total_width/2)) <= continental_width/2
              && abs(y - (total_length/2)) <= continental_length/2)
            resulting_cf = 0;
        }
        else
        {
          if (z >= total_depth-continental_depth
              && abs(x - (total_width/2)) <= continental_width/2
              && abs(y - (total_length/2)) >= continental_length/2)
            resulting_cf = 0;
        }
      }
      else
        if (z >= total_depth-continental_depth
            && abs(x - (total_width/2)) <= continental_width/2)
          resulting_cf = 0;

      //Rotated weak zone
      // if (z >= total_depth-wz_depth-wz_height/2
      //     && z <= total_depth-wz_depth+wz_height/2
      //     && (abs(x-total_width/2-(total_depth-wz_depth-z)/std::tan(M_PI/2-wz_angle*M_PI/180*(total_length/2-y)/(total_length/2)))<=wz_width/2))
      //     resulting_cf = 2;

      //Planar weak zone
      if (z >= total_depth-wz_depth-wz_height/2
          && z <= total_depth-wz_depth+wz_height/2
          && (abs(x-total_width/2-(total_depth-wz_depth-z)/std::tan(M_PI/2-wz_angle*M_PI/180))<=wz_width/2))
          resulting_cf = 2;


      if (resulting_cf == compositional_field) return 1;
      else return 0;



    } //initial_composition


    template <int dim>
    void
    ContinentalCollision<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Continental collision");
        {
          prm.declare_entry ("X extent", "512e3",
                             Patterns::Double (),
                             "Length of the model domain along the X axis. Units: m.");
          prm.declare_entry ("Y extent", "512e3",
                             Patterns::Double (),
                             "Length of the model domain along the Y axis. Units: m.");
          prm.declare_entry ("Z extent", "512e3",
                             Patterns::Double (),
                             "Length of the model domain along the Z axis. Units: m.");
          prm.declare_entry ("Continental depth", "30e3",
                             Patterns::Double (),
                             "Depth of the bottom boundary of the continental lithosphere. Units: m.");
          prm.declare_entry ("Continental length", "256e3",
                             Patterns::Double (),
                             "Length of the bottom boundary of the continental lithosphere. Units: m.");
          prm.declare_entry ("Continental width", "256e3",
                             Patterns::Double (),
                             "Width of the bottom boundary of the continental lithosphere. Units: m.");
          prm.declare_entry ("Oceanic depth", "15e3",
                             Patterns::Double (),
                             "Depth of the bottom boundary of the oceanic lithosphere. Units: m.");
          prm.declare_entry ("Lithospheric depth", "120e3",
                             Patterns::Double (),
                             "Depth of the bottom boundary of the mantle lithosphere. Units: m.");
          prm.declare_entry ("Weak zone depth", "75e3",
                             Patterns::Double (),
                             "Depth of the middle point of the weak zone. Units: m.");
          prm.declare_entry ("Weak zone height", "40e3",
                             Patterns::Double (),
                             "Vertical extent of the weak zone. Units: m.");
          prm.declare_entry ("Weak zone width", "20e3",
                             Patterns::Double (),
                             "Horizontal extent of the weak zone. Units: m.");
          prm.declare_entry ("Weak zone angle", "60",
                             Patterns::Double (),
                             "Angle of the weak zone. Units: degrees.");
          prm.declare_entry ("Configuration", "continent centered",
                         Patterns::Selection("continent centered|continent aside"),
                         "Choose between 2 geometries. Doesn't have effect in 2D");

        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    ContinentalCollision<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Continental collision");
        {
          total_width = Utilities::string_to_double(prm.get ("X extent"));
          total_length = Utilities::string_to_double(prm.get ("Y extent"));
          total_depth = Utilities::string_to_double(prm.get ("Z extent"));
          continental_width = Utilities::string_to_double(prm.get ("Continental width"));
          continental_length = Utilities::string_to_double(prm.get ("Continental length"));
          continental_depth = Utilities::string_to_double(prm.get ("Continental depth"));
          oceanic_depth = Utilities::string_to_double(prm.get ("Oceanic depth"));
          lithospheric_depth = Utilities::string_to_double(prm.get ("Lithospheric depth"));
          wz_width = Utilities::string_to_double(prm.get ("Weak zone width"));
          wz_height = Utilities::string_to_double(prm.get ("Weak zone height"));
          wz_depth = Utilities::string_to_double(prm.get ("Weak zone depth"));
          wz_angle = Utilities::string_to_double(prm.get ("Weak zone angle"));
          configuration = prm.get ("Configuration");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }



    template <int dim>
    void
    ContinentalCollision<dim>::initialize()
    {

    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialComposition
  {
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(ContinentalCollision,
                                               "Continental collision",
                                               "A model in which the composition is chosen constant on "
                                               "all the sides of a box. Additional boundary indicators "
                                               "are added to the lithospheric parts of the vertical boundaries. "
                                               "This model is to be used with the 'Two Merged Boxes' Geometry Model.")
  }
}
