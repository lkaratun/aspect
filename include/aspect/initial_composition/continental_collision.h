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


#ifndef _aspect_initial_composition_continental_collision_h
#define _aspect_initial_composition_continental_collision_h

#include <aspect/initial_composition/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace InitialComposition
  {
    /**
     * A class that implements a composition boundary condition for a
     * box geometry with additional boundary indicators for the
     * lithospheric part of the left and right (and front and back in 3D) boundaries.
     *
     * @ingroup BoundaryCompositions
     */
    template <int dim>
    class ContinentalCollision : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * This function returns constant compositions at the boundaries.
         *
         * @copydoc aspect::BoundaryComposition::Interface::initial_composition()
         */
        virtual
        double initial_composition (const Point<dim> &position,
                                     const unsigned int compositional_field) const;

        /**
         * Declare the parameters this class takes through input files. This
         * class declares the inner and outer boundary compositions.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

        /**
         * This function performs some basic sanity checks on the parameter
         * values previously read from the input file.
         */
        virtual void initialize ();

      private:
        /**
         * The values of the various composition variables on each of the
         * 2*dim+2*(dim-1) boundaries of the box.
         */
        // std::vector<double> composition_values[2*dim+2*(dim-1)];
        //Z direction
        double total_depth;
        //X direction
        double total_width;
        //Y direction
        double total_length;
        double lithospheric_depth;
        double oceanic_depth;
        double continental_depth;
        //y-direction
        double continental_length;
        //x-direction
        double continental_width;

        double wz_width;
        double wz_height;
        double wz_depth;
        double wz_angle;
    };
  }
}


#endif
