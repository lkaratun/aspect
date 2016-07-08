/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef __aspect__velocity_boundary_conditions_plate_tectonics_h
#define __aspect__velocity_boundary_conditions_plate_tectonics_h

#include <aspect/velocity_boundary_conditions/interface.h>
#include <aspect/simulator_access.h>



namespace aspect
{
  namespace VelocityBoundaryConditions
  {
    using namespace dealii;

    /**
     * A class that implements velocity boundary conditions based on a
     * functional description provided in the input file.
     *
     * @ingroup VelocityBoundaryConditionsModels
     */
    template <int dim>
    class plate_tectonics : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        //plate_tectonics ();

        /**
         * Return the boundary velocity as a function of position. For the
         * current class, this function obviously simply returns a zero
         * tensor.
         */
        virtual
        Tensor<1,dim>
        boundary_velocity (const Point<dim> &position) const;
        static
        void
        declare_parameters (ParameterHandler &prm);
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        double vel;
        double angle;
        double wz_width;
        double wz_height;
        double wz_depth;
        double width;
        double depth;
        double total_thickness;
        double lithospheric_thickness;
        double transition_zone;
        double alfa;



    };
  }
}


#endif
