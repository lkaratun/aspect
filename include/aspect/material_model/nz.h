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



#ifndef __aspect__model_nz_h
#define __aspect__model_nz_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model that consists of globally constant values for all
     * material parameters except the density, which depends linearly on the
     * temperature. The model is considered incompressible.
     *
     * This material model implements what the "Simple" model was originally
     * intended to do, before it got too complicated.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class nz : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:

        virtual bool
        viscosity_depends_on (const NonlinearDependence::Dependence dependence) const;

        virtual bool
        density_depends_on (const NonlinearDependence::Dependence dependence) const;

        virtual bool
        compressibility_depends_on (const NonlinearDependence::Dependence dependence) const;

        virtual bool
        specific_heat_depends_on (const NonlinearDependence::Dependence dependence) const;

        virtual bool
        thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const;

        virtual bool is_compressible () const;

        virtual double reference_viscosity () const;

        virtual double reference_density () const;

        virtual void evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                              typename Interface<dim>::MaterialModelOutputs &out) const;
															
				virtual double viscosity_ratio (const double temperature,
																				const double pressure,
																				const std::vector<double> &comp,
																				const SymmetricTensor<2,dim> &strain_rate,
																				const Point<dim> &p) const;															


        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */
        /**
         * Declare the parameters this class takes through input files.
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
         * @}
         */

      private:
        enum averaging_scheme
        {
          harmonic,
          arithmetic,
          geometric,
          maximum_composition
        } viscosity_averaging;

        std::vector<double> compute_volume_fractions(const std::vector<double> &compositional_fields) const;
        double average_value (const std::vector<double> &composition,
                              const std::vector<double> &parameter_values,
                              const enum averaging_scheme &average_type) const;
        double reference_T;
        double eta_min;
        double eta_max;
        double eta_reference;
        std::vector<double> angle_if;
        std::vector<double> cohesion;
        double reference_specific_heat;
        double k_value;
        double composition_viscosity_prefactor;
        double thermal_viscosity_exponent;
        double mantle_viscosity;
        double min_strain_rate;
        double ref_strain_rate;
        std::vector<double> densities;
        std::vector<double> thermal_alpha;
        //new
        std::vector<double> activation_energies;
        std::vector<double> activation_volumes;
        std::vector<double> material_parameters;
        std::vector<double> nvs;
        double num_plastic;

    };

  }
}

#endif
