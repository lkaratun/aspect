/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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

//#include <lithosphere.h>
#include <aspect/material_model/nz.h>
#include <deal.II/base/parameter_handler.h>

#include <math.h>
#include <boost/math/special_functions/fpclassify.hpp> // isnan



using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {
		template <int dim>
    std::vector<double>
    nz<dim>::
    compute_volume_fractions(const std::vector<double> &compositional_fields) const
    {
      std::vector<double> volume_fractions(compositional_fields.size());

      //clip the compositional fields so they are between zero and one
      std::vector<double> x_comp = compositional_fields;
      for ( unsigned int i=0; i < x_comp.size(); ++i)
        x_comp[i] = std::min(std::max(x_comp[i], 0.0), 1.0);

      //sum the compositional fields for normalization purposes
      double sum_composition = 0.0;
      for ( unsigned int i=0; i < x_comp.size(); ++i)
        sum_composition += x_comp[i];

      if (sum_composition)
        for ( unsigned int i=0; i < x_comp.size(); ++i)
          volume_fractions[i] = x_comp[i]/sum_composition;
      else
        for ( unsigned int i=0; i < x_comp.size(); ++i)
          volume_fractions[i]=1/x_comp.size();
      return volume_fractions;
    }

    template <int dim>
    double
    nz<dim>::
    average_value ( const std::vector<double> &composition,
                    const std::vector<double> &parameter_values,
                    const enum averaging_scheme &average_type) const
    {
      double averaged_parameter = 0.0;
      const std::vector<double> volume_fractions = compute_volume_fractions(composition);

      switch (average_type)
        {
          case arithmetic:
          {
            for (unsigned int i=0; i< volume_fractions.size(); ++i)
              averaged_parameter += volume_fractions[i]*parameter_values[i];
            break;
          }
          case harmonic:
          {
            for (unsigned int i=0; i< volume_fractions.size(); ++i)
              averaged_parameter += volume_fractions[i]/(parameter_values[i]);
            averaged_parameter = 1.0/averaged_parameter;
            break;
          }
          case geometric:
          {
            for (unsigned int i=0; i < volume_fractions.size(); ++i)
              averaged_parameter += volume_fractions[i]*std::log(parameter_values[i]);
            averaged_parameter = std::exp(averaged_parameter);
            break;
          }
          case maximum_composition:
          {
            const unsigned int i = (unsigned int)(std::max_element( volume_fractions.begin(),
                                                                    volume_fractions.end() )
                                                  - volume_fractions.begin());
            averaged_parameter = parameter_values[i];
            break;
          }
          default:
          {
            AssertThrow( false, ExcNotImplemented() );
            break;
          }
        }
      return averaged_parameter;
    }


    template <int dim>
    void
    nz<dim>::
    evaluate(const typename Interface<dim>::MaterialModelInputs &in,
             typename Interface<dim>::MaterialModelOutputs &out) const
    {

      //Gas constant
      const double R = 8.314;


      bool crash = false;

      for (unsigned int i=0; i<in.position.size(); ++i)
        {
          const std::vector<double> composition = in.composition[i];
          const std::vector<double> volume_fractions = compute_volume_fractions(composition);
					

          std::vector<double> viscosities(volume_fractions.size());

          //double strainrate_E2 = sqrt(0.5*(in.strain_rate[i][0][0]*in.strain_rate[i][0][0]+in.strain_rate[i][1][1]*in.strain_rate[i][1][1]+in.strain_rate[i][2][2]*in.strain_rate[i][2][2]) + in.strain_rate[i][0][1]*in.strain_rate[i][0][1] + in.strain_rate[i][0][2]*in.strain_rate[i][0][2] + in.strain_rate[i][1][2]*in.strain_rate[i][1][2]);

          //IIE=s11s22+s22s33+s33s11−s12^2−s13^2−s23^2
          double strainrate_E2 = 0;
          //Textbook formulation
          //if (in.strain_rate.size())
          //strainrate_E2 = in.strain_rate[i][0][0]*in.strain_rate[i][1][1]+in.strain_rate[i][1][1]*in.strain_rate[i][2][2]+in.strain_rate[i][2][2]*in.strain_rate[i][0][0]-in.strain_rate[i][0][1]*in.strain_rate[i][0][1]-in.strain_rate[i][0][2]*in.strain_rate[i][0][2]-in.strain_rate[i][1][2]*in.strain_rate[i][1][2];




          if (in.strain_rate.size())
            strainrate_E2 = in.strain_rate[i].norm();

          strainrate_E2 = std::max(strainrate_E2,min_strain_rate);
          //if (strainrate_E2>1e-18)
          //std::cout<<"strainrate_E2: "<<strainrate_E2<<"  ";

          double sigma_y, viscosity_MC, strainrate_mod, viscosity_dislocation_creep, total_viscosity;
          //double h=1e-3;
          //double m=0;


          for (unsigned int j=0; j < volume_fractions.size(); j++)
            {

              //***************Calculate crustal viscosity***************
              //Drucker-Prager yield criterion (by John)
              sigma_y = ( (dim==3)
                          ?
                          ( 6.0 * cohesion[j] * std::cos(angle_if[j]*M_PI/180) + 2.0 * std::max(in.pressure[i],0.0) * std::sin(angle_if[j]*M_PI/180) )
                          / ( std::sqrt(3.0) * (3.0 + std::sin(angle_if[j]*M_PI/180) ) )
                          :
                          cohesion[j] * std::cos(angle_if[j]*M_PI/180) + std::max(in.pressure[i],0.0) * std::sin(angle_if[j]*M_PI/180) );


              //Mohr-Coulomb criterion
              //sigma_y = cohesion[j] + in.pressure[i]*std::tan(angle_if[j]*M_PI/180);
			if (strainrate_E2)
			{
				strainrate_mod=strainrate_E2/ref_strain_rate;
				strainrate_mod=pow(strainrate_mod,(nps[j]-1)/nps[j]);
				strainrate_mod*=ref_strain_rate;
				viscosity_MC = 1/(1/(sigma_y/(2*strainrate_mod)+eta_min)+1/eta_max);
			}
			else
				viscosity_MC=eta_max;

              //Calculate dislocation creep

              //Gerya formula
              double F2 = 1/pow(2,(2*nvs[j]-1)/nvs[j]);
              //F2 = (pow(3, -(1+n_v)/(2*n_v))*pow(2, (1-n_v)/n_v));

              if (strainrate_E2 && in.temperature[i])
                {
                  //Normalize strain rate
                  strainrate_mod=strainrate_E2/ref_strain_rate;
                  strainrate_mod=pow(strainrate_mod,(nvs[j]-1)/nvs[j]);
				  strainrate_mod*=ref_strain_rate;
                  strainrate_mod=1/strainrate_mod;
				  
                  //if (j<3) //crust and weak zone
                  if (in.temperature[i])
                    //Gerya formula
                    viscosity_dislocation_creep =  F2 * pow(material_parameters[j],-1/nas[j]) * strainrate_mod *
                                                   exp((activation_energies[j]+activation_volumes[j]*std::max(in.pressure[i],0.0))/(nts[j]*R*in.temperature[i]));

                  else viscosity_dislocation_creep=eta_max;

                  //const double depth = this->get_geometry_model().depth(in.position[i]);
                  //if (depth<100000 && strainrate_E2>1e-18 && j==1)// && volume_fractions[1]>0.7)//&& in.strain_rate.size() && composition[1]>0.7

                  /*
                  const double x=in.position[i](0);
                  const double y=in.position[i](1);
                  const double z=in.position[i](2);

                  if (in.temperature[i]<1000 && in.temperature[i]>900 && x> 100e3 && x<200e3 && j==3)
                  {
                    // std::cout<<"depth: "<<depth<<"  ";
                    std::cout<<"strainrate_E2: "<<strainrate_E2<<"  ";
					std::cout<<"nvs[j]: "<<nvs[j]<<"  ";
					std::cout<<"strainrate_mod: "<<strainrate_mod<<"  ";
					
                    // std::cout<<"std::tan(angle_if[1]*M_PI/180): "<<std::tan(angle_if[1]*M_PI/180)<<"  ";
                    //std::cout<<"sigma_y: "<<sigma_y<<"  ";
                    //std::cout<<"viscosity_MC: "<<viscosity_MC<<"  ";
                    std::cout<<"temperature: "<<in.temperature[i]<<"  ";
                    std::cout<<"pressure: "<<in.pressure[i]<<"  ";
                    std::cout<<"exp: "<<exp((activation_energies[j]+activation_volumes[j]*in.pressure[i])/(nts[j]*R*in.temperature[i]))<<"  ";
                    //std::cout<<"prefactor: "<<F2 * 1/pow(material_parameters[j],1/nvs[j])<<"  ";
					std::cout<<"preexp term: "<<pow(material_parameters[j],-1/nas[j])<<"  ";
                    std::cout<<"viscosity_dislocation_creep: "<<viscosity_dislocation_creep<<"  ";
                  }
					*/




                }
              else viscosity_dislocation_creep=eta_max;



              //viscosity_dislocation_creep=std::max(std::min(viscosity_dislocation_creep,eta_max),eta_min);
              //viscosity_MC=std::max(std::min(viscosity_MC,eta_max),eta_min);

              //total_viscosity = pow(10,(log10(viscosity_MC)+log10(viscosity_dislocation_creep))/2);


              //double total_viscosity = 2*viscosity_dislocation_creep*viscosity_MC/(viscosity_dislocation_creep+viscosity_MC);

              if (!crash)
                {
					if (boost::math::isnan(viscosity_dislocation_creep) || boost::math::isnan(viscosity_MC))
                    {
                      std::cout<<"viscosity_dislocation_creep or MC is NAN. Visc_disl_creep = "<<viscosity_dislocation_creep<<"  Visc_MC = "<<viscosity_MC<<" \n";
						std::cout<<"strainrate_E2: "<<strainrate_E2<<"  ";
						std::cout<<"strainrate_mod: "<<strainrate_mod<<"  \n";
						std::cout<<"material_parameters[j]: "<<material_parameters[j]<<"  ";
						std::cout<<"F2: "<<F2<<"  ";
						std::cout<<"nas[j]: "<<nas[j]<<"  ";
						std::cout<<"nvs[j]: "<<nvs[j]<<"  ";
						std::cout<<"temperature: "<<in.temperature[i]<<"  ";
						std::cout<<"pressure: "<<in.pressure[i]<<"  \n";
						std::cout<<"prefactor: "<<pow(material_parameters[j],-1/nas[j]) * strainrate_mod<<"  ";
						std::cout<<"exp: "<<exp((activation_energies[j]+activation_volumes[j]*std::max(in.pressure[i],0.0))/(nts[j]*R*in.temperature[i]))<<"  \n\n";
						
                      crash = true;											
					}
                }

              
			  if (boost::math::isnan(viscosity_dislocation_creep) || boost::math::isnan(viscosity_MC))
				  total_viscosity=eta_min;
			  else
				  total_viscosity = std::min(viscosity_dislocation_creep,viscosity_MC);
              viscosities[j] = std::max(std::min(total_viscosity,eta_max),eta_min);

              // if (!crash)
                // {
                  // if (viscosities[j]>=eta_min && viscosities[j]<=eta_max)
                    // {}
                  // else
                    // {
                      // std::cout<<"viscosities[j]: "<<viscosities[j]<<"  ";
                      // crash = true;
                    // }


                // }							

            }






          //***************Calculate total viscosity***************


          double viscosity = 0;
          /*
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
              viscosity += log10(viscosities[j])*volume_fractions[j];

          viscosity = pow(10,viscosity);*/

          viscosity = average_value(composition, viscosities, viscosity_averaging);

          // const double x=in.position[i](0);
          // const double y=in.position[i](1);
          // if (in.temperature[i]<1250 && x> 155e3 && volume_fractions[3]>0.99 && viscosity < 1e22 && x>150e3&&x<250e3 && y<50e3)
            // {

              // std::cout<<"u_crust: "<<volume_fractions[0]<<"  ";
              // std::cout<<"l_crust: "<<volume_fractions[1]<<"  ";
              // std::cout<<"weak_zone: "<<volume_fractions[2]<<"  ";
              // std::cout<<"lith: "<<volume_fractions[3]<<"  ";
              // std::cout<<"sublith: "<<volume_fractions[4]<<"  ";
              // std::cout<<"strainrate_E2: "<<strainrate_E2<<"  ";
              // // std::cout<<"std::tan(angle_if[1]*M_PI/180): "<<std::tan(angle_if[1]*M_PI/180)<<"  ";
              // std::cout<<"sigma_y: "<<sigma_y<<"  ";
              // std::cout<<"viscosity_MC: "<<viscosity_MC<<"  ";
              // std::cout<<"temperature: "<<in.temperature[i]<<"  ";
              // std::cout<<"pressure: "<<in.pressure[i]<<"  ";
              // std::cout<<"exp: "<<exp((activation_energies[3]+activation_volumes[3]*in.pressure[i])/(nvs[3]*R*in.temperature[3]))<<"  ";
              // std::cout<<"prefactor: "<<1/pow(material_parameters[3],1/nvs[3])<<"  ";
              // std::cout<<"viscosity_dislocation_creep: "<<viscosity_dislocation_creep<<"  ";
              // std::cout<<"final visc: "<<viscosity<<"  ";
            // }

          //Crop viscosity
          out.viscosities[i] = std::max(std::min(viscosity,eta_max),eta_min);


          if (out.viscosities[i]>=eta_min && out.viscosities[i]<=eta_max)
            {}
          else
            {
              std::cout<<"out.viscosities[i]: "<<out.viscosities[i]<<"  ";
            }

          double density = 0.0;
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            {
              //not strictly correct if thermal expansivities are different, since we are interpreting
              //these compositions as volume fractions, but the error introduced should not be too bad.
              const double temperature_factor = (1.0 - thermal_alpha[j] * (in.temperature[i] - reference_T[j]));
			  const double pressure_factor = std::exp(reference_compressibility * (in.pressure[i] - this->get_surface_pressure()));
              density += volume_fractions[j] * densities[j] * temperature_factor * pressure_factor;
            }
          out.densities[i] = density;


          double thermal_expansivity = 0.0;
		  for (unsigned int j=0; j < volume_fractions.size(); ++j)
            thermal_expansivity += volume_fractions[j] * thermal_alpha[j];
          out.thermal_expansion_coefficients[i] = thermal_expansivity;
          out.specific_heat[i] = reference_specific_heat;
		  if (t_conductivity)
			  out.thermal_conductivities[i] = t_conductivity;
		  else
			  out.thermal_conductivities[i] = t_diffusivity*reference_specific_heat*density;
          out.compressibilities[i] = reference_compressibility;
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;
        }
    }

    template <int dim>
    double
    nz<dim>::
    viscosity_ratio (const double temperature,
                     const double pressure,
                     const std::vector<double> &composition,
                     const SymmetricTensor<2,dim> &strain_rate,
                     const Point<dim> &p) const
    {
      const std::vector<double> volume_fractions = compute_volume_fractions(composition);
			std::vector<double> viscosities_MC(volume_fractions.size());
			std::vector<double> viscosities_dislocation_creep(volume_fractions.size());
			const double R = 8.314;
			
			double strainrate_E2 = 0;
			strainrate_E2 = strain_rate.norm();
			strainrate_E2 = std::max(strainrate_E2,min_strain_rate);

			double sigma_y, viscosity_MC, strainrate_mod, viscosity_dislocation_creep;

			for (unsigned int j=0; j < volume_fractions.size(); j++)
				{

					//***************Calculate crustal viscosity***************
					//Drucker-Prager yield criterion 
					sigma_y = ( (dim==3)
											?
											( 6.0 * cohesion[j] * std::cos(angle_if[j]*M_PI/180) + 2.0 * std::max(pressure,0.0) * std::sin(angle_if[j]*M_PI/180) )
											/ ( std::sqrt(3.0) * (3.0 + std::sin(angle_if[j]*M_PI/180) ) )
											:
											cohesion[j] * std::cos(angle_if[j]*M_PI/180) + std::max(pressure,0.0) * std::sin(angle_if[j]*M_PI/180) );

					if (strainrate_E2)
					{
						strainrate_mod=strainrate_E2/ref_strain_rate;
						strainrate_mod=pow(strainrate_mod,(nps[j]-1)/nps[j]);
						strainrate_mod*=ref_strain_rate;
						viscosity_MC = 1/(1/(sigma_y/(2*strainrate_mod)+eta_min)+1/eta_max);
					}
					else
						viscosity_MC=eta_max;
					

					//Calculate dislocation creep
					//Gerya formula
					double F2 = 1/pow(2,(2*nvs[j]-1)/nvs[j]);

					if (strainrate_E2 && temperature)
						{
							//Normalize strain rate
							strainrate_mod=strainrate_E2/ref_strain_rate;
							strainrate_mod=pow(strainrate_mod,(nvs[j]-1)/nvs[j]);
							strainrate_mod*=ref_strain_rate;
							strainrate_mod=1/strainrate_mod;
							//Gerya formula
							viscosity_dislocation_creep =  F2 * pow(material_parameters[j],-1/nas[j]) * strainrate_mod *
																						 exp((activation_energies[j]+activation_volumes[j]*std::max(pressure,0.0))/(nts[j]*R*temperature));
						}
					else viscosity_dislocation_creep=eta_max;
					
					
					// total_viscosity = std::min(viscosity_dislocation_creep,viscosity_MC);
					// viscosities[j] = std::max(std::min(total_viscosity,eta_max),eta_min);
					
					viscosities_MC[j]=std::max(std::min(viscosity_MC,eta_max),eta_min);
					viscosities_dislocation_creep[j]=std::max(std::min(viscosity_dislocation_creep,eta_max),eta_min);
				}

			//***************Calculate total viscosity***************
			viscosity_MC = average_value(composition, viscosities_MC, viscosity_averaging);
			//viscosity_MC = eta_max;
			viscosity_dislocation_creep = average_value(composition, viscosities_dislocation_creep, viscosity_averaging);	
			
			return viscosity_MC/viscosity_dislocation_creep; // High return value = viscous behaviour
    }

    template <int dim>
    double
    nz<dim>::
    reference_viscosity () const
    {
      return eta_reference;
    }


    template <int dim>
    double
    nz<dim>::
    reference_density () const
    {
      return densities[0];
    }

    template <int dim>
    bool
    nz<dim>::
    viscosity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      if ((dependence & NonlinearDependence::strain_rate) != NonlinearDependence::none)
        return true;
      if ((dependence & NonlinearDependence::temperature) != NonlinearDependence::none)
        return true;
      if ((dependence & NonlinearDependence::pressure) != NonlinearDependence::none)
        return true;
      if (((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none))
        //&&                (viscosity_difference1 != 0))
        return true;

      // if ((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none)
      // return true;
      return false;
    }

    template <int dim>
    bool
    nz<dim>::
    density_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      // compare this with the implementation of the density() function
      // to see the dependencies
      if ((dependence & NonlinearDependence::temperature) != NonlinearDependence::none)
        return true;
      return false;
    }


    template <int dim>
    bool
    nz<dim>::
    compressibility_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    nz<dim>::
    specific_heat_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    nz<dim>::
    thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      return false;
    }


    template <int dim>
    bool
    nz<dim>::
    is_compressible () const
    {
      return (reference_compressibility != 0);
    }



    template <int dim>
    void
    nz<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("nz model");
        {
          prm.declare_entry ("Densities", "3300.",
                             Patterns::List(Patterns::Double(0)),
                             "List of densities, $\\rho$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: $kg / m^3$");

          prm.declare_entry ("Reference temperature", "293",
                             Patterns::List(Patterns::Double(0)),
                             "The reference temperature $T_0$. The reference temperature is used "
                             "in both the density and viscosity formulas. Units: $K$.");
          prm.declare_entry ("Minimum viscosity", "1e18",
                             Patterns::Double (0),
                             "The minimum possible viscosity $\\eta_0$ value. This viscosity may be "
                             "modified by both temperature and compositional dependencies. Units: $kg/m/s$.");
          prm.declare_entry ("Maximum viscosity", "1e25",
                             Patterns::Double (0),
                             "The maximum possible viscosity $\\eta_0$ value. This viscosity may be "
                             "modified by both temperature and compositional dependencies. Units: $kg/m/s$.");
          prm.declare_entry ("Reference viscosity", "1e23",
                             Patterns::Double (0),
                             "The reference viscosity $\\eta_0$ value. Units: $kg/m/s$.");
          prm.declare_entry ("Angle of internal friction", "15",
                             Patterns::List(Patterns::Double(0)),
                             "Angle of internal friction: used in the brittle failure criterion. Units: degrees.");
          prm.declare_entry ("Cohesion", "1e7",
                             Patterns::List(Patterns::Double(0)),
                             "Cohesion: used in the brittle failure criterion. Units: Pa.");
          prm.declare_entry ("Composition viscosity prefactor", "1.0",
                             Patterns::Double (0),
                             "A linear dependency of viscosity on the first compositional field. "
                             "Dimensionless prefactor. With a value of 1.0 (the default) the "
                             "viscosity does not depend on the composition. See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\xi$ there.");
          prm.declare_entry ("Thermal viscosity exponent", "0.0",
                             Patterns::Double (0),
                             "The temperature dependence of viscosity. Dimensionless exponent. "
                             "See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\beta$ there.");
          prm.declare_entry ("Thermal conductivity", "2.7",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$. "
                             "Units: $W/m/K$.");
		  prm.declare_entry ("Thermal diffusivity", "0.8e-6",
							 Patterns::Double (0),
							 "The value of the thermal conductivity $k$. "
							 "Units: $W/m/K$.");
          prm.declare_entry ("Reference specific heat", "1250",
                             Patterns::Double (0),
                             "The value of the specific heat $cp$. "
                             "Units: $J/kg/K$.");
          prm.declare_entry ("Thermal expansivities", "3.5e-5",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal expansivities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: $1 / K$");
          prm.declare_entry ("Mantle lithosphere density", "3300",
                             Patterns::Double(),
                             "If compositional fields are used, then one would frequently want "
                             "to make the density depend on these fields. In this simple material "
                             "model, we make the following assumptions: if no compositional fields "
                             "are used in the current simulation, then the density is simply the usual "
                             "one with its linear dependence on the temperature. If there are compositional "
                             "fields, then the density only depends on the first one in such a way that "
                             "the density has an additional term of the kind $+\\Delta \\rho \\; c_1(\\mathbf x)$. "
                             "This parameter describes the value of $\\Delta \\rho$. Units: $kg/m^3/\\textrm{unit "
                             "change in composition}$.");

          prm.declare_entry ("Mantle viscosity", "1e20",
                             Patterns::Double(),
                             "Mantle layer viscosity");


          //new
          prm.declare_entry ("Activation energies", "500",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation energies, $E_a$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: $J / mol$");
          prm.declare_entry ("Activation volumes", "1e-5",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation volumes, $E_v$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: $J / mol$");

          prm.declare_entry ("Stress exponents", "3.5",
                             Patterns::List(Patterns::Double(0)),
                             "List of stress exponents, $n_p$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: None");
          prm.declare_entry ("Temperature exponents", "3.5",
                             Patterns::List(Patterns::Double(0)),
                             "List of temerature exponents, $n_t$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: None");														 
          prm.declare_entry ("Plastic exponents", "10",
                             Patterns::List(Patterns::Double(0)),
                             "List of plastic stress exponents, $n_p$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: None");														 
          prm.declare_entry ("Material exponents", "3.5",
                             Patterns::List(Patterns::Double(0)),
                             "List of material stress exponents, $n_p$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: None");
							 
          prm.declare_entry ("Material parameters", "1e-20",
                             Patterns::List(Patterns::Double(0)),
                             "List of material parameter values. Units: Pa^-n*s^-1");
          prm.declare_entry ("Number of compositional fields with plastic rheology", "2",
                             Patterns::Double(),
                             "Number of compositional fields with plastic rheology, used in calculation of viscosity");
          prm.declare_entry ("Viscosity averaging scheme", "harmonic",
                             Patterns::Selection("arithmetic|harmonic|geometric|maximum composition"),
                             "When more than one compositional field is present at a point "
                             "with different viscosities, we need to come up with an average "
                             "viscosity at that point.  Select a weighted harmonic, arithmetic, "
                             "geometric, or maximum composition.");

          prm.declare_entry ("Minimum strain rate", "1e-25",
                             Patterns::Double (0),
                             "Minimum strain rate for viscosity stabilization"
                             "Unitless");
          prm.declare_entry ("Reference strain rate", "1e-15",
                             Patterns::Double (0),
                             "Reference strain rate for viscosity stabilization"
                             "Unitless");
          prm.declare_entry ("Reference compressibility", "4e-12",
                             Patterns::Double (0),
                             "The value of the reference compressibility. "
                             "Units: $1/Pa$.");							 

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    nz<dim>::parse_parameters (ParameterHandler &prm)
    {
      const unsigned int n_fields = this->n_compositional_fields();
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("nz model");
        {
          
          eta_min                    = prm.get_double ("Minimum viscosity");
          eta_max                    = prm.get_double ("Maximum viscosity");
          eta_reference              = prm.get_double ("Reference viscosity");
          composition_viscosity_prefactor = prm.get_double ("Composition viscosity prefactor");
          thermal_viscosity_exponent = prm.get_double ("Thermal viscosity exponent");
          t_conductivity             = prm.get_double ("Thermal conductivity");
		  t_diffusivity              = prm.get_double ("Thermal diffusivity");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          mantle_viscosity           = prm.get_double ("Mantle viscosity");
          num_plastic                = prm.get_double ("Number of compositional fields with plastic rheology");
          min_strain_rate            = prm.get_double ("Minimum strain rate");
          ref_strain_rate            = prm.get_double ("Reference strain rate");
		  reference_compressibility  = prm.get_double ("Reference compressibility");
          // if (thermal_viscosity_exponent!=0.0 && reference_T == 0.0)
            // AssertThrow(false, ExcMessage("Error: Material model simple with Thermal viscosity exponent can not have reference_T=0."));

          std::vector<double> x_values;

          //------FROM MD ---------------

          // Parse densities
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Densities")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of density list must be either one, or n_compositional_fields"));
          if (x_values.size() == 1)
            densities.assign( n_fields , x_values[0]);
          else
            densities = x_values;

          // Parse activation energies
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Activation energies")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of activation energy list must be either one, or n_compositional_fields"));
          if (x_values.size() == 1)
            activation_energies.assign( n_fields , x_values[0] );
          else
            activation_energies = x_values;

          // Parse activation volumes
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Activation volumes")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of activation volumes list must be either one, or n_compositional_fields"));
          if (x_values.size() == 1)
            activation_volumes.assign( n_fields , x_values[0] );
          else
            activation_volumes = x_values;

          // Parse thermal expansivities
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Thermal expansivities")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of thermal expansivity list must be either one, or n_compositional_fields"));
          if (x_values.size() == 1)
            thermal_alpha.assign( n_fields , x_values[0]);
          else
            thermal_alpha = x_values;
		
          // Parse stress exponents
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Stress exponents")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of nv list must be either one, or n_compositional_fields"));
          if (x_values.size() == 1)
            nvs.assign( n_fields , x_values[0]);
          else
            nvs = x_values;

          // Parse temeparture exponents
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Temperature exponents")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of Temperature exponents list must be either one, or n_compositional_fields"));
          if (x_values.size() == 1)
            nts.assign( n_fields , x_values[0]);
          else
            nts = x_values;					
					
		  // Parse plastic exponents
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Plastic exponents")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of np list must be either one, or n_compositional_fields"));
          if (x_values.size() == 1)
            nps.assign( n_fields , x_values[0]);
          else
            nps = x_values;
		
		  // Parse material exponents
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Material exponents")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of na list must be either one, or n_compositional_fields"));
          if (x_values.size() == 1)
            nas.assign( n_fields , x_values[0]);
          else
            nas = x_values;		

          // Parse material parameters
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Material parameters")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of Material parameters list must be either one, or n_compositional_fields"));
          if (x_values.size() == 1)
            material_parameters.assign( n_fields , x_values[0]);
          else
            material_parameters = x_values;


          //-----------END FROM MD ----------------

          // Parse angles of internal friction
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Angle of internal friction")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of Angle of internal friction list must be either one, or n_compositional_fields"));
          if (x_values.size() == 1)
            angle_if.assign( n_fields , x_values[0]);
          else
            angle_if = x_values;


          // Parse cohesions
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Cohesion")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of cohesion list must be either one, or n_compositional_fields"));
          if (x_values.size() == 1)
            cohesion.assign( n_fields , x_values[0]);
          else
            cohesion = x_values;

          //Parse averaging scheme
          if (prm.get ("Viscosity averaging scheme") == "harmonic")
            viscosity_averaging = harmonic;
          else if (prm.get ("Viscosity averaging scheme") == "arithmetic")
            viscosity_averaging = arithmetic;
          else if (prm.get ("Viscosity averaging scheme") == "geometric")
            viscosity_averaging = geometric;
          else if (prm.get ("Viscosity averaging scheme") == "maximum composition")
            viscosity_averaging = maximum_composition;
          else
            AssertThrow(false, ExcMessage("Not a valid viscosity averaging scheme"));
					
					//Parse reference temperature
					x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Reference temperature")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of Reference temperature list must be either one, or n_compositional_fields"));
          if (x_values.size() == 1)
            reference_T.assign( n_fields , x_values[0]);
          else
            reference_T = x_values;
					
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(nz,
                                   "nz",
                                   "A material model that has constant values "
                                   "for all coefficients but the density and viscosity. The defaults for all "
                                   "coefficients are chosen to be similar to what is believed to be correct "
                                   "for Earth's mantle. All of the values that define this model are read "
                                   "from a section ``Material model/Rigidplastic model'' in the input file, see "
                                   "Section~\\ref{parameters:Material_20model/Rigidplastic_20model}."
                                   "\n\n"
                                   "This model uses the following set of equations for the two coefficients that "
                                   "are non-constant: "
                                   "\\begin{align}"
                                   "  \\eta(p,T,\\mathfrak c) &= \\tau(T) \\zeta(\\mathfrak c) \\eta_0, \\\\"
                                   "  \\rho(p,T,\\mathfrak c) &= \\left(1-\\alpha (T-T_0)\\right)\\rho_0 + \\Delta\\rho \\; c_0,"
                                   "\\end{align}"
                                   "where $c_0$ is the first component of the compositional vector "
                                   "$\\mathfrak c$ if the model uses compositional fields, or zero otherwise. "
                                   "\n\n"
                                   "The temperature pre-factor for the viscosity formula above is "
                                   "defined as "
                                   "\\begin{align}"
                                   "  \\tau(T) &= H\\left(e^{\\beta (T-T_0)/T_0}\\right),"
                                   "  \\qquad\\qquad H(x) = \\begin{cases}"
                                   "                            10^{-2} & \\text{if}\\; x<10^{-2}, \\\\"
                                   "                            x & \\text{if}\\; 10^{-2}\\le x \\le 10^2, \\\\"
                                   "                            10^{2} & \\text{if}\\; x>10^{2}, \\\\"
                                   "                         \\end{cases}"
                                   "\\end{align} "
                                   "where $\\beta$ corresponds to the input parameter ``Thermal viscosity exponent'' "
                                   "and $T_0$ to the parameter ``Reference temperature''. If you set $T_0=0$ "
                                   "in the input file, the thermal pre-factor $\\tau(T)=1$."
                                   "\n\n"
                                   "The compositional pre-factor for the viscosity is defined as "
                                   "\\begin{align}"
                                   "  \\zeta(\\mathfrak c) &= \\xi^{c_0}"
                                   "\\end{align} "
                                   "if the model has compositional fields and equals one otherwise. $\\xi$ "
                                   "corresponds to the parameter ``Composition viscosity prefactor'' in the "
                                   "input file."
                                   "\n\n"
                                   "Finally, in the formula for the density, $\\Delta\\rho$ "
                                   "corresponds to the parameter ``Density differential for compositional field 1''."
                                   "\n\n"
                                   "Note that this model uses the formulation that assumes an incompressible "
                                   "medium despite the fact that the density follows the law "
                                   "$\\rho(T)=\\rho_0(1-\\beta(T-T_{\\text{ref}}))$. "
                                   "\n\n"
                                   "\\note{Despite its name, this material model is not exactly ``simple'', "
                                   "as indicated by the formulas above. While it was originally intended "
                                   "to be simple, it has over time acquired all sorts of temperature "
                                   "and compositional dependencies that weren't initially intended. "
                                   "Consequently, there is now a ``simpler'' material model that now fills "
                                   "the role the current model was originally intended to fill.}")
  }
}
