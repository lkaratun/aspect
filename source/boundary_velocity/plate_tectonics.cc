#include <aspect/boundary_velocity/plate_tectonics.h>
#include <math.h>

namespace aspect
{
  namespace BoundaryVelocity
  {
    template <typename T> int sgn(T val) {
      return (T(0) < val) - (val < T(0));
    }

    template <int dim>
    Tensor<1,dim>
    plate_tectonics<dim>::
    boundary_velocity (const Point<dim> &p) const
    {
      Tensor<1,dim> velocity;
      for (int i = 0; i<dim; ++i)
        velocity[i]=0;
      //coordinates
      const double x=p(0);
      const double y=p(1);
      const double z=p(dim-1);

      const double cm=0.01;
      double vel_x,vel_y,vel_z;
      if (dim==3)
      {
       vel_x=vel*cm*cos(angle*M_PI/180);
       vel_y=vel*cm*sin(angle*M_PI/180);
     }
     else
       vel_x=vel*cm;


    if (transition_type == "linear")
    {
     const double slm=total_thickness-lithospheric_thickness;

	  //influx through crust+lm
     double vel_x_in=vel_x;
      //outflux through sublith. mantle (excluding transition zone)
	  //double vel_x_out=-(vel_x_in*mantle_thickness-vel_x_in*sqrt(mantle_thickness*mantle_thickness-2*transition_zone*lithospheric_thickness+transition_zone*transition_zone))/transition_zone;


	  //Solving quadratic equation for vel_out
     double a = -transition_zone+2*slm;
     double b = 2*slm*vel_x_in-2*transition_zone*vel_x_in-2*lithospheric_thickness*vel_x_in;
     double c = -2*vel_x_in*vel_x_in*lithospheric_thickness-vel_x_in*vel_x_in*transition_zone;
     double D = b*b-4*a*c;
     double x1 = (-b + sqrt(D))/(2*a);
     double x2 = (-b - sqrt(D))/(2*a);
     double vel_x_out =0;


     if (x1*x2 > 0) std::cout<<std::endl<<"BOTH ROOTS FOR VELOCITY OUTFLOW EQUATION HAVE THE SAME SIGN"<<std::endl;
     else if (x1>0) vel_x_out = x1;
     else if (x2>0) vel_x_out = x2;
     else std::cout<<std::endl<<"BOTH ROOTS FOR VELOCITY OUTFLOW EQUATION ARE 0 OR NEGATIVE"<<std::endl;


     double transition_point = transition_zone * (vel_x_in/vel_x_out)/2;
     double transition_zone_mantle=abs(transition_zone*vel_x_out/vel_x_in);

     double vel_x_in_left = vel_x_in * 2*(1-influx_assymetry);
     double vel_x_in_right = vel_x_in * 2*influx_assymetry;
     double vel_x_out_left = vel_x_out * 2*(1-outflux_assymetry);
     double vel_x_out_right = vel_x_out * 2*outflux_assymetry;



	  //Proper b.c. rewritten
      //const double alfa = 45*M_PI/180; //Alpine fault angle

      //In- and out-flux for side faces
		if (x<width/2-wz_width/2)//left side
    {
      if (z > total_thickness-lithospheric_thickness)
       velocity[0]=vel_x_in_left;
     else if (z>total_thickness-lithospheric_thickness-transition_zone)
     {
       velocity[0]=(vel_x_in_left+vel_x_out_left)*(1-(total_thickness-lithospheric_thickness-z)/transition_zone) - vel_x_out_left;
     }
     else
       velocity[0]=-vel_x_out_left;


   }

      else if (x>width/2+wz_width/2) //right side-
      {
        if (z > total_thickness-lithospheric_thickness)
         velocity[0]=-vel_x_in_right;
       else if (z>total_thickness-lithospheric_thickness-transition_zone)
       {
         velocity[0]=-((vel_x_in_right+vel_x_out_right)*(1-(total_thickness-lithospheric_thickness-z)/transition_zone) - vel_x_out_right);
       }
       else
         velocity[0]=vel_x_out_right;




     }


     if (dim==3)
     {
       if (x<width/2-wz_width/2)
       {
        velocity[1]=vel_y;
        velocity[1]=-vel_y;
      }
      else
      {
        velocity[1]=-vel_y;
        velocity[1]=vel_y;
      }

    }
  }// linear transition

  else if (transition_type == "sin")
  {
    //Using depth starting from 1/4 to 3/4 as transition zone
    if (z > total_thickness*3/4)
      //Using sgn to make it work for both right and left walls in 1 line
      velocity[0] = vel_x*sgn(width/2 - x);
    else if (z < total_thickness*1/4)
      velocity[0] = -vel_x*sgn(width/2 - x);
    else
      velocity[0] = vel_x*sgn(width/2 - x) * std::sin((z - total_thickness/2) * M_PI / (total_thickness/2));















     if (dim==3)
     {
       if (x<width/2-wz_width/2)
       {
        velocity[1]=vel_y;
        velocity[1]=-vel_y;
      }
      else
      {
        velocity[1]=-vel_y;
        velocity[1]=vel_y;
      }



  }

	  //Debug
	   /*
	   if (x==0)
	   {
		//std::cout<< "vel_x_in_left="<<vel_x_in_left<<" vel_x_out_left="<<vel_x_out_left;


		std::ofstream myfile("vel", std::ios::app);
		myfile <<z <<" "<<velocity[0]<<std::endl;
		myfile.close();
		}
		*/


		// //Proper b.c. with a mistake
      // //const double alfa = 45*M_PI/180; //Alpine fault angle

      // //In- and out-flux for side faces
		// if (x<width/2-wz_width/2)//left side
        // {
          // if (z > total_thickness-lithospheric_thickness+transition_zone)
		  // {
			  // velocity[0]=vel_x_in_left;
			  // //velocity[0] *= 2*(1-influx_assymetry);
		  // }
          // else if (z>total_thickness-lithospheric_thickness-transition_zone_mantle)
		  // {
			  // velocity[0]=vel_x_out_left+(vel_x_in_left-vel_x_out_left)*(z-(total_thickness-lithospheric_thickness-transition_zone_mantle))/(transition_zone+transition_zone_mantle);
			  // //velocity[0] *= 2*(1-influx_assymetry);
		  // }
          // else
		  // {
			  // velocity[0]=vel_x_out_left;
			  // //velocity[0] *= 2*(1-outflux_assymetry);
		  // }


        // }

      // else if (x>width/2+wz_width/2) //right side-
        // {
          // if (z > total_thickness-lithospheric_thickness+transition_zone)
		  // {
			  // velocity[0]=-vel_x_in_right;
			  // //velocity[0] *= 2*influx_assymetry;
		  // }
          // else if (z>total_thickness-lithospheric_thickness-transition_zone_mantle)
		  // {
			  // velocity[0]=-(vel_x_out_right+(vel_x_in_right-vel_x_out_right)*(z-(total_thickness-lithospheric_thickness-transition_zone_mantle))/(transition_zone+transition_zone_mantle));
			  // //velocity[0] *= 2*influx_assymetry;
		  // }
          // else
		  // {
            // velocity[0]=-vel_x_out_right;
			// //velocity[0] *= 2*outflux_assymetry;
		  // }




		// }


			// if (dim==3)
				// {
					// if (x<width/2-wz_width/2)
						// {
						// velocity[1]=vel_y;
						// velocity[1]=-vel_y;
						// }
					// else
					// {
						// velocity[1]=-vel_y;
						// velocity[1]=vel_y;
					// }

				// }





      /*
      if (abs(velocity[1])<=1e-3)
      {}
      else
        std::cout<<"y velocity: "<<velocity[1]<<"  ";
      */


      // TEST
      /*
      velocity[0]=0;
      if (x<width/2-wz_width/2)//left side
      {
          // if (z > total_thickness-lithospheric_thickness+transition_zone)
              // velocity[0]=1;
          // else if (z>total_thickness-lithospheric_thickness-transition_zone_mantle)
              // velocity[0]=1;
          // else
            if (z > total_thickness/2)
              velocity[0]=1;
            else
              velocity[0]=-1;
            velocity[0]=vel_x_in*(z-depth/2)/(depth/2);
      }

      else if (x>width/2+wz_width/2) //right side-
      {
          // if (z > total_thickness-lithospheric_thickness+transition_zone)
              // velocity[0]=-1;
          // else if (z>total_thickness-lithospheric_thickness-transition_zone_mantle)
              // velocity[0]=-1;
          // else
            if (z > total_thickness/2)
              velocity[0]=-1;
            else
              velocity[0]=1;
            velocity[0]=-vel_x_in*(z-depth/2)/(depth/2);
      }

          velocity[1]=0;
          velocity[2]=0;
      */
      // END TEST

      // strike-slip for side faces
      // if (x==0)
      // velocity[1]=vel_y;
      // if (x==width)
      // velocity[1]=-vel_y;


      /*
      //Front and back faces: inclined fault
      if (x<width/2-wz_depth/tan(alfa)-wz_width/2) //left part
          velocity[1]=vel_y;
      else if (x>width/2+wz_depth/tan(alfa)+wz_width/2) //right part
          velocity[1]=-vel_y;
      else if (z<depth-wz_depth+(width/2-x)*tan(alfa)-wz_width*tan(alfa)/2) //middle part
          velocity[1]=vel_y;
      else if (z>depth-wz_depth+(width/2-x)*tan(alfa)+wz_width*tan(alfa)/2) //middle part
          velocity[1]=-vel_y;
      */



      //Front and back faces: vertical fault
      /*
      if (z>=total_thickness-lithospheric_thickness) //apply y-velocity to lithosphere only
      {
          if (x<width/2-wz_width*2) //left part
              velocity[1]=vel_y;
          else if (x>width/2+wz_width*2) //right part
              velocity[1]=-vel_y;
      }
      */

      /*
      //B.c. v.2.0: outflux through the bottom
      if (x==0)
      {
        velocity[0]=vel_x_in;
        velocity[1]=vel_y;
      }
      if (x==width)
      {
        velocity[0]=-vel_x_in;
        velocity[1]=-vel_y;
      }
      if (dim == 3)
      {
        velocity[2]=0;
        if (z==0) velocity[2] = -vel_x_in*2*depth/width;
      }
      if (dim ==2)
      {
        velocity[1]=0;
        if (y==0) velocity[1] = -vel_x_in*2*depth/width;
      }
      */

      /*
      //B.c. v.2.1: linearly varying outflux through the bottom
      if (x==0)
      {
        velocity[0]=vel_x_in*z/depth;
        velocity[1]=vel_y;
      }
      if (x==width)
      {
        velocity[0]=-vel_x_in*z/depth;
        velocity[1]=-vel_y;
      }
      if (dim == 3)
      {
        //velocity[2]=0;
        if (z==0) velocity[2] = -vel_x_in*2*depth/width*    (-abs(x-width/2)+width/2)/(width/2);
      }

      */


      // //B.c. v.1.5: in/outflux through sidewalls, constant gradient
      // if (x==0)
      // {
        // velocity[0]=vel_x_in*(z-depth/2)/(depth/2);
        // if (dim==3) velocity[1]=vel_y;
      // }
      // if (x==width)
      // {
        // velocity[0]=-vel_x_in*(z-depth/2)/(depth/2);
        // if (dim==3) velocity[1]=-vel_y;
      // }
      // if (dim == 3)
      // {
        // velocity[2]=0;
        // //if (z==0) velocity[2] = -vel_x_in*2*depth/width*    (-abs(x-width/2)+width/2)/(width/2);
      // }



    if (this->convert_output_to_years())
      return velocity / year_in_seconds;
    else
      return velocity;
  }
}




    template <int dim>
  void
  plate_tectonics<dim>::declare_parameters (ParameterHandler &prm)
  {
    prm.enter_subsection("Boundary velocity model");
    {
      prm.enter_subsection("Plate tectonics");
      {
        prm.declare_entry ("Obliquity angle", "0",
         Patterns::Double(0),
         "Angle of velocity obliquity. 0 is pure convergence, 90 is pure strike-slip motion."
         "Units: degrees");

        prm.declare_entry ("Velocity", "1",
         Patterns::Double (0),
         "Total velocity value. Units: $cm/yr$ or $cm/s$.");
        prm.declare_entry ("Weak zone width", "40e3",
         Patterns::Double (0),
         "Weak zone width. Units: $m$.");
        prm.declare_entry ("Weak zone height", "120e3",
         Patterns::Double (0),
         "Weak zone height. Units: $m$.");
        prm.declare_entry ("Weak zone depth", "60e3",
         Patterns::Double (0),
         "Weak zone depth. Units: $m$.");
        prm.declare_entry ("Total width", "1024e3",
         Patterns::Double (0),
         "Total model width. Units: $m$.");
        prm.declare_entry ("Total height", "512e3",
         Patterns::Double (0),
         "Total model height. Units: $m$.");
        prm.declare_entry ("Lithospheric thickness", "120e3",
         Patterns::Double (0),
         "Lithospheric thickness. Units: $m$.");
        prm.declare_entry ("Transition zone height", "150e3",
         Patterns::Double (0),
         "Transition zone for velocity on side faces (left and right). Units: $m$.");
        prm.declare_entry ("Weak zone angle", "45",
         Patterns::Double (0),
         "Angle of the weak zone (Alpine fault). 0 is horizontal, 90 is vertical Units: $deg$.");
        prm.declare_entry ("Influx assymetry", "0.5",
         Patterns::Double (0),
         "0 is influx from the left side, 0.5 is symmetric, 1 is from the right side");
        prm.declare_entry ("Outflux assymetry", "0.5",
         Patterns::Double (0),
         "0 is outflux from the left side, 0.5 is symmetric, 1 is from the right side");
        prm.declare_entry ("Transition type", "linear",
         Patterns::Selection("linear|sin"),
         "Shape of the velocity within the transition zone from outflux to influx");

      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }



    template <int dim>
  void
  plate_tectonics<dim>::parse_parameters (ParameterHandler &prm)
  {
    prm.enter_subsection("Boundary velocity model");
    {
      prm.enter_subsection("Plate tectonics");
      {
        vel                         = prm.get_double ("Velocity");
        angle                       = prm.get_double ("Obliquity angle");
        alfa                        = prm.get_double ("Weak zone angle");
        wz_width                    = prm.get_double ("Weak zone width");
        wz_height                   = prm.get_double ("Weak zone height");
        wz_depth                    = prm.get_double ("Weak zone depth");
        width                       = prm.get_double ("Total width");
        depth                       = prm.get_double ("Total height");
        total_thickness             = prm.get_double ("Total height");
          lithospheric_thickness      = prm.get_double ("Lithospheric thickness");//+prm.get_double ("Transition zone height");
          transition_zone             = prm.get_double ("Transition zone height");
          influx_assymetry	          = prm.get_double ("Influx assymetry");
          outflux_assymetry	          = prm.get_double ("Outflux assymetry");
          transition_type             = prm.get("Transition type");
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
  namespace BoundaryVelocity
  {
    ASPECT_REGISTER_BOUNDARY_VELOCITY_MODEL(plate_tectonics,
     "plate tectonics",
     "My b.c.")
  }
}

