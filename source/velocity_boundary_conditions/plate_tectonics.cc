#include <aspect/velocity_boundary_conditions/plate_tectonics.h>
#include <math.h>

namespace aspect
{
  namespace VelocityBoundaryConditions
  {
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
      double vel_x=vel*cm*cos(angle*M_PI/180);
      double vel_y=vel*cm*sin(angle*M_PI/180);
      double vel_z=0;

      //total_depth=512e3, total_width=1536e3, crustal_depth=20e3, lithospheric_depth=120e3, wz_width=40.e3, wz_height=120.e3, wz_depth=60e3
//        const double wz_width=40.e3;
//        const double wz_height=120.e3;
//        const double wz_depth=60.e3;
//        const double width=1024.e3;
//        const double depth=512.e3;
//        const double total_thickness=512e3;
//        double lithospheric_thickness=120e3;
      const double mantle_thickness=total_thickness-lithospheric_thickness;
//        const double transition_zone=50e3;
      //lithospheric_thickness+=transition_zone; //move the upper boundary of the transition zone to the base of lithosphere

      double vel_x_in=vel_x;
      double vel_x_out=-(vel_x_in*mantle_thickness-vel_x_in*sqrt(mantle_thickness*mantle_thickness-2*transition_zone*lithospheric_thickness+transition_zone*transition_zone))/transition_zone;
      //std::cout<<"vel_x_out: "<<vel_x_out<<"  ";
      double transition_zone_mantle=abs(transition_zone*vel_x_out/vel_x_in);





      //const double alfa = 45*M_PI/180; //Alpine fault angle

      //In- and out-flux for side faces
				if (x<width/2-wz_width/2)//left side
        {
          if (z > total_thickness-lithospheric_thickness+transition_zone)
            velocity[0]=vel_x_in;
          else if (z>total_thickness-lithospheric_thickness-transition_zone_mantle)
            velocity[0]=vel_x_out+(vel_x_in-vel_x_out)*(z-(total_thickness-lithospheric_thickness-transition_zone_mantle))/(transition_zone+transition_zone_mantle);
          else
            velocity[0]=vel_x_out;
        }

      else if (x>width/2+wz_width/2) //right side-
        {
          if (z > total_thickness-lithospheric_thickness+transition_zone)
            velocity[0]=-vel_x_in;
          else if (z>total_thickness-lithospheric_thickness-transition_zone_mantle)
            velocity[0]=-(vel_x_out+(vel_x_in-vel_x_out)*(z-(total_thickness-lithospheric_thickness-transition_zone_mantle))/(transition_zone+transition_zone_mantle));
          else
            velocity[0]=-vel_x_out;
        }



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

      /*
      //B.c. v.1.5: in/outflux through sidewalls, constant gradient
      if (x==0)
      {
        velocity[0]=vel_x_in*(z-depth/2)/(depth/2);
        velocity[1]=vel_y;
      }
      if (x==width)
      {
        velocity[0]=-vel_x_in*(z-depth/2)/(depth/2);
        velocity[1]=-vel_y;
      }
      if (dim == 3)
      {
        velocity[2]=0;
        //if (z==0) velocity[2] = -vel_x_in*2*depth/width*    (-abs(x-width/2)+width/2)/(width/2);
      }
      */


      if (this->convert_output_to_years())
        return velocity / year_in_seconds;
      else
        return velocity;
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

          prm.declare_entry ("Velocity", "2",
                             Patterns::Double (0),
                             "Total velocity value. Units: $m/s$.");
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
          prm.declare_entry ("Transition zone height", "50e3",
                             Patterns::Double (0),
                             "Transition zone for velocity on side faces (left and right). Units: $m$.");
          prm.declare_entry ("Weak zone angle", "45",
                             Patterns::Double (0),
                             "Angle of the weak zone (Alpine fault). 0 is horizontal, 90 is vertical Units: $deg$.");

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
          lithospheric_thickness      = prm.get_double ("Lithospheric thickness")+prm.get_double ("Transition zone height");
          transition_zone             = prm.get_double ("Transition zone height");
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
  namespace VelocityBoundaryConditions
  {
    ASPECT_REGISTER_VELOCITY_BOUNDARY_CONDITIONS(plate_tectonics,
                                                 "plate tectonics",
                                                 "My b.c.")
  }
}

