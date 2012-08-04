#ifndef _FLOW_H
#define _FLOW_H

#include <vector>

namespace flow {

  namespace gammaconst {
    
    const double gamma = 1.4;
    const double g1 = (gamma - 1.0)/(2.0*gamma);
    const double g2 = (gamma + 1.0)/(2.0*gamma);
    const double g3 = 2.0*gamma/(gamma - 1.0);
    const double g4 = 2.0/(gamma - 1.0);
    const double g5 = 2.0/(gamma + 1.0);
    const double g6 = (gamma - 1.0)/(gamma + 1.0);
    const double g7 = (gamma - 1.0)/2.0;
    const double g8 = gamma - 1.0;

  }
  
  class StateVector {
    private:
      std::vector<double> state;            // Conservative formulation
    
    public:
      StateVector(const double& density,
		  const double& u_velocity,
		  const double& pressure);
      virtual ~StateVector();
    
      const double getDensity() const;
      const double getUVelocity() const;
      const double getUMomentum() const;
      const double getPressure() const;
      const double getTotalEnergy() const;

      void setPrimitives(const double& density, 
			 const double& u_velocity, 
			 const double& pressure);
      void setConservatives(const double& density, 
			    const double& u_momentum, 
			    const double& total_energy);
  };
  
}

#endif 
