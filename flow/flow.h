#ifndef _FLOW_H
#define _FLOW_H

namespace flow {
  
  class StateVector {
    private:
      std::vector<double> state(3);            // Conservative formulation
    
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

      void setPrimitives(const double& density, const double& u_velocity, const double& pressure);
      void setConservatives(const double& density, const double& u_momentum, const double& total_energy);
  };

  class GammaConstants {
    private:
      double gamma, g1, g2, g3, g4, g5, g6, g7, g8;
      GammaConstants(const GammaConstants& rhs);
      GammaConstants& operator=(const GammaConstants& rhs);
      
    public:
      explicit GammaConstants(const double& _gamma);
      ~GammaConstants();
      
      const double getGamma() const;
      const double getG1() const;
      const double getG2() const;
      const double getG3() const;
      const double getG4() const;
      const double getG5() const;
      const double getG6() const;
      const double getG7() const;
      const double getG8() const;
  };
  
}

#endif 
