#ifndef _FLOW_H
#define _FLOW_H

namespace flow {
  
  class StateVector {
    private:
      std::vector<double> state(3);
    
    public:
      explicit StateVector(const double& density,
	   		   const double& velocity,
			   const double& pressure);
      virtual ~StateVector();
    
      const double getDensity() const;
      const double getVelocity() const;
      const double getPressure() const; 
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
