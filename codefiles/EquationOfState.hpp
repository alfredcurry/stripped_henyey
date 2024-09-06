#ifndef _EoS_CALC_H
#define _EoS_CALC_H

#include <vector>
#include <Eigen/Dense>
#include <string>
#include "StructStruct.hpp"
#include <boost/math/interpolators/pchip.hpp>

/// Class for evaluation the equation of state. EoS type determined at compile time, with certain parameter tweaks available at run-time
class EoS{
    public:
        template <class Real> /// Computes the values of the EoS at each point in the planet
        void compute_values(struct P_P *planet );
        template <class Real>  /// Computes the derivatives the EoS values at each point in the planet
        void compute_derivs(struct P_P *planet );
        template <class Real>  /// Computes the values of the density only
        void compute_density_only( struct P_P *planet );
        EoS(double *_variables);
        std::string EoStype;

    private:
        int J;
        std::vector<double> variables;
};


#endif