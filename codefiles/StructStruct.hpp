#ifndef _STRUCTURE_STRUCT__H
#define _STRUCTURE_STRUCT__H

#include <iostream>
#include <Eigen/Dense>

struct P_P{
    /**
     * @brief Structure for holding the various properties of the solution. 'y' is an array containing the solutions to the dependent variables, where there are 'I' variables and 'J' cells. 'm' is the indendent variable. 'ym' and 'mm' are the values halfway between the 'm' points
     */

    Eigen::ArrayXd m , mm , dm , dmm, rho, scale , Cp , del_rho , nabla , T0 , P0 , g , H;
    Eigen::ArrayXd diffrhoP, diffrhoT , diffdel_rhoP, diffdel_rhoT , diffCpP , diffCpT , Bavg , GM_R4 , diffGM_R4 ;
    Eigen::ArrayXXd y , dy , ym ;
    Eigen::ArrayXXd diffnabla;
    double dt, M, n , mu_m , Tsurf[7] , Psurf[7];
    char Tbound , Pbound , nabtype;
    // Add more useful arrays - derivative properties from y,m, but be sure to also add to 'resize' and the consructor below!

    void resize(int J) //Resizes the various arrays if one wishes to change the number of cells
    {
        int I = scale.size();
        m.resize(J) ; mm.resize(J) ; dm.resize(J-1) ; dmm.resize(J-1); rho.resize(J); Cp.resize(J) ; del_rho.resize(J) ; nabla.resize(J); T0.resize(J) ; P0.resize(J) ; H.resize(J-1) ; GM_R4.resize(J-1) ; diffGM_R4.resize(J-1) ; g.resize(J) ; y.resize(J,I) ; dy.resize(J-1,I) ; Bavg.resize(J-1) ; 
        diffnabla.resize(J,I) ; diffrhoP.resize(J); diffrhoT.resize(J) ; diffdel_rhoP.resize(J); diffdel_rhoT.resize(J); diffCpP.resize(J); diffCpT.resize(J) ; 
    }
    P_P(int I, int J):
        m(J) , mm(J) , dm(J-1) , dmm(J-1), rho(J), scale(I), Cp(J) , del_rho(J) , nabla(J), y(J,I) , dy(J-1,I) , ym(J,I) , Bavg(J-1) , GM_R4(J-1) , diffGM_R4(J-1) , H(J-1) , diffnabla(J,I) , diffrhoP(J), diffrhoT(J) , diffdel_rhoP(J), diffdel_rhoT(J), diffCpP(J), diffCpT(J) 
        {}
};

#endif