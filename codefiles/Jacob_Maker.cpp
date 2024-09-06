#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "StructStruct.hpp"
#include "Jacob_Maker.hpp"
#include "PhysicalConsts.h"
#include "adiabatnab.hpp"
#include <stdexcept>

using namespace Eigen;

double JacobMaker::Gvec(struct P_P *planet , EoS *eos , double *Gv ){
            Map<Matrix<double , Dynamic , Dynamic , RowMajor> > G(Gv, J , I);

            preprest(planet, eos);

            //Insert your own discretised equations here!
            G.block(0,0,J-1,1) = planet->dy.col(0) + planet->dmm*planet->GM_R4;            

            G.block(0,1,J-1,1) = planet->dy.col(1) - planet->dm * 3/(4*M_PI*planet->rho.head(J-1) * pow(planet->scale(1),3)/planet->M);   

            G.block(0,2,J-1,1) = planet->dy.block(0,2,J-1,1) + planet->dmm.head(J-1) * planet->GM_R4.head(J-1) * planet->Bavg.head(J-1) ;

            G.block(0,3,J-1,1) = planet->dy.col(3) -  planet->M /planet->scale(3) * planet->dm * ( planet->H - planet->Cp.head(J-1) * (planet->scale(2) * planet->y.block(0,2,J-1,1)-planet->T0.head(J-1))/planet->dt + planet->del_rho.head(J-1) * (planet->scale(0)*planet->y.block(0,0,J-1,1)-planet->P0.head(J-1))/planet->dt ) ;          
                
            G(J-1 , 0) = planet->y(J-1,0) - (planet->Psurf[0])/planet->scale(0);
            G(J-1 , 1) = planet->y(J-1,1) - planet->y(J-2,1);
            G(J-1 , 2) = planet->y(J-1,2) - (planet->Tsurf[0])/planet->scale(2); 
            G(J-1 , 3) = planet->y(J-1,3) - planet->y(J-2,3);
                     
            return G.squaredNorm();
        }
       
void JacobMaker::GJacob(struct P_P *planet , EoS *eos , double *l , double *d , double *u ){
            int step = I*I;
                    std::cout << "before derivs" << std::endl;
            
            prepderivs(planet , eos);


            //inner boundary

            mat_t Gl0(l , I , I);
            mat_t Gd0(d , I , I);
            mat_t Gu0(u , I , I);
            
            Gl0 = MatrixXd::Zero(I,I);

            //Insert your own inner derivatives here
            Gd0 = Gdiff(planet->y.row(0), planet->y.row(1) , planet->scale , planet->M , planet->P0(0), planet->T0(0) , planet->dm(0) , planet->dmm(0), planet->rho(0) , planet->dt , planet->del_rho(0) , planet->Cp(0) ,  planet->nabla(0) , planet->GM_R4(0) , planet->diffGM_R4(0) , planet->Bavg(0) , planet->diffnabla.block(0,0,2,I), planet->diffCpP(0) , planet->diffCpT(0) , planet->diffrhoP(0) , planet->diffrhoT(0) , planet->diffdel_rhoP(0) , planet->diffdel_rhoT(0) );
           
            Gu0 = MatrixXd::Zero(I,I);
            Gu0(0,0) = 1 ;
            
            Gu0(2,0) = planet->dmm(0) * planet->GM_R4(0) * 0.5 * (-planet->y(1,2)*planet->nabla(1)/pow(planet->y(1,0),2) + planet->y(1,2)/planet->y(1,0) * planet->diffnabla(1,0));
            Gu0(2,1) = planet->dmm(0) * planet->GM_R4(0) * 0.5 * planet->y(0,2)/planet->y(0,0) * 0.5 * planet->diffnabla(0+1,1);
            Gu0(2,2) = 1 + planet->dmm(0) * planet->GM_R4(0) * 0.5 * (planet->nabla(1)+planet->y(1,2)*planet->diffnabla(1,2) )/planet->y(1,0) ; 
            Gu0(2,3) = planet->dmm(0) * planet->GM_R4(0) * 0.5 * planet->y(1,2)/planet->y(1,0) * 0.5 * planet->diffnabla(0+1,3);
                    std::cout << "afer more derivs" << std::endl;

            //non boundary diff_Gs            
            for(int j=1 ; j < J-1 ; j++){
                //leaving empty step at start for l
                mat_t Gl(l + j*step , I , I);
                mat_t Gd(d + j*step, I , I);
                mat_t Gu(u + j*step , I , I);
            
                //Insert your own derivatives here
                Gd = Gdiff(planet->y.row(j), planet->y.row(j+1) , planet->scale , planet->M , planet->P0(j), planet->T0(j) , planet->dm(j) , planet->dmm(j) , planet->rho(j) , planet->dt , planet->del_rho(j) , planet->Cp(j) ,  planet->nabla(j) , planet->GM_R4(j) , planet->diffGM_R4(j) , planet->Bavg(j) , planet->diffnabla.block(j,0,2,I), planet->diffCpP(j) , planet->diffCpT(j) , planet->diffrhoP(j), planet->diffrhoT(j), planet->diffdel_rhoP(j) , planet->diffdel_rhoT(j) );
                
                Gu = MatrixXd::Zero(I,I);
                Gl = MatrixXd::Zero(I,I);

                Gu(0,0) = 1 ;
                Gl(1,1) = -3*pow(planet->y(j-1,1) , 2);
                Gl(3,3) = -1;
               
                Gu(2,0) = planet->dmm(j) * planet->GM_R4(j) * 0.5 * (-planet->y(j+1,2)*planet->nabla(j+1)/pow(planet->y(j+1,0),2) + planet->y(j+1,2)/planet->y(j+1,0) * planet->diffnabla(j+1,0));
                Gu(2,1) = planet->dmm(j) * planet->GM_R4(j) * 0.5 * planet->y(j,2)/planet->y(j,0) * 0.5*planet->diffnabla(j+1,1);
                Gu(2,2) = 1 + planet->dmm(j) * planet->GM_R4(j) * 0.5 * (planet->nabla(j+1)+planet->y(j+1,2)*planet->diffnabla(j+1,2)  )/planet->y(j+1,0) ; 
                Gu(2,3) = planet->dmm(j) * planet->GM_R4(j)  * 0.5 * planet->y(j,2)/planet->y(j,0) * 0.5*planet->diffnabla(j+1,3);
                
                Gl(2,1) = planet->dmm(j) * planet->GM_R4(j) * 0.5 * planet->y(j,2)/planet->y(j,0) * 0.5*planet->diffnabla(j,1);
                Gl(2,3) = planet->dmm(j) * planet->GM_R4(j) * 0.5 * planet->y(j,2)/planet->y(j,0) * 0.5*planet->diffnabla(j,3);
                
              
            }

            //outer boundary - note will obviously change as outer BC change to more realistic
            mat_t GlJ(l + (J-1)*step , I , I);
            mat_t GdJ(d + (J-1)*step , I , I);
            mat_t GuJ(u + (J-1)*step , I , I);
            
            GdJ = MatrixXd::Zero(I,I);
            GlJ = MatrixXd::Zero(I,I);
            GuJ = MatrixXd::Zero(I,I);

            //Insert your own outer derivatives here
            GdJ(0,0) = 1 ;
            GdJ(0,1) = - ( planet->Psurf[1+1*1]) *planet->scale(1)/planet->scale(0);
            GdJ(0,2) = - (+ planet->Psurf[1+1*2]) *planet->scale(2)/planet->scale(0);

            GdJ(1,1) = 1;
            
            GdJ(2,0) = 0;
            GdJ(2,1) = - ( planet->Tsurf[1+1*1])*planet->scale(1)/planet->scale(2); 
            GdJ(2,2) = 1  ;
            GdJ(2,3) = - (  planet->Tsurf[1+1*3])*planet->scale(3)/planet->scale(2); 

            GdJ(3,3) = 1;

            GlJ(0,1) = - 0;
            GlJ(0,3) = - 0;

            GlJ(2,1) = -0; 
            GlJ(2,3) = - 0; 

            GlJ(1,1) = -1;
            GlJ(3,3) = -1;
                        
        }

JacobMaker::JacobMaker(int _I, int _J)
        :
            I(_I), J(_J) 
        {};
    
MatrixXd JacobMaker::Gdiff( const ArrayXd y  , const ArrayXd yplus ,  const ArrayXd scale , double M , double P0,  double T0, const double dm, const double dmm, const double rho , const double dt, const double del_rho , const double Cp , const double nabla , double GM_R4 , double diffGM_R4 ,  double Bavg , ArrayXXd diffnabla , double diffCpP , double diffCpT , double diffrhoP , double diffrhoT , double diffdel_rhoP , double diffdel_rhoT ){
            Matrix<double, Dynamic, Dynamic,Eigen::RowMajor> G0;
            G0 = MatrixXd::Zero(I,I);

            G0(0,0) = -1.0;
            G0(0,1) = dmm * diffGM_R4;

            G0(1,0) = dm * 3/(4*M_PI*pow(rho,2)* pow(scale(1),3)/M) * diffrhoP*scale(0);
            G0(1,1) = 3 * pow(y(1) , 2);
            G0(1,2) = dm * 3/(4*M_PI*rho*rho * pow(scale(1),3)/M) * diffrhoT*scale(2);
           
            G0(2,0) = dmm * GM_R4 * 0.5 * (-y(2)*nabla/pow(y(0),2) + y(2)/y(0) * diffnabla(0,0));
            G0(2,1) = dmm * ( diffGM_R4 * Bavg + GM_R4* 0.5 * (y(2)/y(0) * 0.5*diffnabla(0,1) + yplus(2)/yplus(0)*0.5*diffnabla(1,1))); 
            G0(2,2) = -1 + dmm * GM_R4 * 0.5 * (nabla/y(0) + y(2)/y(0) * diffnabla(0,2));
            G0(2,3) = dmm * GM_R4 * 0.5 * ( y(2)/y(0)*0.5*diffnabla(0,3) + yplus(2)/yplus(0)*0.5*diffnabla(1,3));

            G0(3,0) = -M*dm/dt / scale(3) * scale(0) * ( - diffCpP*(scale(2)*y(2)-T0) +  del_rho + diffdel_rhoP * (scale(0)*y(0) - P0) )  ; 
            G0(3,2) = -M*dm/dt / scale(3) * scale(2) * ( -Cp - diffCpT*(scale(2)*y(2)-T0) + diffdel_rhoT * (scale(0)*y(0)-P0)  ) ;
            G0(3,3) = 1;
            return G0;
        }
        
void JacobMaker::preprest( struct P_P *planet , EoS *eos ){
            
        //differences
            planet->dy.col(0) = diff(planet->y.col(0));
            planet->dy.col(1) = diffminus(pow(planet->y.block(0,1,J-1,1), 3));
            planet->dy(0,1) = planet->dy(0,1) ;
            planet->dy.col(2) = diff(planet->y.col(2));
            planet->dy.col(3) = diffminus(planet->y.block(0,3,J-1,1));  
            planet->dy(0,3) = planet->dy(0,3) ; 
            //dimensional versions defined at cell centres (i.e. where P,T,rho are) - used for nabla calucluations
            planet->ym.col(0) = planet->y.col(0) * planet->scale(0);
            planet->ym.block(0,1,J-1,1) = (planet->y.block(0,1,J-1,1) - diffminus(planet->y.block(0,1,J-1,1))/2.0 )*planet->scale(1);
            planet->ym(J-1,1) = planet->y(J-1,1)*planet->scale(1);
           
            planet->ym.col(2) = planet->y.col(2) * planet->scale(2);
            planet->ym.block(0,3,J-1,1) = (planet->y.block(0,3,J-1,1) - planet->dy.col(3)/2.0 ) * planet->scale(3);
            planet->ym(J-1,3) = planet->y(J-1,3)*planet->scale(3);
        
            planet->g = G_Newt * planet->M * planet->mm /(planet->ym.col(1)*planet->ym.col(1));
        
            eos->compute_values<ArrayXd>( planet);
            std::cout << "nabtype " << planet->nabtype << std::endl;
            if (planet->nabtype == 'i'){
                planet->nabla = ArrayXd::Zero(planet->nabla.size());    
            }else if (planet->nabtype == 'a'){
                planet->nabla = Nab_Ad<ArrayXd>::nab(planet->ym.col(0) , planet->ym.col(2) , planet->del_rho , planet->Cp );
            }else{
                throw std::invalid_argument("Invalid nabtype. Valid types: 'i' or 'a'");
            }
           
            for( int j=0; j<J-1; j++){ //auxilliary parameter B average because of cell locations
                
                planet->Bavg(j) = 0.5*( planet->y(j,2)*planet->nabla(j)/planet->y(j,0) + planet->y(j+1,2)*planet->nabla(j+1)/planet->y(j+1,0) );
            }
            
            double G_Ns = pow(planet->M,2)*G_Newt/(planet->scale(0)*pow(planet->scale(1),4)); //Scaled G for P & T eqns

            planet->GM_R4.head(J-2) = planet->m.segment(1,J-2) * G_Ns /(4*M_PI*pow(planet->y.block(0,1,J-2,1),4));
            planet->GM_R4(J-2) = (planet->m(J-2) + planet->m(J-1))/2 * G_Ns /(4*M_PI*pow((planet->y(J-2,1)+planet->y(J-3,1))/2,4));
        }
        
void JacobMaker::prepderivs(struct P_P *planet , EoS *eos  ){

            //compute derivatives of various
            if( planet->Pbound == 'f' ){
                planet->Psurf[1+1*1] = 0;
		        planet->Psurf[1+1*2] = 0;
            }
            if( planet->Tbound == 'f' ){
                planet->Tsurf[1+1*1] = 0;
               planet->Tsurf[1+1*3] = 0;
            }
            planet->diffGM_R4.head(J-2) = - 4/planet->y.block(0,1,J-2,1) * planet->GM_R4.head(J-2);
            planet->diffGM_R4(J-2) = - 2*2/(planet->y(J-2,1)+planet->y(J-3,1))* planet->GM_R4(J-2);
            eos->compute_derivs<ArrayXd>( planet) ;

            if (planet->nabtype == 'i'){
                planet->diffnabla =0;
            }else
            if(planet->nabtype == 'a'){
                planet->diffnabla.col(0) = Nab_Ad<ArrayXd>::DiffnabP(planet->ym.col(0) , planet->ym.col(2) , planet->del_rho , planet->Cp ,  planet->diffdel_rhoP, planet->diffCpP) ;
                planet->diffnabla.col(1) = 0;
                planet->diffnabla.col(2) = Nab_Ad<ArrayXd>::DiffnabT(planet->ym.col(0) , planet->ym.col(2) , planet->del_rho , planet->Cp ,  planet->diffdel_rhoT, planet->diffCpT);
                planet->diffnabla.col(3) = 0;
            }

        }
        
ArrayXd JacobMaker::diffminus(const ArrayXd x){
            Eigen::ArrayXd dx(x.size());
            dx(0) = x(0);
            for(int j = 1; j < x.size() ; j++){
                dx(j) = x(j) - x(j-1);
            }
            return dx;
        }
ArrayXd JacobMaker::diffplus(const ArrayXd x){
            Eigen::ArrayXd dx(x.size());
            dx(x.size()-1) = -x(x.size()-1);
            for(int j = 0; j < x.size()-1 ; j++){
                dx(j) = x(j+1) - x(j);
            }
            return dx;
        }
ArrayXd JacobMaker::diff(const ArrayXd x){
            Eigen::ArrayXd dx(x.size()-1);
            for(int j = 0; j < x.size()-1 ; j++){
                dx(j) = x(j+1) - x(j);
            }
            return dx;
        }

ArrayXd JacobMaker::massbins(int J){
            ArrayXd m(J);
            int J_bound = J/10;
            m.head(J-2*J_bound) = ArrayXd::LinSpaced(J-2*J_bound,0,0.985);
            m.segment(J-2*J_bound, J_bound) = ArrayXd::LinSpaced(J_bound,0.985+0.02/J_bound,0.998);
            m.tail(J_bound) = ArrayXd::LinSpaced(J_bound,0.998+0.002/J_bound, 1);
            return m;
        }
void JacobMaker::prepm(struct P_P *planet ){
            //input is the 'Mj', find the 'Mj-1/2' and dMs
            int J = planet->m.size();
            planet->dm = JacobMaker::diff(planet->m);
            planet->mm.head(J-1) = planet->m.tail(J-1) - planet->dm/2 ;
            planet->mm(J-1) = planet->m(J-1);
            //if(planet->M_core > 0){
		//planet->mm(0) = planet->M_core/planet->M;
		//}
		planet->dmm = JacobMaker::diff(planet->mm);
        }
ArrayXd JacobMaker::power_law_m(int J , double X , double power){
    ArrayXd m = ArrayXd::LinSpaced(J , 0 , 1);
    m = pow(m , power); 
    m = m*X + (1-X);
    return m;
}
ArrayXd JacobMaker::power_law_Plus(int J , double X , double power , double mid){
    ArrayXd m(J) ;
    int Jfrac = J/2;
    m.head(J-Jfrac) = ArrayXd::LinSpaced(J-Jfrac ,1-X, 1-X + X*mid );
    m.tail(Jfrac + 1) = ArrayXd::LinSpaced(Jfrac +1, 0 , 1);
    m.tail(Jfrac + 1) = pow(1-m.tail(Jfrac+1) , power); 
    m.tail(Jfrac + 1) = 1 - m.tail(Jfrac+1)*X*(1-mid);
    
    return m;
}
