#include <iostream>
#include <Eigen/Dense>
#include "StructStruct.hpp"
#include "PhysicalConsts.h"
#include <boost/math/interpolators/pchip.hpp>
#include <vector> 
#include "Print_and_save.hpp"
#include "Interpolator.hpp"


using namespace Eigen;
class stepper{
    public:
        Saver print;
        
        void timestepper(struct P_P *planet , EoS *eos , int tmax , double yrs_max , double acc , char *num , std::string &_folder){
            /** Step forward in time for ascribed steps / time limit  */ 
            // This is just the bare bones of what it should have
        
            double time0,dt0;
            t=0;
            time = 0;
            
            for(t=1; t<=tmax ; t++){
		        time0 = time;
                dt0 = planet->dt;
                
                time = time0 + planet->dt;
                
                con.converge(planet , eos , acc);                
                  
                print.print_accuracies(planet , t , time , con );
                //std::cout << planet->y.col(3).transpose() << std::endl;
                if( *num != '0'){
                     print.save_bulk(planet , eos , t , time );
                }
                if(con.maxinacc > acc){
                    std::cout << "\nHigh error!" <<std::endl;
                    
             	}
                if (t%1000 == 0 || t == tmax ){
                    print.save_values(*num , planet , eos  , time , f_time , con);
                } 
                if (time/yrs_to_s >= yrs_max || planet->M/M_Earth < 1e-7){
                    print.save_values(*num , planet , eos  , time , f_time , con);
                    break;
                }
            }

        }

        
        void Merge_Regrid(struct P_P *planet , struct  EoS *eos ){
            /** Regrid by merging and splitting cells */

            
        }

        void static_struct(struct P_P *planet , EoS *eos ,  double acc){
            /** Compute the structure assuming no time evolution (t-dependent terms -> 0) */

            planet->dt = 1e308;
		    f_time = 0;
            planet->T0 = planet->scale(2) * planet->y.col(2);
            planet->P0 = planet->scale(0) * planet->y.col(0);

            con.converge(planet , eos , acc);
            print.print_accuracies(planet , t , time , con );
            

            if(con.maxinacc > acc){
                print.save_values('1' , planet, eos , time , f_time , con);
                print.Not_Converged(planet , eos , con);
            }
        }
        
        void step_chooser(struct P_P *planet , double f ){
            
        }

    stepper(int _I, int _J):
        I(_I) , J(_J) , con(_I,_J) , print(_I,_J)
        {};
    
        Converger con;

    private:
        int I , J , t;
        double time , f_time , Linit;
        
       
};
