#include <iostream>
#include <Eigen/Dense>
#include "StructStruct.hpp"
#include "PhysicalConsts.h"
#include "converge.hpp"

using namespace Eigen;

class Saver{ /// Class for saving output values, and printing key values
    public:
        std::string bulkfile , massfile;
        
        void print_accuracies(struct P_P *planet , int &t , double &time , Converger &con ){
            //print out some key values
            std::cout << "M " << planet->M/M_Earth<< " M total " << std::endl;
            if( planet->dt <1e308 ){
                std::cout << "t: " << t << " time: " << time/yrs_to_s << " yrs ";
            }
            
            std::cout << planet->scale.transpose() << " Psurf " << planet->Psurf[0] << " P_bl " << planet->ym(J-1,0) << " Tsurf " << planet->Tsurf[0] <<" T_bl " << planet->ym(J-1,2) ;
            
            std::cout << std::endl;
            
            std::cout << "No. iterations: " << con.n << " of a maximum of: "<< con.nmax << std::endl;
             
            std::cout << "Max innacuracy: " << con.maxinacc << ", Mean Accuracy: " << con.meaninacc << "\n" << std::endl; 

        }

        void save_values(char num , struct P_P *planet , EoS *eos , double &time , double &f_time , Converger &con ){
            /// Save the internal structure
            char timestring[100];

            int a ;
            if(f_time == 0){
                a = snprintf(timestring, 100 , ".txt" );
            }

            std::ofstream fout;
            /// 'Allthing' contains most of the key parameters, plus the errors on P,r,T,L
            /// Fls contains the contribtions to the heat fluxes
            /// 'others' contains some more quantities

            

                Allthing = Eigen::ArrayXXd::Zero(J,I + I + 3);
                Fls = Eigen::ArrayXXd::Zero(J,3);
                
                Allthing.block(0,I+3, J,I) = con.ErrorArray ; 
            

            Allthing.col(0) = planet->m;
            Allthing.block(0,1,J,I) = planet->y;
            Allthing.block(0, I+1,J,1) = planet->rho;
            Allthing.block(0,I+2,J,1) = planet->nabla;
            
            if (num == '0'){
                    std::cout << "no file saved" <<std::endl;
            }else{
                   std::cout << "num " << num << std::endl;
                    fout.open(folder+"/results"+num+basicfile+timestring);
                    std::cout << folder+"/results"+num+basicfile+timestring << std::endl;
                    //fout << "m/M $p/P_C$ r/R rho G1 G2 \n";
                        fout << Allthing ;
                    fout.close();
                    fout.open(folder+"/scale"+num+basicfile+timestring);
                        fout << planet->M <<"\n"<< planet->scale ;
                        
                    fout.close();
                   
                    std::cout << "files saved\n" <<std::endl;
            }
        }
     
        void save_bulk( struct P_P *planet , EoS *eos , int &t , double &time ){
            /// Save various bulk quantities, including the position of melt fronts
            std::ofstream fout;    

            if (t == 0){
                fout.open(massfile );
            }else{
                fout.open(massfile, std::ios::app);
                fout << std::endl;
            }
            fout.precision(8);
            fout << time/yrs_to_s << " " ;
            fout.close();   


            if (t == 0){
                fout.open(bulkfile );
            }
            else{
                fout.open(bulkfile, std::ios::app);
                fout << std::endl;
            }
            fout.precision(8);
            fout << time/yrs_to_s << " " << planet->scale.transpose() << " " << planet->Psurf[0] << " " << planet->Tsurf[0] << " " ;
                
            int J = planet->m.size() , j;            

		        double gamma = 1 + 1.0/planet->n; 
                
                fout << planet->dt/yrs_to_s << " " << t << " " << J;
            fout.close();   
        }
        
        void set_files(struct P_P *planet ,  char num , double &f_time , std::string &_folder ,  std::string EoStype , double &Linit){
            /// Set the names of files for the particular case chosen
            char Tstring[100] , Pstring[100] , Resstring[100] , condstring[100] , corestring[100] , nstring[100] , massstring[100];
            int a;

            folder = _folder;
            std::cout << "EoS type: " << EoStype << std::endl;
            if(EoStype == "Ideal" || EoStype == "Ideal_Tab"){
                a = snprintf(nstring, 100 , "n%.1f_%.3g" , planet->n , planet->mu_m) ;
                a = snprintf(massstring, 100 , "_" );  

            }
            if(f_time == 0){
                a = snprintf(Resstring, 100 , "Res_%.1d" , J );
            }else{
                a = snprintf(Resstring, 100 , "Res_%.1d_ftime_%.2g" , J , f_time);
            }
                
                if(planet->Tbound == 'b'){ 
                    a = snprintf(Tstring, 100 , "L%.2g_", planet->scale(3) ) ;
                }else 
                if(planet->Tbound == 'f'){
                    a = snprintf(Tstring, 100 , "Ts%.4g" , planet->Tsurf[0] );  
                    if(planet->nabtype == 'i' || planet->nabtype == 'a' ){
                        //do nothing
                    }else{
                        a += snprintf(Tstring+a, 100-a , "_L%.2g" , planet->scale(3) );
                    }
                }else{
                    throw std::invalid_argument("No file format found for T boundary" );
                }

                if(planet->Pbound == 'f'){
                    a = snprintf(Pstring, 100 , "Ps%.2g" , planet->Psurf[0] );  
                }
                a = 0;

                    a = snprintf(corestring+a, 100-a , "_" );
                
            if(planet->nabtype == 'a'){
                a = snprintf(condstring, 100 , "AD" );
            }else if(planet->nabtype == 'i'){
                a = snprintf(condstring, 100 , "ISO" );
            }
            
            bulkfile = folder+"/bulk"+num+nstring+massstring+Resstring+corestring+condstring+Tstring+Pstring+".txt";
            massfile = folder+"/massloss"+num+nstring+massstring+Resstring+corestring+condstring+Tstring+Pstring+".txt";
            
            basicfile = std::string(nstring)+massstring+Resstring+corestring+condstring+Tstring+Pstring;
            std::cout << "files set" << std::endl;

        }
        void Not_Converged(struct P_P *planet ,  EoS *eos , Converger &con  ){
            /// If the convergence algorithm has failed, abort the code and save the final values, so they can be used to figure out why
           {
                Allthing = Eigen::ArrayXXd::Zero(J,I + I + 3);
                
                Allthing.block(0,I+3, J,I) = con.ErrorArray ; 
            }
            Allthing.col(0) = planet->m;
            Allthing.block(0,1,J,I) = planet->y;
            Allthing.block(0, I+1,J,1) = planet->rho;
            Allthing.block(0,I+2,J,1) = planet->nabla;

            std::cout << Allthing << std::endl;
            
            int t_broken = 10;
            double time_broken = 0;
            save_bulk( planet , eos , t_broken , time_broken );
            throw std::invalid_argument("Too high an error!" );

        }
        void L_contributions(struct P_P *planet , char *timestring){
            /// Calculates various contributions to the lumminosty (for debugging)
            std::ofstream fout;
            
                Eigen::ArrayXXd Lcontribs(J,5);
                Lcontribs.col(0) = ((planet->y.col(0)*planet->scale(0)-planet->P0)/planet->P0) ;
                Lcontribs.col(1) = ((planet->y.col(2)*planet->scale(2)-planet->T0)/planet->T0);
                Lcontribs.col(2) = planet->dy.col(3);
                Lcontribs.col(3) = (planet->M /planet->scale(3) * planet->dm *  planet->H);
                Lcontribs.col(4) = (planet->M /planet->scale(3) * planet->dm *( - planet->Cp * (planet->scale(2) * planet->y.col(2)-planet->T0)/planet->dt + planet->del_rho * (planet->scale(0)*planet->y.col(0)-planet->P0)/planet->dt ));

            
        }

        Saver(int &_I, int &_J):
            I(_I) , J(_J)
        {       }
        void resize(int &_J){
            J = _J;
        }
    private:
        Eigen::ArrayXXd Fls , Allthing , others , debug_quantities , composition_array;
        std::string basicfile, folder ;
        int I , J ;
        
};