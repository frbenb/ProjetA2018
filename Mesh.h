#ifndef MESH_H_
#define MESH_H_

#include <string>
#include "math.h"
#include "InitialSystem.h"

using namespace std;

class Mesh {
	public:
            Mesh();
            ~Mesh();

            //double geometryCalculation();
            void read_su2(string filename);
            void read_tecplot(string filename);
            void write_tecplot(string filename);

            void iterate_pseudo_timestep(int level, int nstage);

            void timestep();
                
            void iteratestep();

            void save_w0();

            void spectral_radius(int level);

            void residual(int level, double beta, int istage,int dissip);

            void eflux();

            void dflux(int level, int beta);

            void dflux2(int level, int beta);

            void update_solution();

            void update_boundary();
	
	private:
            InitialSystem* NSC_;
            unsigned int imax_, jmax_;          //imax, jmax
            unsigned int imaxGhost_, jmaxGhost_;//himax, hjmax
            unsigned int nbKnots_;        	    //nbpt??? - total number of knots
            unsigned int numberOfCells_;	    //total number of cells
            unsigned int rimax_, rjmax_;        //rimax, rjmax - restricted domain cv imax,jmax
            unsigned int inci_, incj_;          //inci, incj - address  increments in i,j
            
            double** x_;                //x, y - mesh coordinates
            double** y_;
            double** cellArea_;         //area
            double** normal_i_x_; 		//six, siy - face i projections
            double** normal_i_y_;
            double** normal_j_x_;       //sjx, sjy - face j projections
            double** normal_j_y_;

            double** rho_;              //ro, uu, vv, pp - primitive variables cell centered
            double** u_;
            double** v_;
            double** p_;
            double** rho_nodes_;        //rocv, uucv, vvcv, ppcv - primitive variables cell-vertex
            double** u_nodes_;
            double** v_nodes_;
            double** p_nodes_;
            double** rho_0_;            //ro0, ru0, rv0, re0 - conservative variables rk(0)
            double** u_0_;
            double** v_0_;
            double** p_0_;

            double** residualInviscid_rho_; //Ri_ro, Ri_uu, Ri_vv, Ri_pp - inviscid Residual
            double** residualInviscid_u_;
            double** residualInviscid_v_;
            double** residualInviscid_p_;
            double** residualDissip_rho_;   //???
            double** residualDissip_u_;
            double** residualDissip_v_;
            double** residualDissip_p_;
            double** tmp_rho_;              //tmp_ro, tmp_uu, tmp_vv, tmp_pp - artificial Residual
            double** tmp_u_;
            double** tmp_v_;
            double** tmp_p_;
            double** deltaT_;               //dt
            double** speci_;                //speci, specj
            double** specj_;
};

#endif


