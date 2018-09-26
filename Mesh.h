#ifndef MESH_H_
#define MESH_H_

#include <string>
#include "InitialSystem.h"

using namespace std;

class Mesh {
	public:
            Mesh(unsigned int imax, unsigned int jmax, unsigned int itl, unsigned int itu, InitialSystem* NSC);
            ~Mesh();

            const void print();

            //double geometryCalculation();
            void read_su2(string filename);
            void read_tecplot(string filename);
            void write_tecplot(string filename);
	
	private:
            InitialSystem* initSyst_;
            unsigned int imax_, jmax_, itl_, itu_;
            unsigned int imaxGhost_, jmaxGhost_;
            unsigned int nbKnots_;        	//total number of knots
            unsigned int numberOfCells_;	//total number of cells
            unsigned int rimax_, rjmax_;       /*restricted domain cv imax,jmax*/ // What is this?
            unsigned int inci_, incj_;
            double** x_;
            double** y_;
            double** cellArea_;
            double** normal_i_x_; 			//non normalized
            double** normal_i_y_;
            double** normal_j_x_;
            double** normal_j_y_;

            double** rho_;
            double** u_;
            double** v_;
            double** p_; 							 //primitives variables cell centered
            double** rho_nodes_; 
            double** u_nodes_;
            double** v_nodes_;
            double** p_nodes_;   //vertexes
            double** rho_0_;
            double** u_0_;
            double** v_0_;
            double** p_0_;  				 //intialisation at rk(0)

            double** residualInviscid_rho_;
            double** residualInviscid_u_;
            double** residualInviscid_v_;
            double** residualInviscid_p_;
            double** residualDissip_rho_;
            double** residualDissip_u_;
            double** residualDissip_v_;
            double** residualDissip_p_;
            double** tmp_rho_;
            double** tmp_u_;
            double** tmp_v_;
            double** tmp_p_; //What is this?
            double** deltaT_;
            double** speci_;
            double** specj_;
};

#endif