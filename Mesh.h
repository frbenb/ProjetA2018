#ifndef MESH_H_
#define MESH_H_

#include <vector>
#include <iostream>

using namespace std;

class Mesh {

	public :
	void Mesh();
	~ Mesh();

	double geometryCalculation());
	void read_su2(string filename);
	void read_tecplot(string filename);
	void write_tecplot(string filename);

	
	private :
	unsigned int imax_;
	unsigned int jmax_;
	unsigned int imaxGhost_;
	unsigned int jmaxGhost_;
	unsigned int nbKnots_;        	//total number of knots
	unsigned int numberOfCells_;	//total number of cells
	unsigned int inci_;
	unsigned int incj_;
	double** x_;
	double** y_;
	double** cellsAera_;
	double** normal_i_x_; 			//non normalized
	double** normal_i_y_;
	double** normal_j_x_;
	double** normal_j_y_;

	double** rho_, u_, v_, p_; 							 //primitives variables cell centered
	double** rho_nodes_, u_nodes_, v_nodes_, p_nodes_;   //cell vertexed
	double** rho_0_, u_0_, v_0_, p_0_;  				 //intialisation at rk(0)

	double** residualInviscid_rho_, residualInviscid_u_, residualInviscid_v_, residualInviscid_p_;
	double** residualDissip_rho_, residualDissip_u_, residualDissip_v_, residualDissip_p_;
	double** tmp_rho_, tmp_u, tmp_v, tmp_p; //on sait pas cest quoi
	double** deltaT_;
	double** speci_, specj_;

}
