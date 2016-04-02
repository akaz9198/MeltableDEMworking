#ifndef MESH_H
#define MESH_H
//

//

typedef QList <Cbox> oneDarray;
typedef QList <oneDarray> twoDarray;
typedef QList <twoDarray> threeDarray;

class Cbox 
{
	public:

	bool am_I_bottom;
	bool am_I_top;
	QList <Cparticle *> part; 
	QList <Cbox *> contact_box;
};

class Cmesh  
{
public:
	Cvector step;  							/**<Size of each box*/
//	threeDarray box; 						/**<3D array of box*/ 
	class Cbox box[50][50][30];
	int N[DIM];								/**<number of box in each direction, for convinience only*/
	Cmesh(Cvector, double, Ccell);
//	void get_neigbour(Ccell&);
};


#endif
