//======================================
// Class contact data array and functions
// Functions details are in "aggregation.c"
//======================================

//YG
class Caggregate
{
	public:
	int id;
	Cvector X;             /**< Mass Center position. */ 
	Cvector V;             /**< Aggregation velocity.*/ 
	Cvector A;             /**< Aggregation acceleration. */ 
	Cvector Ome;           /**< Aggregation angular velocity.*/
	Cvector OmeDot;        /**< Aggregation angular acceleration. */ 
	Cvector Fsum;          /**< Sum of forces the aggregation is subjected to.*/
	Cvector Gsum;          /**< Sum of moment the aggregation is subjected to.*/
	
	double m;              /**< Aggregation mass. */
	Cmatrix J;              /**< Aggregation moment of inertia. Tensor */
	
	QList <Cparticle *> Pa;	/**<List of particle as member of this aggregation*/
	QList <Cvector> dX;
//	QList <Ccontact *> contact;	/**<List of contact form this aggregation*/
	
	Ccell *cell;
	
	// Functions
	void predictor(double dt, double dt2_on_2);	/**< Integrate the velocity and the acceleration to get the new postion. \warning Depend on wether there is inertia or not.*/
	void corrector(double dt_on_2,Ccell &); 	/**< Get the new velocity and acceleration \warning Only with inertia. */	

	Caggregate();
	void empty();
	
	void merge(Cparticle &);	
	void merge(Caggregate &);
	void breakage();
	void update(Ccell &c);
	
	friend ofstream &operator<<(ofstream &,Caggregate);		/**< Print the particle's data into a file.*/ 
	friend ifstream & operator>>(ifstream &,Caggregate &);	/**< Read the particle's data from a file.  */ 
};