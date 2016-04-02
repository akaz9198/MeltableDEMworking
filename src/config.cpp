#include "cdinit.cpp"

void Cconfig::iterate(double time_step)
{ 	
	dt=time_step;
	dt_on_2 = dt/2.;
	dt2_on_2=dt*dt/2.;
	
	predictor();			//motion integration
	
	update_contact();		//find contact and get the force and torque
	sum_force(); 			//sum the force moment of each contact on particle   
	if(simule_thermal_conduction) sum_heat(); 					//Heat transfer 
	
	cell.rigid_velocity *= 0.0;
	corrector();			//acceleration of particles according to the sum of force/moment they experience
	cell.rigid_velocity /= parameter.total_mass;

// Velocity offset by rigid motion
#pragma omp parallel for num_threads(NTHREADS)	// YG, MPI
	for(int ip=0; ip< P.size(); ip++) 
		P[ip].V -= cell.rigid_velocity;
}

void Cconfig::predictor()
{ 
	cell.predictor(dt,dt2_on_2);//move the cell

#pragma omp parallel for num_threads(NTHREADS)	// YG, MPI
	for(int  ip=0; ip< P.size();ip++) //move the particles
	{ 	
		P[ip].predictor(dt,dt2_on_2);
		P[ip].set_me_in_main_cell(cell);//set the particle back in the cell if it has got out  
	}
	
//	for(int ig=0;ig<G.size();ig++) G[ig].predictor(dt,dt2_on_2);
	
	
	if(simule_thermal_expansion){
#pragma omp parallel for num_threads(NTHREADS)	// YG, MPI		
		for(int  ip=0; ip< P.size();ip++)  P[ip].expand_radius(dt);}
}    

void Cconfig::corrector()
{ 
cell.corrector(dt_on_2);

#pragma omp parallel for num_threads(NTHREADS)	// YG, MPI
for(int ip=0; ip< P.size();ip++)	
		P[ip].corrector(dt_on_2,cell); 
		
// for(int ig=0;ig<G.size();ig++) G[ig].corrector(dt_on_2,cell);
	 
}

void  Cconfig::sum_heat()
{
	
	double sum_T=0;
	for(int ip=0; ip< P.size();ip++) sum_T += P[ip].T*P[ip].m;
	parameter.average_temperature = sum_T / parameter.total_mass;
	
	PN=0.0; PS=0.0; PT=0.0; PR=0.0;
	for(int ic=0;ic<C.size();ic++)
	{
		PN += C[ic].production_normal;
		PS += C[ic].production_slide;
		PT += C[ic].production_twist;
		PR += C[ic].production_rolling;
	}
	
	heat_in = cell.shear_stress_in * cell.shear_rate + cell.normal_stress_in * cell.dilat_rate;
	heat_in *= cell.L.x[0]* cell.L.x[1]* cell.L.x[2];
	heat_out = 0.0;
	
#pragma omp parallel for num_threads(NTHREADS)	// YG, MPI
	for(int ip=0; ip< P.size();ip++)
	{ 
		P[ip].phi_ext=0;
		P[ip].production=0;
		P[ip].phi =  P[ip].phi_ext; //initialisation for each particles; 0 by default
	}

	for(int ic=0;ic<C.size();ic++) //sum over each interaction
	{
		P[C[ic].A].phi+= C[ic].phi ;
		P[C[ic].B].phi-= C[ic].phi ;
		
		if(C[ic].Flag_Boundary){
			double dPhi = -2.*C[ic].conductivity*C[ic].a *(parameter.average_temperature-20.0)*fabs(C[ic].dX.x[1])/cell.L.x[1]; //exp
			P[C[ic].A].phi += dPhi;
			P[C[ic].B].phi += dPhi;
			heat_out += 2.0*dPhi;
		}

		if(simule_thermal_production)
		{
			P[C[ic].A].production+= C[ic].production/2.; 
			P[C[ic].B].production+= C[ic].production/2.; 		
		}
	}
	
#pragma omp parallel for num_threads(NTHREADS)
	for(int ip=0; ip< P.size();ip++) P[ip].Tdot =  (P[ip].phi +P[ip].production) /(P[ip].c*P[ip].m);
}


void Cconfig::sum_force()
{
 	Cmatrix stress;

	Cvector g;//vector gravity
	if(cell.boundary=="WALL_INCLINED"){
	g.x[0]=  cell.gravity*sin(PI/180 * cell.slope );
	g.x[1]= -cell.gravity*cos(PI/180 * cell.slope );
	g.x[2]=0;
	}//	else g is null
	
#pragma omp parallel for num_threads(NTHREADS)	// YG, MPI
	for(int ip=0; ip< P.size();ip++) 
	{
		P[ip].Fsum = (g*P[ip].m);
		P[ip].Gsum*=0;
	}
	
// YG, no MPI
	for(int ic=0;ic<C.size();ic++)
	{
		Cvector dx;
		P[C[ic].A].Fsum += C[ic].F;
		P[C[ic].B].Fsum -= C[ic].F;
		P[C[ic].A].Gsum +=  (C[ic].RA^C[ic].F) +C[ic].G;
		P[C[ic].B].Gsum -=  (C[ic].RB^C[ic].F) +C[ic].G;

		dx = C[ic].dX;

		stress+=   (C[ic].F| dx);
	}
	

	if(PSEUDO_2D)	stress/=(cell.L.x[0]*cell.L.x[1]);  
	else	stress/=(cell.L.x[0]*cell.L.x[1]*cell.L.x[2]);  

	//cell.stress = stress.symetric();  
	cell.stress =stress;
	
	
	if(cell.boundary!="WALL_INCLINED") return;
	cell.normal_stress_in=0;
	cell.shear_stress_in=0;
	for(int ic=0;ic<C.size();ic++)
		{
			if(P[C[ic].A].AM_I_BOUNDARY==-1||P[C[ic].A].AM_I_BOUNDARY==-2){ cell.normal_stress_in+=C[ic].F.x[1]; cell.shear_stress_in+=C[ic].F.x[0];}
			if(P[C[ic].B].AM_I_BOUNDARY==-1||P[C[ic].B].AM_I_BOUNDARY==-2){ cell.normal_stress_in-=C[ic].F.x[1]; cell.shear_stress_in-=C[ic].F.x[0];}
		}
	cell.normal_stress_in/=(cell.L.x[0]*	cell.L.x[2]);cell.shear_stress_in/=(cell.L.x[0]*	cell.L.x[2]);
	cell.normal_stress_in*=-1;//get it positive in compression
}

void Cconfig::update_contact()
{
	renew_contact_list();

#pragma omp parallel for num_threads(NTHREADS)	// YG, MPI 
	for(int ic=0;ic<C.size();ic++)
	{
		C[ic].EVALE_Geo();
		C[ic].relative_velocity();
		C[ic].increment_force(dt);
	}
	
	if(simule_thermal_conduction){
#pragma omp parallel for num_threads(NTHREADS)	// YG, MPI	
	for(int ic=0;ic<C.size();ic++)	C[ic].EVALE_heat_flow();
	}
}


void Cconfig::renew_contact_list()//Find contact via the mesh method
{
// SEG FAULT: destructor of old mesh?

	//FIND NEIGHBOURS BY CELL METHOD
	Cmesh mesh(cell.L, 1.2*parameter.Dmax,cell); 
	
// YG, no MPI here
	for( int ip=0;ip<P.size();ip++)	P[ip].set_in_box(mesh,cell); //set particles in box
		
#pragma omp parallel for num_threads(NTHREADS)	// YG, MPI
	for( int ip=0;ip<P.size();ip++)	P[ip].get_neighbour(cell);//guess what
		
// YG, Be careful for the following section		
//	for(int it=0; it<NTHREADS; it++){ CThread[it].clear();} // done below
		
	int in, ip, c, tid; 
	bool exist_contact;
	Ccontact cont;
#pragma omp parallel for private(in, c, exist_contact,cont, tid) schedule(dynamic) num_threads(NTHREADS)	// YG, MPI testing
	for(ip=0;ip<P.size();ip++)	// Problem with this loop!! Inverse the loop, different results!!
		{
		for(in=0; in<P[ip].neighbour.size();in++)
			{
				exist_contact=false;
				for(c=0; c<P[ip].contact.size();c++)
					if(P[ip].contact[c]->pB == P[ip].neighbour[in] || P[ip].contact[c]->pA == P[ip].neighbour[in] )
						{	exist_contact= true; break; }//check if the contact exists
				if(!exist_contact)//if not, create a new contact
				{ 	
					tid = omp_get_thread_num();
					Ccontact cont(&P[ip], P[ip].neighbour[in], &cell, &parameter);	
					if( cont.AM_I_CONTACTING() )
					{
//#pragma omp critical 
//						{C<<cont;}
						CThread[tid]<<cont; //don't do it if particles are not contacting 
					}													
				}
			}
		}
	for(int it=0; it<NTHREADS; it++){C+=CThread[it]; CThread[it].clear();}
	
// YG, no MPI here	
	for(int ic=0; ic < C.size(); ic++) //delete contact within the list
		if (!C[ic].AM_I_CONTACTING()) 
			{ 
				if(ic< C.size()-1 ) 
				{
				C[ic]=C[C.size()-1]; 
				ic--; 
				}
			C.pop_back(); 
			}

#pragma omp parallel for num_threads(NTHREADS)	// YG, MPI   
	for(int ic=0; ic < C.size(); ic++)  C[ic].age+=dt;//get contact older
#pragma omp parallel for num_threads(NTHREADS)	// YG, MPI		
	for(int ip=0;ip<P.size();ip++) P[ip].contact.clear();//rebuilt the pointer list
// YG, no MPI here
    for(int ic=0; ic < C.size(); ic++) {
    	C[ic].pA->contact.push_back(&C[ic]);
//    	C[ic].pB->contact.push_back(&C[ic]); // YG, added
		}
		
}




 void Cconfig::energy(void)
{ 
 // double Etot,Eela,Eforce,Emoment, Ekin,Erot,Etrans;
   Etot=0;
   Eela=0;Eforce=0;Emoment=0;
   Ekin=0;Erot=0; Etrans=0;
   
 /* foreach(Ccontact cont, C)
  {
//    Eforce += 2./5. *cont.fn*cont.deltaN + 0.5 * cont.ft*cont.ft/(cont.a*parameter.MODULE_T);
    
   cont.gt=cont.Gt.NORM(); 
   cont.gn=con.Gn.NORM();

  //  Emoment+= 0.5/(cont.a*cont.a*cont.a)*(cont.gn*cont.gn/parameter.MODULE_N +  cont.gt*cont.gt/parameter.MODULE_T);
  }
   Eela=Eforce+Emoment;
  
  foreach(Cparticle p, P)
  {
    Cvector V;
    V = p.V;
    V.x[0]-=cell.shear_rate*p.X.x[1];
    Etrans+=0.5*p.m*(V*V);
    Erot +=0.5* 2./5.*p.m* p.R*p.R*(p.Ome*p.Ome);
 
  }
  */
  Ekin=Erot+Etrans;
  Etot = Ekin+Eela;
}



void Cconfig::fread(Cin_out where_to_read)
{
	where_to_read.set_file_name();
	where_to_read.check_file();

    P.clear(); C.clear();

	ifstream file;

	file.open(where_to_read.save_file[0].c_str());  // read particles
	while(1<2)
	{
		Cparticle part;
		file>>part;
		if(!file.eof()) P.push_back(part);
		else break;
	}
	file.close();
	
	file.open(where_to_read.save_file[1].c_str());  // read contact
	while(1<2)
	{
		Ccontact cont;
		file>> cont;
		if(!file.eof()) C.push_back(cont);
		else break;
	}
	file.close();


	file.open(where_to_read.save_file[2].c_str());//save parameter
	file>>t;
	file>>parameter; 
	file.close();

	file.open(where_to_read.save_file[3].c_str()); //save cell
	file>>cell; 
	file.close();
	
	
	//bluild the pointers
	for(int ip=0;ip<P.size();ip++)//cell knows which particle is a plan
	{
		if(P[ip].AM_I_BOUNDARY==-2 ) cell.plan_bottom=&P[ip];
		if(P[ip].AM_I_BOUNDARY==2 )  cell.plan_top=&P[ip];
		P[ip].id=ip;
	}

	for(int ic=0;ic<C.size();ic++)
		{ 
			C[ic].pA = &P[C[ic].A]; 
			C[ic].pB = &P[C[ic].B];
			C[ic].cell = &cell;
			C[ic].parameter = &parameter;
			P[C[ic].A].contact.push_back(&C[ic]);
		} 

	update_particle();
	parameter.dimensionless_number(cell,P);	
	iterate(0.0);//use here to recover all the data that haven't been saved, but which derive from saved data

}

void Cconfig::fprint(Cin_out where_to_save)
{
	ofstream file;
	
	where_to_save.set_file_name();
	where_to_save.check_file();

	file.open(where_to_save.save_file[0].c_str());  // save particles
	foreach(Cparticle part, P) file<<part;
	file.close();

	file.open(where_to_save.save_file[1].c_str()); // save contact
	foreach(Ccontact cont, C) file<<cont;
	file.close();

	file.open(where_to_save.save_file[2].c_str()); //save parameter
	file<<t<<"\t";
	file<<parameter;
	file.close();
	
	file.open(where_to_save.save_file[3].c_str()); //save cell
	file<<cell; 
	file.close();

	//YG
//	file.open(where_to_save.save_file[4].c_str(), std::ios_base::app); 
//	file<<t<<"\t"<<cell.Xshift<<"\t"<<cell.L.x[1]<<"\t"<<cell.stress.x[1][1]<<"\t"<<cell.stress.x[1][0]<<"\t"
//	<<parameter.average_temperature<<"\t"<<heat_in<<"\t"<<heat_out<<"\t"<<PN<<"\t"<<PS<<"\t"<<PT<<"\t"<<PR<<endl; 
//	file.close();
	
//int who= RUSAGE_SELF;
//struct rusage usage;
//struct rusage *puse=&usage;
//getrusage(who,puse);
	
//	file.open(where_to_save.save_file[5].c_str(), std::ios_base::app); 
//	file<<t<<"\t"<<puse->ru_maxrss<<"\t"<<puse->ru_ixrss<<"\t"<<puse->ru_idrss<<"\t"<<puse->ru_isrss<<"\t"
//	<<puse->ru_minflt<<"\t"<<puse->ru_majflt<<endl;
}


void  Cconfig::Evale_conductivity_tensor() /**< Measure the conductivity tensor */
{
	Cmatrix K,H; //effective conductivity, convectivity

	K*=0.;
/*	foreach(Ccontact *cont, C)    K+= cont->alpha*(2.*cont->a*(cont->dX*cont->dX)) ;
	if(PSEUDO_2D)cell.K =K/(cell.L.x[0]*cell.L.x[1]);
	else cell.K = K/(cell.L.x[0]*cell.L.x[1]*cell.L.x[2]);
*/
	H*=0.;
	foreach(Cparticle part, P) 
	{ Cvector deltaV;
	deltaV=part.V;
	deltaV.x[0]-= cell.shear_rate*part.X.x[1];
	H+= ( deltaV|part.X*part.c*part.m); //*
	}

	if(PSEUDO_2D) cell.H =H/(cell.L.x[0]*cell.L.x[1]);
	else cell.H = H/(cell.L.x[0]*cell.L.x[1]*cell.L.x[2]);

	cell.production=0;
/*	foreach(Ccontact *cont, C) cell.production+=cont->production;
	if(PSEUDO_2D)  cell.production/=(cell.L.x[0]*cell.L.x[1]);
	else cell.production/=(cell.L.x[0]*cell.L.x[1]*cell.L.x[2]);
*/
}




