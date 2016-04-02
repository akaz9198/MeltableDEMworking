// place your code here
Caggregate::Caggregate(){}

void Caggregate::merge(Cparticle &p)
{
	p.member = 1;
	p.G = this;
	Pa.push_back(&p);
}

void Caggregate::merge(Caggregate &g)
{
	for(int ip=0;ip<g.Pa.size();ip++) 
	{
		g.Pa[ip]->G = this;
		Pa.push_back(g.Pa[ip]);	
	}
	g.empty();
}

void Caggregate::breakage()
{
	double Fcritical =1.0;
	for(int ip=0;ip<Pa.size();ip++){
		double F = (Pa[ip]->Fsum*dX[ip]) /dX[ip].NORM();
		if(F>Fcritical)
		{
			Pa[ip]->member = 0;
			Pa[ip]=Pa[Pa.size()-1];
			dX[ip]=dX[Pa.size()-1];
			ip--; 
			Pa.pop_back();
			dX.pop_back();
		}
	if(Pa.size()==1){
		Pa[0]->member = 0;
		Pa.pop_back();
		dX.pop_back();
	}
	}
}

void Caggregate::update(Ccell &c) // update mass center, etc.
{
	cell = &c;
	
	dX.clear();
	m = 0.0; J*=0.0; X*=0.0;
	Cvector dx;
	
	for(int ip=0;ip<Pa.size();ip++)
	{ 
//		Pa[ip]->G = this;
		m += Pa[ip]->m;
		X += Pa[ip]->X *Pa[ip]->m;
	}
	X /= m;							// Center of mass
	for(int ip=0;ip<Pa.size();ip++)
	{
		dx = Pa[ip]->X -X;
		cell->rescale(dx);
		dX<<dx;						// position of particle inside aggregation	
		J += MI*Pa[ip]->J + (MI*dx.NORM()*dx.NORM() - (dx|dx))*Pa[ip]->m; // moment of inertia	
	}
}

void Caggregate::predictor(double dt, double dt2_on_2)
{
	Cvector deltax, dome;
	deltax = (V*dt)+ (A*dt2_on_2);
	dome = (OmeDot*dt);
	X += deltax;
	V += (A*dt);
	Ome+= dome;
// To Do
	for(int ip=0;ip<Pa.size();ip++)
	{
		Pa[ip]->X += deltax + (dome^dX[ip]);
		Pa[ip]->V = V + (OmeDot^dX[ip]);
		Pa[ip]->Ome += dome;
		
		if(id==0) Pa[ip]->Tdot += 5.0;		//exp
		Pa[ip]->T+= Pa[ip]->Tdot*dt;
	}	
}

void Caggregate::corrector(double dt_on_2,Ccell &cell)
{
	Cvector Ap,OmeDotp;
	Ap=A;
	OmeDotp=OmeDot;
	
	A = Fsum/m;
	V += (A-Ap)*dt_on_2;
//	OmeDot *= 0.0;
	OmeDot= Gsum*J.INVERSE();
	Ome += (OmeDot-OmeDotp)*dt_on_2;
// To Do	
	for(int ip=0;ip<Pa.size();ip++)
	{
		Pa[ip]->OmeDot = OmeDot;
		Pa[ip]->Ome = Ome;
		Pa[ip]->A = A;		
		Pa[ip]->V = V+ (OmeDot^dX[ip]);
	}
	
}

ofstream & operator<<(ofstream &file,Caggregate g)
{
	for(int ip=0;ip<g.Pa.size();ip++) file<<g.Pa[ip]->id<<"\t";
//	file<<g.J.x[1][1]<<"\t"<<g.J.x[1][2]<<"\t"<<g.Ome.x[1]<<"\t"<<g.Pa[0]->J<<"\t"<<g.m;
	file<<endl;	//new line
 	return file;
}

ifstream & operator>>(ifstream &file,Caggregate &g)
{
 	return file;
}

void Caggregate::empty()
{
	Pa.clear();
	dX.clear();
}