#include "Halide.h"
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h> 
using namespace Halide;

//Func buildRbfFunc(std::string call, Func rij, Func scalefunc, Expr rin, Expr rmax, int pdegree, int K, )

 void buildSnap(Func &  rbf, Func xij, Func besselparams, Expr rin, Expr rcut,
	       Expr besseldegree, Expr inversedegree, Expr nbseelpars, Expr npairs)
{
  Func abf("snap_abf");
  Func rbfp("snap_rbfr");

  Var np("nij"), bfp("snap_bfp"), bfd("snap_bfd"), ibfp("snap_abfp");

  Expr one = Expr((double) 1.0);
  Expr rmax = rcut - rin;
  Expr dij = xij(np);
  Expr r = dij - rin;        
  Expr y = r/rmax;      
  Expr y2 = y*y;
  Expr y3 = one - y2*y;
  Expr y4 = y3*y3 + Expr((double)1e-6);
  Expr y5 = sqrt(y4);
  Expr y6 = exp(-one/y5);
  Expr fcut = y6/exp(-one);

  Expr alpha = max(besselparams(bfp), Expr((double) 1e-3));
  Expr x =  (one - exp(-alpha*r/rmax))/(one-exp(-alpha));


  Expr a = (bfd+1)*Expr((double)M_PI);
  Expr b = (sqrt(2 * one/(rmax))/(bfd+1));

  rbfp(bfp, bfd, np) =  b*fcut*sin(a*x)/r;
  rbfp.compute_root();
  abf(ibfp, np) = fcut/pow(dij, ibfp + 1);
  abf.compute_root();

  rbf(bfp, bfd, np)= Expr((double) 0.0);
  RDom rpp(0, npairs, 0, besseldegree, 0, nbseelpars);
  rbf(rpp.z, rpp.y, rpp.x) = rbfp(rpp.z, rpp.y, rpp.x);

  RDom rp(0, npairs, 0, inversedegree);
  rbf(nbseelpars, rp.y, rp.x) = abf(rp.y, rp.x);
  //  return rbf;
}

void buildRBF(Func & rbf, Func & drbf, Func & abf, Func & dabf,
	      Func xij, Func besselparams, Expr rin, Expr rmax,
	      Expr bdegree, Expr adegree, Expr nbparams, Expr npairs,
	      Var bfi, Var bfp, Var np, Var dim)
{

  Expr one = Expr((double) 1.0);
  Expr zero = Expr((double) 0.0);
  Expr xij1 = rbf(np, 0);
  Expr xij2 = rbf(np, 1);
  Expr xij3 = rbf(np, 2);

  Expr s = xij1*xij1 + xij2*xij2 + xij3*xij3;
  Expr dij = sqrt(s);
  Expr dr1 = xij1/dij;    
  Expr dr2 = xij2/dij;    
  Expr dr3 = xij3/dij;    

  Expr r = dij - rin;        
  Expr y = r/rmax;    
  Expr y2 = y*y;
  Expr y3 = one - y2*y;
  Expr y4 = y3*y3 + Expr((double) 1e-6);
  Expr y5 = sqrt(y4); //pow(y4, 0.5);
  Expr y6 = exp(-one/y5);
  Expr y7 = pow(y4, 1.5);
  Expr fcut = y6/exp(-one);
  Expr dfcut = ((3.0/(rmax*exp(-one)))*(y2)*y6*(y*y2 - one))/y7;

  Expr alpha = min(Expr((double)1e-3), besselparams(bfp));
  Expr x =  (one - exp(-alpha*r/rmax))/(one-exp(-alpha));
  Expr dx = (alpha/rmax)*exp(-(alpha*r/rmax))/(one - exp(-alpha));

  Expr a = (bfd + 1) * M_PI;
  Expr b = sqrt(2 * one/rmax)/(bdf + 1);

  rbf(bfp, bfd, np) = b * fcut * sin(a*x)/r;
  rbf.bound(bfp, 0, nbparams);
  rbf.bound(bfd, 0, bdegree);
  rbf.bound(np, 0, npairs);
  Expr drbfdr = b*(dfcut*sin(a*x)/r - fcut*sin(a*x)/(r*r) + a*cos(a*x)*fcut*dx/r);
  drbf(bfp, bfi, np, dim) = (rbf(np, dim)/dij) * drbfdr;
  drbf.bound(dim, 0, 3);
  drbf.bound(bfp, 0, nbparams);
  drbf.bound(bfi, 0, bdegree);
  drbf.bound(np, 0, npairs);

  abf(bfi, np) = fcut/pow(dij, bfi+one);
  abf.bound(bfi, 0, adegree);
  abf.bound(np, 0, npairs);
  Expr drbfdr_a = dfcut/a - (bfi+1.0)*fcut/(a*dij);
  dabf(bfi, np, dim) = (rbf(np, dim)/dij) * drbfdr_a;
  dabf.bound(dim, 0, 3);
  dabf.bound(bfi, 0, adegree);
  dabf.bound(np, 0, npairs);
}

void buildStructureMatMul(Func & energyij, Func & forceij,
			  Func rbf, Func abf, Func drbf, Func dabf, Func Phi1, Func Phi2,
			  Expr bdegree, Expr adegree, Expr tdegree, Expr nbparams, Expr npairs,
			  Var bfi, Var bfa, Var bfp, Var np, Var dim){
  //Multiply atom * basisfunction  by basisfunction * rbf
  Expr zero = Expr((double) 0.0);

  energy(bfi, np)= zero;
  fenergy(bfi, np, dim) = zero;
  
  energy.bound(bfi, 0, tdegree);
  energy.bound(np, 0, npairs);
  fenergy.bound(bfi, 0, tdegree);
  fenergy.bound(np, 0, npairs);
  fenergy.bound(dim, 0, 3);

  RDom rbf(0, bdegree, 0, nbparams);
  //  RDom drbf(0, bdegree, 0, nbparams, 0, 3);
  energy(bfi, np) += rbf(rbf.y, rbf.x, np) * Phi1(rbf.y, rbf.x, bfi);//ordering here is questionable
  fenergy(bfi,np, dim) += drbf(rbf.y, rbf.x, np, dim) * Phi1(rbf.y, rbf.x, bfi);//ordering here is questionable
  energy.updates(0).bound(bfi, 0, bdegree);
  fenergy.updates(0).bound(bfi, 0, bdegree);

  RDom abf(0, adegree);
  energy(bfa, np) += abf(abf.x, np) * Phi2(abf.x, bfa);//ordering here is questionable
  fenergy(bfa, np, dim) += abf(abf.x, np, dim) * Phi2(abf.x, bfa);//ordering here is questionable
  energy.updates(1).bound(bfa, 0, bdegree);//check this
  fenergy.updates(1).bound(bfa, 0, bdegree);//check this
  
}



void buildPodTally2b(Func & eatom, Func & fatom,
		     Func eij, Func fij, Func ai, Func aj, Func ti, Func tj, Func interaction,
		     Expr npairs, Expr natom, Expr nelements, Expr nelementCombos, Expr nbf,
		     Var np, Var atom, Var bf, Var dim, Var elem, var inter
		     )
{
  Expr one = Expr((double) 1.0);
  Expr zero = Expr((double) 0.0);
  eatom(bf, atom, np) = zero;
  fatom(bf, atom, np, dim) = zero;

  eatom.bound(atom, 0, natom);
  eatom.bound(np, 0, npairs);
  eatom.bound(bf, 0, nbf);
  
  fatom.bound(atom, 0, natom);
  fatom.bound(bf, 0, nbf);
  fatom.bound(np, 0, npairs);
  fatom.bound(dim, 0, 3);

  RDom r(0, npairs, 0, nbf);

  Expr i1 = clamp(ai(r.x), 0, natom);
  Expr j1 = clamp(aj(r.x), 0, natom);
  Expr typei = clamp(ti(r.x) - 1, 0, nelements);
  Expr typej = clamp(tj(r.x) - 1, 0, nelements);

  Expr inter_ij = clamp(interaction(typei, typej) - 1, 0, nelementCombos);
  eatom(r.y, inter_ij, i1) += eij(r.y, r.x);
  fatom(r.y, inter_ij, i1, dim) += fij(r.y, r.x, dim);
  fatom(r.y, inter_ij, j1, dim) -= fij(r.y, r.x, dim);

}

void buildNeighPairs(Func & outputs, Func & vectors,
		     Func pairlist, Func pairnumsum, Func atomtype, Func alist, Func atompos,
		     Expr natom, Expr dim, Expr nmax, Expr npairs,
		     Var atom, Var d, Var nm, Var np, Var numOuts){
  
  outputs(np, numOuts) = mux(numOuts,{-1, -1, -1, -1});
  outputs.bound(np, 0, npairs);
  outputs.bound(numOuts, 0, 4);
  
  vectors(np, d) = Expr((double) 0.0);
  vectors.bound(d, 0, dim);
  vectors.bound(np, 0, npairs);


  RDom r(0, natom, 0, nmax);
  r.where(r.y < pairnumsum(r.x + 1) && r.y >= pairnumsum(r.x));

  Expr jacc = clamp(pairlist(r.y), 0, npairs);
  Expr att = clamp(alist(jacc), 0, natom); 
  Expr att_tt = atomtype(att); 
  outputs(r.x, numOuts) = mux(numOuts, {r.x, att, atomtype(r.x), att_tt});
  vectors(r.x, d) = atompos(r.y, d) - atompos(r.x, d);
  
}


void buildPod1Body(Func & eatom, Func atomtype,
		   Expr nelements, Expr natom,
		   Var i, Var m){
  eatom(i, m) = select(atomtype(i) == m, (Expr((double)1.0)), (Expr((double)0.0)));
  eatom.bound(i, 0, natom);
  eatom.bound(m, 0, nelements);
}

void buildPod1Body_p(Func & eatom, Func & fatom,
		     Func atomtype,
		     Expr nelements, Expr natom, Expr dim,
		     Var i, Var m){
  eatom(i, m) = select(atomtype(i) == m, (Expr((double)1.0)), (Expr((double)0.0)));
  eatom.bound(i, 0, natom);
  eatom.bound(m, 0, nelements);
  fatom(i) = (Expr((double)0.0));
  fatom.bound(i, 0, dim * natom * nelements);
}



class pod1 : public Halide::Generator<pod1> {
public:
  //Func pairnumsum, Func pairlist,
  //				       Expr NPairs, Expr NAtoms, Expr NMax, Expr dim,
  //				       Func atomtype, Func alist, Func atompos
  
  Input<Buffer<int>> pairlist{"pairlist", 1};
  Input<Buffer<int>> pairnumsum{"pairnumsum", 1};
  Input<Buffer<int>> atomtype{"atomptype", 1};
  Input<Buffer<int>> alist{"alist", 1};
  Input<Buffer<double>> atompos{"atompos", 2};

  Output<Buffer<double>> rij{"rij", 2};
  Output<Buffer<int>> meta{"meta", 2};

  Pipeline pipeline;

  GeneratorParam<int> NMax{"NMax", 100};
  GeneratorParam<int> NTypes{"NTypes", 3};
  //  GeneratorParam<int> M{"M", 50};
  //  GeneratorParam<int> I{"I", 5000};
  


  void generate (){

    Expr NPairs = pairlist.dim(0).max();
    Expr NAtoms = atomtype.dim(0).max();
    alist.dim(0).set_bounds(0, NAtoms);
    atompos.dim(0).set_bounds(0, NAtoms);
    atompos.dim(1).set_bounds(0, 3);
    Func np1_data, np1_vecs;

    Var i("i"), j("j");

    //    rij(i, j)= np1_vecs(i, j);
    //    meta(i, j) = np1_data(i, j);

  }

};


class snapshot : public Halide::Generator<snapshot> {
public:
  //Func pairnumsum, Func pairlist,
  //				       Expr NPairs, Expr NAtoms, Expr NMax, Expr dim,
  //				       Func atomtype, Func alist, Func atompos
  
  Input<Buffer<double>> xij{"xij", 1};
  Input<Buffer<double>> besselparams{"besselparams_buf", 1};
  Input<double> rin{"rin", 1};
  Input<double> rcut{"rcut", 1};
  Input<int> inversedegree_pre{"inversedegree", 1};
  Input<int> bessel_degree_pre{"radialdegree", 1};
  Input<int> nbssselpars{"besselparams", 1};
  Input<int> npairs{"npairs", 1};
  Output<Buffer<double>> rbf{"rbf", 3};


  void generate (){
    xij.dim(0).set_bounds(0, npairs);
    besselparams.dim(0).set_bounds(0, nbssselpars);
    //    Expr npairs = xij.dim(0).max();
    //    Expr nbseelparams = besselparams.dim(0).max();
    //    Expr totbesseldegree = rbf.dim(1).max();
    //    Expr nbesselpars = rbf.dim(2).max();
    //    Expr inversedegree = max(min(totbesseldegree, inversedegree_pre), 0);
    //    Expr radialdegree = max(min(totbesseldegree - inversedegree, totbesseldegree), 0);


    Func temp;
    buildSnap(temp, xij, besselparams, rin, rcut, bessel_degree_pre, inversedegree_pre, nbssselpars, npairs);
    Var a, b,c;
    rbf(a,b,c) = temp(a,b,c);
    rbf.dim(2).set_bounds(0, npairs);
    rbf.dim(1).set_bounds(0, max(inversedegree_pre, bessel_degree_pre));
    rbf.dim(0).set_bounds(0, nbssselpars + 1);
    
  }

};
  

HALIDE_REGISTER_GENERATOR(pod1, pod1);
HALIDE_REGISTER_GENERATOR(snapshot, snapshot);

// td::tuple<Func, Func> buildNeighPairs(std::string call, Func pairnumsum, Func pairlist,
// 				       Expr NPairs, Expr NAtoms, Expr NMax, Expr dim, Expr NTypes,
// 				       Func atomtype, Func alist, Func atompos
// 		     ){



//   Func outputs(call + "_bnp");
//   Func vectors(call + "_bnpd");

//   Var ii("ii_" + call), c("c_" + call), j("j_" + call);
//   //initialize ouputs
//   outputs(j, c) = mux(c, {-1, -1, -1, -1});//ai, aj, ti, tj
//   outputs.bound(c, 0, 4);
//   outputs.bound(j, 0, NPairs);

//   vectors(j, c) = Expr((double) 0.0);
//   vectors.bound(c, 0, dim);
//   vectors.bound(j, 0, NPairs);


//   RDom r(0, NAtoms, 0, NMax);
//   r.where(r.y < pairnumsum(r.x + 1) && r.y >= pairnumsum(r.x)); //iteration

//   Expr jacc = clamp(pairlist(r.y), 0, NPairs);
//   Expr att = clamp(alist(jacc), 0, NAtoms); 
//   Expr att_tt = atomtype(att); 
//   outputs(r.x, c) = mux(c, {r.x, att, atomtype(r.x), att_tt});

//   vectors(r.x, c) = atompos(r.y, c) - atompos(r.x, c);

//   return std::make_tuple(outputs, vectors);
  
// }
